#include "rosic_WaveTable.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveTable::WaveTable()
{
  startPhase                  = 0.0;
  rightPhaseOffset            = 0.0;
  fullWaveWarp                = 0.0;
  halfWaveWarp                = 0.0;
  rangeMin                    = -1.0;
  rangeMax                    = +1.0;
  spectralContrast            = 1.0;
  spectralSlope               = 0.0;
  evenOddRatio                = 0.5; 
  phaseScale                  = 1.0;
  phaseShift                  = 0.0;
  evenOddPhaseShift           = 0.0;
  stereoPhaseShift            = 0.0;
  evenOddStereoPhaseShift     = 0.0;
  combHarmonic                = 1.0;
  //combOffset              = 0.0;
  combAmount                  = 0.0;
  lowestHarmonicToKeep        = 0; 
  highestHarmonicToKeep       = tableLength/2;
  timeReverse                 = false;
  polarityInvert              = false;
  meanRemove                  = false;
  normalize                   = false;
  fitToRange                  = false;
  autoUpdate                  = true;

  slewRateLimiter.setSampleRate(tableLength);
  slewRateLimiter.setAttackTime(0.0);
  slewRateLimiter.setReleaseTime(0.0);

  fourierTransformer.setBlockSize(tableLength);
  fourierTransformer.setDirection(FourierTransformerRadix2::FORWARD);
  fourierTransformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);

  //fillWithZeros();
  //renderWaveform(); // will initially create a sine-wave
  updateBuffers();
}

WaveTable::~WaveTable()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void WaveTable::getWaveform(double *targetBufferL, double *targetBufferR, 
                            int targetBufferLength)
{
  if( targetBufferLength == tableLength )
  {
    RAPT::rsArrayTools::copy(waveBufferL, targetBufferL, tableLength);
    RAPT::rsArrayTools::copy(waveBufferR, targetBufferR, tableLength);
  }
  else
  {
    RAPT::rsArrayTools::copyBufferWithLinearInterpolation(waveBufferL, tableLength, targetBufferL, 
      targetBufferLength);
    RAPT::rsArrayTools::copyBufferWithLinearInterpolation(waveBufferR, tableLength, targetBufferR, 
      targetBufferLength);
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void WaveTable::updateBuffers()
{
  waveformRenderer.renderWaveForm(prototypeBufferL, tableLength);
  RAPT::rsArrayTools::copy(prototypeBufferL, prototypeBufferR, tableLength);
  // \todo: let the waveformRenderer directly create a stereo waveform

  int    n, k;
  double nw;
  double tmpL[tableLength];
  double tmpR[tableLength];

  //-----------------------------------------------------------------
  // time-domain manipulations:

  for(n=0; n<tableLength; n++)
  {
    nw      = warpPhaseIndex(fmod((double)n, (double)tableLength));
    tmpL[n] = getPrototypeValueAt(0, nw);
    tmpR[n] = getPrototypeValueAt(1, nw);
  }
  if( timeReverse == true )
  {
    RAPT::rsArrayTools::reverse(tmpL, tableLength);
    RAPT::rsArrayTools::reverse(tmpR, tableLength);
  }

  //-----------------------------------------------------------------
  // spectral manipulations (magnitude and phase):

  Complex *spectrumL = new Complex[tableLength];
  Complex *spectrumR = new Complex[tableLength];
  double  *magL      = new double[tableLength/2];
  double  *magR      = new double[tableLength/2];
  double  *phsL      = new double[tableLength/2];
  double  *phsR      = new double[tableLength/2];

  // obtain magnitude and phase:
  RAPT::rsArrayTools::convertBuffer(tmpL, spectrumL, tableLength);
  RAPT::rsArrayTools::convertBuffer(tmpR, spectrumR, tableLength);
  fourierTransformer.setDirection(FourierTransformerRadix2::FORWARD);
  fourierTransformer.transformComplexBufferInPlace(spectrumL);
  fourierTransformer.transformComplexBufferInPlace(spectrumR);
  fourierTransformer.getRealSignalMagnitudesAndPhases(tmpL, magL, phsL);
  fourierTransformer.getRealSignalMagnitudesAndPhases(tmpR, magR, phsR);
 
  double maxMag = RAPT::rsArrayTools::maxValue(magL, tableLength/2);
  maxMag = RAPT::rsMax(maxMag, RAPT::rsArrayTools::maxValue(magR, tableLength/2)); // maximum magnitude

  double contrastNormalizer; 
  if( maxMag != 0.0 )
    contrastNormalizer = maxMag / pow(maxMag, spectralContrast);
  else 
    contrastNormalizer = 1.0;
  double slopeNormalizer    = 1.0;
  if( spectralSlope > 0.0 )
    slopeNormalizer = 1.0 / RAPT::rsDbToAmp( 0.5*spectralSlope*log2(0.5*tableLength) )   ;
      // good for waveforms with 1/n falloff of the harmonics

  double phi              = RAPT::rsDegreeToRadiant(phaseShift);
  double phiEvenOdd       = RAPT::rsDegreeToRadiant(0.5*evenOddPhaseShift);
  double phiStereo        = RAPT::rsDegreeToRadiant(0.5*stereoPhaseShift);
  double phiEvenOddStereo = RAPT::rsDegreeToRadiant(0.5*evenOddStereoPhaseShift);
  double evenAmp, oddAmp, phiL, phiR;
  if( evenOddRatio > 0.5 )
  {
    oddAmp  = 1.0;
    evenAmp = 2.0 - 2.0*evenOddRatio;
  }
  else
  {
    evenAmp = 1.0;
    oddAmp  = 2.0 * evenOddRatio;
  }
  double weight = 1.0;  // this is used for the magnitude weighting

  // the special case k=0 (DC and Nyquist-frequency) must be treated separately:
  weight   = contrastNormalizer*slopeNormalizer*evenAmp;
  if( lowestHarmonicToKeep > 0 )
    weight = 0.0;
  magL[0] *= weight;
  magR[0] *= weight;
  weight   = contrastNormalizer*slopeNormalizer*evenAmp*RAPT::rsDbToAmp(spectralSlope*log2(tableLength/2));
  if( highestHarmonicToKeep < tableLength/2 )
    weight = 0.0;
  phsL[0] *= weight;
  phsR[0] *= weight;

  // loop over the bins:
  for(k=1; k<tableLength/2; k++)   
  {
    // apply contrast function to spectral magnitudes:
    magL[k] = pow(magL[k], spectralContrast);
    magR[k] = pow(magR[k], spectralContrast);

    // calculate weight for the magnitude at this bin:
    weight = contrastNormalizer * slopeNormalizer * RAPT::rsDbToAmp(spectralSlope*log2(k));
    if( RAPT::rsIsEven(k) )
    {
      weight *= evenAmp;
      phiL    = phiEvenOdd - phiStereo - phiEvenOddStereo;
      phiR    = phiEvenOdd + phiStereo + phiEvenOddStereo;
    }
    else
    {
      weight *= oddAmp;
      phiL    = -phiEvenOdd - phiStereo + phiEvenOddStereo;
      phiR    = -phiEvenOdd + phiStereo - phiEvenOddStereo;
    }
    if( k < lowestHarmonicToKeep || k > highestHarmonicToKeep )
      weight = 0.0;

    // apply magnitude weighting:
    magL[k] *= weight;
    magR[k] *= weight;

    // apply phase modifications:
    phsL[k] = phaseScale*phsL[k] + phi + phiL;
    phsR[k] = phaseScale*phsR[k] + phi + phiR;
  }

  fourierTransformer.setDirection(FourierTransformerRadix2::INVERSE);
  fourierTransformer.getRealSignalFromMagnitudesAndPhases(magL, phsL, tmpL);
  fourierTransformer.getRealSignalFromMagnitudesAndPhases(magR, phsR, tmpR);

  delete[] spectrumL;
  delete[] spectrumR;
  delete[] magL;
  delete[] magR;
  delete[] phsL;
  delete[] phsR;

  //-----------------------------------------------------------------
  // filtering/smoothing manipulations:

  double combOffset = 0.5*tableLength/ combHarmonic;
  RAPT::rsArrayTools::addCircularShiftedCopy(tmpL, tableLength, combOffset, 0.01*combAmount);
  RAPT::rsArrayTools::addCircularShiftedCopy(tmpR, tableLength, combOffset, 0.01*combAmount);

  RAPT::rsArrayTools::copy(tmpL, waveBufferL, tableLength);
  RAPT::rsArrayTools::copy(tmpR, waveBufferR, tableLength);
  const int numWarmUpCycles = 3;
  slewRateLimiter.reset();
  for(int c=0; c<numWarmUpCycles+1; c++)
  {
    for(n=0; n<tableLength; n++)
      tmpL[n] = slewRateLimiter.getSample(waveBufferL[n]);
  }
  slewRateLimiter.reset();
  for(int c=0; c<numWarmUpCycles+1; c++)
  {
    for(n=0; n<tableLength; n++)
      tmpR[n] = slewRateLimiter.getSample(waveBufferR[n]);
  }

  //-----------------------------------------------------------------
  // range manipulations:

  if( polarityInvert == true )
  {
    RAPT::rsArrayTools::scale(tmpL, tmpL, tableLength, -1.0);
    RAPT::rsArrayTools::scale(tmpR, tmpR, tableLength, -1.0);
  }
  if( meanRemove == true )
  {
    RAPT::rsArrayTools::removeMean(tmpL, tableLength);
    RAPT::rsArrayTools::removeMean(tmpR, tableLength);
  }
  if( normalize == true )
  {
    RAPT::rsArrayTools::normalize(tmpL, tableLength, 1.0);
    RAPT::rsArrayTools::normalize(tmpR, tableLength, 1.0);
  }
  if( fitToRange == true )
  {
    RAPT::rsArrayTools::fitIntoRange(tmpL, tableLength, rangeMin, rangeMax);
    RAPT::rsArrayTools::fitIntoRange(tmpR, tableLength, rangeMin, rangeMax);
  }

  //-----------------------------------------------------------------
  // start-phase adjustment:

  double shift = -tableLength* startPhase                  /360.0;
  RAPT::rsArrayTools::circularShiftInterpolated(tmpL, tableLength, shift);
  shift        = -tableLength*(startPhase+rightPhaseOffset)/360.0;
  RAPT::rsArrayTools::circularShiftInterpolated(tmpR, tableLength, shift);

  //-----------------------------------------------------------------
  // manipulations done - copy temporary buffers into member-buffers:

  RAPT::rsArrayTools::copy(tmpL, waveBufferL, tableLength);
  RAPT::rsArrayTools::copy(tmpR, waveBufferR, tableLength);

  waveBufferL[tableLength] = waveBufferL[0]; // for interpolation
  waveBufferR[tableLength] = waveBufferR[0]; 
}

void WaveTable::fillWithZeros()
{
  for(int i=0; i<tableLength; i++)
    prototypeBufferL[i] = prototypeBufferR[i] = 0.0;
  updateBuffers();
}

double WaveTable::getPrototypeValueAt(int channel, double phaseIndex)
{
  int    i = floorInt(phaseIndex);                 // integer part
  double f = phaseIndex - i;                       // fractional part
  int    j = RAPT::rsWrapAround(i+1, tableLength); // next sample position

  double x0, x1;
  if( channel == 0 )
  {
    x0 = prototypeBufferL[i];
    x1 = prototypeBufferL[j];
  }
  else
  {
    x0 = prototypeBufferR[i];
    x1 = prototypeBufferR[j];
  }

  return (1.0-f)*x0 + f*x1;                // linear interpolation
}

double WaveTable::warpPhaseIndex(double unwarpedIndex)
{
  /*
  double a = timeWarp;
  double x = unwarpedIndex / (double) prototypeWaveNumSamples;
  double y = (x-a*x) / (1.0-2.0*a*x+a);
  return y * (double) prototypeWaveNumSamples;
  */
  double tmp = unwarpedIndex / (double) tableLength;   // in 0...+1
  while( tmp < 0.0 )
    tmp += 1.0;
  while( tmp >= 1.0 )
    tmp -= 1.0;

  double a   = fullWaveWarp;
  double b   = pow(20.0, halfWaveWarp);

  tmp = 2.0*tmp-1.0;                            // in -1...+1
  tmp = RAPT::rsSign(tmp) * pow(fabs(tmp), b);  // in -1...+1
  tmp = (tmp-a) / (1.0-a*tmp);                  // in -1...+1
  tmp = 0.5*(tmp+1);                            // in  0...+1

  return tmp * (double) tableLength;

  //double x1 = unwarpedIndex / (double) prototypeWaveNumSamples;
}











