//#include "rosic_MipMappedWaveTableStereo.h"
//using namespace rosic;

MipMappedWaveTableStereo::MipMappedWaveTableStereo()
{
  // set up the fourier-transformers and interpolator:
  forwardTransformer.setBlockSize(tableLength);
  forwardTransformer.setDirection(FourierTransformerRadix2::FORWARD);
  forwardTransformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);

  inverseTransformer.setBlockSize(tableLength);
  inverseTransformer.setDirection(FourierTransformerRadix2::INVERSE);
  inverseTransformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);

  //interpolator.setInterpolationMethod(Interpolator::HERMITE_4P_3O);

  sampleName                  = NULL;
  prototypeWave[0]            = NULL;
  prototypeWave[1]            = NULL;
  prototypeWaveNumSamples     = 0;
  channelSumAndDifferenceMode = false;
  timeReversed                = false;
  polarityInverted            = false;
  fullWavePhaseWarp           = 0.0;
  halfWavePhaseWarp           = 0.0;
  combHarmonic                = 1.0;
  combOffset                  = 0.0;
  combAmount                  = 0.0;
  evenOddRatio                = 0.5;
  highestHarmonicToKeep       = 1024;
  lowestHarmonicToKeep        = 1;
  spectralContrast            = 1.0;
  spectralSlope               = 0.0;
  phaseScale                  = 1.0;
  phaseShift                  = 0.0;
  evenOddPhaseShift           = 0.0;
  stereoPhaseShift            = 0.0;
  evenOddStereoPhaseShift     = 0.0;
  autoReRenderMipMap          = true;

  // initialize the buffers:
  fillWithAllZeros();
}

MipMappedWaveTableStereo::~MipMappedWaveTableStereo()
{
  if( sampleName != NULL )
    delete[] sampleName;
  if( prototypeWave[0] != NULL )
    delete[] prototypeWave[0];
  if( prototypeWave[1] != NULL )
    delete[] prototypeWave[1];
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void MipMappedWaveTableStereo::setWaveform(double** newWaveForm, int lengthInSamples)
{
  mutex.lock();

  // copy the data into internal arrays:
  if( lengthInSamples != prototypeWaveNumSamples )
  {
    prototypeWaveNumSamples = lengthInSamples;
    if( prototypeWave[0] != NULL )
      delete[] prototypeWave[0];
    if( prototypeWave[1] != NULL )
      delete[] prototypeWave[1];
    prototypeWave[0] = new double[prototypeWaveNumSamples];
    prototypeWave[1] = new double[prototypeWaveNumSamples];
  }
  for(int n=0; n<prototypeWaveNumSamples; n++)
  {
    prototypeWave[0][n] = newWaveForm[0][n];
    prototypeWave[1][n] = newWaveForm[1][n];
  }

  // generate the multisample:
  if( autoReRenderMipMap == true )
    renderMipMap();

  mutex.unlock();
}

void MipMappedWaveTableStereo::setWaveform(float** newWaveForm, int lengthInSamples)
{
  mutex.lock();

  // typecast and copy the data into internal arrays:
  if( lengthInSamples != prototypeWaveNumSamples )
  {
    prototypeWaveNumSamples = lengthInSamples;
    if( prototypeWave[0] != NULL )
      delete[] prototypeWave[0];
    if( prototypeWave[1] != NULL )
      delete[] prototypeWave[1];
    prototypeWave[0] = new double[prototypeWaveNumSamples];
    prototypeWave[1] = new double[prototypeWaveNumSamples];
  }
  for(int n=0; n<prototypeWaveNumSamples; n++)
  {
    prototypeWave[0][n] = (double) newWaveForm[0][n];
    prototypeWave[1][n] = (double) newWaveForm[1][n];
  }
  if( autoReRenderMipMap == true )
    renderMipMap();

  mutex.unlock();
}

void MipMappedWaveTableStereo::setSampleName(char *newSampleName)
{
  mutex.lock();

  // free old and allocate new memory for the name:
  if( sampleName != NULL )
  {
    delete[] sampleName;
    sampleName = NULL;
  }
  if( newSampleName != NULL )
  {
    int newLength    = (int) strlen(newSampleName);
    sampleName       = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      sampleName[c] = newSampleName[c];
  }

  mutex.unlock();
}

void MipMappedWaveTableStereo::setTimeReverse(bool shouldBeReversed)
{
  timeReversed = shouldBeReversed;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setPolarityInversion(bool shouldBeInversed)
{
  polarityInverted = shouldBeInversed;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setFullWavePhaseWarp(double newWarpCoefficient)
{
  fullWavePhaseWarp = newWarpCoefficient;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setHalfWavePhaseWarp(double newWarpCoefficient)
{
  halfWavePhaseWarp = newWarpCoefficient;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setCombHarmonic(double newHarmonic)
{
  combHarmonic = newHarmonic;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

/*
void MipMappedWaveTableStereo::setCombOffset(double newCombOffset)
{
  combOffset = newCombOffset;
  if( autoReRenderMipMap == true )
    renderMipMap();
}
*/

void MipMappedWaveTableStereo::setCombAmount(double newCombAmount)
{
  combAmount = newCombAmount;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setChannelSumAndDifferenceMode(bool shouldUseSumAndDifference)
{
  channelSumAndDifferenceMode = shouldUseSumAndDifference;
}

void MipMappedWaveTableStereo::setEvenOddRatio(double newRatio)
{
  evenOddRatio = newRatio;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setHighestHarmonicToKeep(int newHighestHarmonicToKeep)
{
  highestHarmonicToKeep = newHighestHarmonicToKeep;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setLowestHarmonicToKeep(int newLowestHarmonicToKeep)
{
  lowestHarmonicToKeep = newLowestHarmonicToKeep;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setSpectralContrast(double newContrast)
{
  spectralContrast = newContrast;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setSpectralSlope(double newSlope)
{
  spectralSlope = newSlope;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setEvenOddPhaseShift(double newPhaseShift)
{
  evenOddPhaseShift = newPhaseShift;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setPhaseScale(double newPhaseScale)
{
  phaseScale = newPhaseScale;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setPhaseShift(double newPhaseShift)
{
  phaseShift = newPhaseShift;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setStereoPhaseShift(double newPhaseShift)
{
  stereoPhaseShift = newPhaseShift;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setEvenOddStereoPhaseShift(double newPhaseShift)
{
  evenOddStereoPhaseShift = newPhaseShift;
  if( autoReRenderMipMap == true )
    renderMipMap();
}

void MipMappedWaveTableStereo::setAutomaticMipMapReRendering(bool shouldAutomaticallyReRender)
{
  autoReRenderMipMap = shouldAutomaticallyReRender;
}

//-------------------------------------------------------------------------------------------------
// Inquiry:
/*
double** MipMappedWaveTableStereo::getPrototypeWaveform()
{
  return prototypeWave;
}
*/

int MipMappedWaveTableStereo::getPrototypeNumSamples() /*const*/
{
  mutex.lock();
  int result = prototypeWaveNumSamples;
  mutex.unlock();

  return result;
}

void MipMappedWaveTableStereo::copyDataTo(double* buffer, int channel, int level) /*const*/
{
  RAPT::rsAssert(channel >= 0 && channel < 2);
  RAPT::rsAssert(level   >= 0 && level   < getNumLevels());
  mutex.lock();
  for(int i = 0; i < getTableLength(); i++)
    buffer[i] = tableSet[level][i][channel];
  mutex.unlock();
}


/*
double MipMappedWaveTableStereo::getFullWavePhaseWarp()
{
  return fullWavePhaseWarp;
}

double MipMappedWaveTableStereo::getHalfWavePhaseWarp()
{
  return halfWavePhaseWarp;
}

bool MipMappedWaveTableStereo::isMipMapAutoReRenderingActive()
{
  return autoReRenderMipMap;
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

void MipMappedWaveTableStereo::renderMipMap()
{
  if( prototypeWave[0] == NULL || prototypeWave[1] == NULL )
    return;

  mutex.lock();

  int t, k, n; // indices for the channel, table and position

  // Allocate some memory for temporary storage:
  Complex* inSpectrumL  = new Complex[prototypeWaveNumSamples];
  Complex* inSpectrumR  = new Complex[prototypeWaveNumSamples];
  Complex* tmpSpectrumL = new Complex[tableLength];
  Complex* tmpSpectrumR = new Complex[tableLength];
  double*  tmpTableL    = new double[tableLength];
  double*  tmpTableR    = new double[tableLength];

  // Copy the prototype waveform into the complex array which is going to used for the forward FFT,
  // thereby apply time-domain manipulations of the waveform:
  double signFactor = 1.0;
  double readPhase;
  if( polarityInverted )
    signFactor = -1.0;  

  double combOffset          = 180.0 / combHarmonic;
  double combOffsetInSamples = prototypeWaveNumSamples * combOffset / 360.0;
  //double combOffsetNormalized = combOffset / 360.0; // 0...1
  if( timeReversed )
  {
    for(n=0; n<prototypeWaveNumSamples; n++)
    {
      readPhase       = warpPhaseIndex((double) (prototypeWaveNumSamples-n-1));
      inSpectrumL[n]  = Complex( signFactor*getPrototypeValueAt(0, readPhase) );
      inSpectrumR[n]  = Complex( signFactor*getPrototypeValueAt(1, readPhase) );

      //readPhase      += combOffsetInSamples;
      readPhase       = warpPhaseIndex((double) (prototypeWaveNumSamples-n-combOffsetInSamples-1));
      if( readPhase >= prototypeWaveNumSamples )
        readPhase -= prototypeWaveNumSamples;
      inSpectrumL[n] += 0.01*combAmount * Complex( signFactor*getPrototypeValueAt(0, readPhase) );
      inSpectrumR[n] += 0.01*combAmount * Complex( signFactor*getPrototypeValueAt(1, readPhase) );
    }
  }
  else
  {
    for(n=0; n<prototypeWaveNumSamples; n++)
    {
      readPhase       = warpPhaseIndex((double) n);
      inSpectrumL[n]  = Complex( signFactor*getPrototypeValueAt(0, readPhase) );
      inSpectrumR[n]  = Complex( signFactor*getPrototypeValueAt(1, readPhase) );

      //readPhase      += combOffsetInSamples;
      readPhase       = warpPhaseIndex((double) (n+combOffsetInSamples));
      if( readPhase >= prototypeWaveNumSamples )
        readPhase -= prototypeWaveNumSamples;
      inSpectrumL[n] += 0.01*combAmount * Complex( signFactor*getPrototypeValueAt(0, readPhase) );
      inSpectrumR[n] += 0.01*combAmount * Complex( signFactor*getPrototypeValueAt(1, readPhase) );
    }
  }

  // Set the forward transfrom object to the size of the new data-array and transform the
  // prototypes into the frequency domain:
  forwardTransformer.setBlockSize(prototypeWaveNumSamples);
  forwardTransformer.transformComplexBufferInPlace(inSpectrumL);
  forwardTransformer.transformComplexBufferInPlace(inSpectrumR);

  // Calculate the number of FFT-bins to fill (including the redundant bins):
  int m = RAPT::rsMin(prototypeWaveNumSamples, tableLength);
  if( RAPT::rsIsOdd(m) )
    m -= 1;

  // Convert the real/imaginary representation of the spectrum into magnitude/phase and find the
  // peak magnitude (required array size would be m/2+1, however allocating arrays of size m seems
  // to be much better performance wise):
  double  maxMagnitude = 0.0;
  double* magSpectrumL = new double[m];
  double* magSpectrumR = new double[m];
  double* phsSpectrumL = new double[m];
  double* phsSpectrumR = new double[m];
  for(k=0; k <= m/2; k++)
  {
    magSpectrumL[k] = inSpectrumL[k].getRadius();
    magSpectrumR[k] = inSpectrumR[k].getRadius();
    phsSpectrumL[k] = inSpectrumL[k].getAngle();
    phsSpectrumR[k] = inSpectrumR[k].getAngle();
    if( magSpectrumL[k] > maxMagnitude )
      maxMagnitude = magSpectrumL[k];
    if( magSpectrumR[k] > maxMagnitude )
      maxMagnitude = magSpectrumR[k];
  }

  // Calculate some weighting factors:
  double contrastNormalizer = maxMagnitude / pow(maxMagnitude, spectralContrast);

  double slopeNormalizer = 1.0;
  if( spectralSlope > 0.0 )
    slopeNormalizer = 1.0 / RAPT::rsDbToAmp( 0.5*spectralSlope*log2(0.5*m) )   ;
    // this normalizer is good for waveforms with 1/n falloff of the harmonics

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

  tmpSpectrumL[0] = inSpectrumL[0];
  tmpSpectrumR[0] = inSpectrumR[0];
  double weight   = 1.0;
  for(k=1; k <= m/2; k++)
  {
    // Apply contrast function to spectral magnitudes:
    magSpectrumL[k] = pow(magSpectrumL[k], spectralContrast);
    magSpectrumR[k] = pow(magSpectrumR[k], spectralContrast);

    // Calculate weight for the magnitude at this bin:
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

    // For debug - should hold when no modifications are applied:
    //RAPT::rsAssert(weight == 1 && phi == 0 && phiL == 0 && phiR == 0);

    // Apply magnitude weighting:
    magSpectrumL[k] *= weight;
    magSpectrumR[k] *= weight;

    // Apply phase modifications:
    phsSpectrumL[k] = phaseScale*phsSpectrumL[k] + phi + phiL;
    phsSpectrumR[k] = phaseScale*phsSpectrumR[k] + phi + phiR;

    // Establish the complex spectral value from the magnitude and phase:
    tmpSpectrumL[k].setRadiusAndAngle(magSpectrumL[k], phsSpectrumL[k]);
    tmpSpectrumR[k].setRadiusAndAngle(magSpectrumR[k], phsSpectrumR[k]);

    // Convert the L/R input data into sum and difference, if this option is chosen:
    if( channelSumAndDifferenceMode )
    {
      Complex tmpS    = tmpSpectrumL[k]+tmpSpectrumR[k];
      Complex tmpD    = tmpSpectrumL[k]-tmpSpectrumR[k];
      tmpSpectrumL[k] = tmpS;
      tmpSpectrumR[k] = tmpD;
    }

    // Symmetrize:
    tmpSpectrumL[tableLength-k] = tmpSpectrumL[k];
    tmpSpectrumR[tableLength-k] = tmpSpectrumR[k];
  }

  // Progressively truncate spectrum and render waveform via iFFT:
  int highBin = tableLength/2;
  int lowBin  = tableLength/4;
  for(t=0; t<numTables; t++)
  {
    // Write the inverse fourier transform of the current spectrum into the temporary tables:
    inverseTransformer.transformSymmetricSpectrum(tmpSpectrumL, tmpTableL);
    inverseTransformer.transformSymmetricSpectrum(tmpSpectrumR, tmpTableR);

    // Copy the temporary table into the member, thereby typecast:
    for(k=0; k<tableLength; k++)
    {
      tableSet[t][k][0] = (float) tmpTableL[k];
      tableSet[t][k][1] = (float) tmpTableR[k];
    }

    // Repeat some samples form the beginning for the interpolator:
    tableSet[t][tableLength+0][0] = tableSet[t][0][0];
    tableSet[t][tableLength+1][0] = tableSet[t][1][0];
    tableSet[t][tableLength+2][0] = tableSet[t][2][0];
    tableSet[t][tableLength+3][0] = tableSet[t][3][0];

    tableSet[t][tableLength+0][1] = tableSet[t][0][1];
    tableSet[t][tableLength+1][1] = tableSet[t][1][1];
    tableSet[t][tableLength+2][1] = tableSet[t][2][1];
    tableSet[t][tableLength+3][1] = tableSet[t][3][1];

    // Truncate the spectrum for the next iteration:
    for(k = lowBin+1; k <= highBin; k++)
    {
      tmpSpectrumL[k]             = 0.0;
      tmpSpectrumR[k]             = 0.0;
      tmpSpectrumL[tableLength-k] = 0.0;
      tmpSpectrumR[tableLength-k] = 0.0;
    }
    highBin /= 2;
    lowBin  /= 2;
  }

  // Free temporary allocated memory:
  delete[] inSpectrumL;
  delete[] inSpectrumR;
  delete[] tmpSpectrumL;
  delete[] tmpSpectrumR;
  delete[] magSpectrumL;
  delete[] magSpectrumR;
  delete[] phsSpectrumL;
  delete[] phsSpectrumR;
  delete[] tmpTableL;
  delete[] tmpTableR;

  mutex.unlock();
}

double MipMappedWaveTableStereo::getPrototypeValueAt(int channel, double phaseIndex)
{
  int    i  = floorInt(phaseIndex);
  double f  = phaseIndex - i;
  double x0 = prototypeWave[channel][i];
  double x1;
  if( (i+1) < prototypeWaveNumSamples )
    x1 = prototypeWave[channel][i+1];
  else
    x1 = prototypeWave[channel][0]; // wraparound for last value

  return (1.0-f)*x0 + f*x1;         // linear interpolation
}

double MipMappedWaveTableStereo::warpPhaseIndex(double unwarpedIndex)
{
  /*
  double a = timeWarp;
  double x = unwarpedIndex / (double) prototypeWaveNumSamples;
  double y = (x-a*x) / (1.0-2.0*a*x+a);
  return y * (double) prototypeWaveNumSamples;
  */
  double tmp = unwarpedIndex / (double) prototypeWaveNumSamples;   // in 0...+1
  while( tmp < 0.0 )
    tmp += 1.0;
  while( tmp >= 1.0 )
    tmp -= 1.0;

  double a   = fullWavePhaseWarp;
  double b   = pow(20.0, halfWavePhaseWarp);

  tmp = 2.0*tmp-1.0;                            // in -1...+1
  tmp = RAPT::rsSign(tmp) * pow(fabs(tmp), b);  // in -1...+1
  tmp = (tmp-a) / (1.0-a*tmp);                  // in -1...+1
  tmp = 0.5*(tmp+1);                            // in  0...+1

  return tmp * (double) prototypeWaveNumSamples;

  //double x1 = unwarpedIndex / (double) prototypeWaveNumSamples;
}

void MipMappedWaveTableStereo::fillWithAllZeros()
{
  int c, t, n; // indices fo table and position
  for(c=0; c<numChannels; c++)
  {
    for(t=0; t<numTables; t++)
    {
      for(n=0; n<tableLength+4; n++)
      {
        tableSet[t][n][c] = 0.0;
      }
    }
  }
}



