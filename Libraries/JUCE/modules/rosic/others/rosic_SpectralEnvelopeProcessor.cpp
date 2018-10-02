//#include "rosic_SpectralEnvelopeProcessor.h"
//#include "rosic_Plotter.h"  // for test/debug
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SpectralEnvelopeProcessor::SpectralEnvelopeProcessor(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: SpectralProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  dBMode      = false;
  sampleRate  = 44100.0;
  transformer.setBlockSize(paddingFactor*blockSize);
  spectralSmoother.setAttackTime(0.0);
  spectralSmoother.setReleaseTime(1000*200.0);  
    // these 'times' are actually frequency intervals in milli-Hertz
}

//-------------------------------------------------------------------------------------------------
// internal fucntions:

void SpectralEnvelopeProcessor::processSpectrum(Complex *spectrum, int spectrumSize)
{
  int k;
  double *mag = new double[spectrumSize];  // magnitudes
  double *env = new double[spectrumSize];  // spectral envelope
  for(k=0; k<spectrumSize; k++)   
    mag[k] = spectrum[k].getRadius();

  estimateSpectralEnvelope(mag, env, spectrumSize);

  for(k=0; k<spectrumSize; k++)
    spectrum[k] /= env[k];

  // do something with the envelope...

  for(k=0; k<spectrumSize; k++)
    spectrum[k] *= env[k];

  delete[] mag;
  delete[] env;
}

void SpectralEnvelopeProcessor::estimateSpectralEnvelope(double *mag, double *env, int spectrumSize)
{
  int    k;
  int    N           = 2*spectrumSize;
  double freqSpacing = sampleRate/N;
  spectralSmoother.setSampleRate(1.0/freqSpacing);

  if( dBMode == false )
  {
    spectralSmoother.reset();
    for(k=0; k<spectrumSize; k++)
      env[k] = spectralSmoother.getSample(mag[k]);
    spectralSmoother.reset();
    for(k=spectrumSize-1; k>=0; k--)
      env[k] = spectralSmoother.getSample(env[k]);
  }
  else
  {
    for(k=0; k<spectrumSize; k++)
      env[k] = RAPT::rsAmpToDb(mag[k]);
    spectralSmoother.reset();
    for(k=0; k<spectrumSize; k++)
      env[k] = spectralSmoother.getSample(env[k]);
    spectralSmoother.reset();
    for(k=spectrumSize-1; k>=0; k--)
      env[k] = spectralSmoother.getSample(env[k]);
    for(k=spectrumSize-1; k>=0; k--)
      env[k] = dB2amp(env[k]);
  }

  /*
  // plotting - only for debug/check:
  static int numCalls = 1;
  if( numCalls == 30 )
  {
    // plot the magnitude spectrum and spectral envelope:
    double *freqs = new double[spectrumSize];
    for(k=0; k<spectrumSize; k++)
      freqs[k] = k*freqSpacing;
    Plotter::plotData(spectrumSize/4, freqs, mag, env);
  }
  numCalls++;  
  */
}


