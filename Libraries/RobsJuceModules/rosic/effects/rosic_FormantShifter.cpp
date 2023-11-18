//#include "rosic_FormantShifter.h"
////#include "rosic_Plotter.h"  // for test/debug
//using namespace rosic;

//=================================================================================================
// class FormantShifter:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FormantShifter::FormantShifter(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: SpectralEnvelopeProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  scale              = 1.0;
  offset             = 0.0;
  energyCompensation = true;
}

//-------------------------------------------------------------------------------------------------
// internal fucntions:

void FormantShifter::processSpectrum(Complex *spectrum, int spectrumSize)
{
  int k;
  double *mag   = new double[spectrumSize];  // magnitudes
  double *env   = new double[spectrumSize];  // spectral envelope
  double energy = 0.0;
  for(k=0; k<spectrumSize; k++)   
  {
    mag[k]  = spectrum[k].getRadius();
    energy += mag[k]*mag[k];
  }

  estimateSpectralEnvelope(mag, env, spectrumSize);

  // take out the spectral envelope:
  double c = 0.0000001; // to avoid division by zero
  for(k=0; k<spectrumSize; k++)
    spectrum[k] /= (env[k] + c);

  // modify the spectral envelope:
  double *tmp = new double[spectrumSize];  // temporary spectral envelope storage
  for(k=0; k<spectrumSize; k++)
    tmp[k] = env[k];

  // this loop may be optimized ... pre-calculations of divisisors..
  for(k=0; k<spectrumSize; k++)
  {
    double fw    = k*sampleRate / (2*spectrumSize);   // frequency to write into
    double fr    = (fw-offset)  / scale;              // frequency to read from      
    double kr    = (2*spectrumSize*fr) / sampleRate;  // bin index to read from
    int    kInt  = floorInt(kr);
    double kFrac = kr-kInt;
    if( kInt >= 0 && kInt < spectrumSize-1 )
      env[k] = (1.0-kFrac) * (tmp[kInt]) + kFrac * (tmp[kInt+1]);
    else
      env[k] = 0.0;
  }

  // re-apply the (modified) spectral envelope:
  for(k=0; k<spectrumSize; k++)
    spectrum[k] *= (env[k] + c);

  // restore the original energy of the spectrum, if so desired:
  if( energyCompensation == true )
  {
    double newEnergy = 0.0;
    for(k=0; k<spectrumSize; k++)   
      newEnergy += spectrum[k].re * spectrum[k].re + spectrum[k].im * spectrum[k].im;
    double compensator = sqrt(energy / (newEnergy+c));
    for(k=0; k<spectrumSize; k++)
      spectrum[k] *= compensator;
  }

  /*
  // plotting - only for debug/check:
  static int numCalls = 1;
  if( numCalls == 30 )
  {
    // plot the magnitude spectrum and spectral envelope:
    double *freqs = new double[spectrumSize];
    for(k=0; k<spectrumSize; k++)
      freqs[k] = k*sampleRate / (2*spectrumSize);
    Plotter::plotData(spectrumSize/4, freqs, mag, tmp, env);
  }
  numCalls++;  
  */

  delete[] mag;
  delete[] env;
  delete[] tmp;

  // ToDo:
  // -Get rid of the memory allocation to make the implementation realtime ready.
}



//=================================================================================================
// class FormantShifterStereo:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FormantShifterStereo::FormantShifterStereo(int maxBlockSize, int maxOverlapFactor, 
                                           int maxPaddingFactor)
: shifterL(maxBlockSize, maxOverlapFactor, maxPaddingFactor),
  shifterR(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  dry = 0.0;
  wet = 1.0;
}

FormantShifterStereo::~FormantShifterStereo()
{

}