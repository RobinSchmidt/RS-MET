//#include "rosic_SpectralFilter.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SpectralFilter::SpectralFilter(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: SpectralProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  lowerCutoff = 250.0;
  upperCutoff = 4000.0;
  sampleRate  = 44100.0;
  mode        = BANDPASS;
  transformer.setBlockSize(paddingFactor*blockSize);
}

//-------------------------------------------------------------------------------------------------
// internal fucntions:

void SpectralFilter::processSpectrum(Complex *spectrum, int spectrumSize)
{
  int k;
  int N     = spectrumSize;    // these are the N positive frequencies only
  int loBin = (int)ceil( 2*spectrumSize*lowerCutoff/sampleRate); 
  int hiBin = (int)floor(2*spectrumSize*upperCutoff/sampleRate);
  if( mode == BANDPASS )
  {
    if( loBin >= 0 )
      spectrum[0].re = 0.0;   // set DC to 0
    if( hiBin <= N )
      spectrum[0].im = 0.0;   // set Nyquist to 0
    for(k=1; k<loBin; k++)
      spectrum[k] = 0.0;
    for(k=hiBin; k<N; k++)
      spectrum[k] = 0.0;
  }
  else if( mode == BANDREJECT )
  {
    if( loBin <= 0 )
    {
      spectrum[0].re = 0.0;
      loBin = 1;
    }
    if( hiBin >= N )
    {
      spectrum[0].im = 0.0;
      hiBin = N-1;
    }
    for(k=loBin; k<=hiBin; k++)
      spectrum[k]   = 0.0;
  }
}

/*

ToDo:
-Let the user select lower and upper transition widths
-For these transition bands, let the user select a transition shape (linear, cosine, cubic, 
 quartic, quintic, etc.). Maybe the shape can be selected independently for the region near the
 top and near the bottom, i.e. an asymmetric shape.
-Maybe let the user select a "bleed through" amount. This would give rise to shelving filters.

*/



