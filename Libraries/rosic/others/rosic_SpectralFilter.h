#ifndef rosic_SpectralFilter_h
#define rosic_SpectralFilter_h

// rosic-indcludes:
#include "rosic_SpectralProcessor.h"

namespace rosic
{

  /**

  This class implements a filter that operates in the spectral (i.e. FFT) domain.

  */

  class SpectralFilter : public SpectralProcessor
  {

  public:

    enum modes
    {
      BYPASS = 0,
      BANDPASS,
      BANDREJECT
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum 
    overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input, 
    an overlap and zero-padding of 2 is usually a good choice. */
    SpectralFilter(int maxBlockSize, int maxOverlapFactor = 4, int maxPaddingFactor = 4); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the samplerate in Hz. */
    void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; }

    /** Sets the lower cutoff frequency in Hz. */
    void setLowerCutoff(double newCutoff) { lowerCutoff = newCutoff; }

    /** Sets the upper cutoff frequency in Hz. */
    void setUpperCutoff(double newCutoff) { upperCutoff = newCutoff; }

    /** Sets the mode of the filter. @see modes. */
    void setMode(int newMode) { mode = newMode; }

    //=============================================================================================

  protected:

    /** This function does the actual filtering of the spectrum. */
    virtual void processSpectrum(Complex *spectrum, int spectrumSize);

    double lowerCutoff, upperCutoff, sampleRate;
    int    mode;

  };

} // end namespace rosic

#endif // rosic_SpectralFilter_h
