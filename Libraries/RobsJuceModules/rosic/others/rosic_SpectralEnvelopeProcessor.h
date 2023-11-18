#ifndef rosic_SpectralEnvelopeProcessor_h
#define rosic_SpectralEnvelopeProcessor_h

//// rosic-indcludes:
//#include "rosic_SpectralProcessor.h"
//#include "rosic_SlewRateLimiter.h"

namespace rosic
{

  /**

  This class serves as basclass for all sorts of spectral processors that need to work with the 
  spectral envelope. ...TBC...

  */

  class SpectralEnvelopeProcessor : public SpectralProcessor
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum 
    overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input, 
    an overlap and zero-padding of 2 is usually a good choice. */
    SpectralEnvelopeProcessor(int maxBlockSize, int maxOverlapFactor = 4, 
      int maxPaddingFactor = 4); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the samplerate in Hz. */
    void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; }

    /** Sets the frequency interval in which the spectral smoother falls down to 1/e - this 
    determines the frequency resolution of the estimated spectral envelope.*/
    void setSmoothingInterval(double newInterval) 
    { spectralSmoother.setReleaseTime(1000.0*newInterval); }

    //=============================================================================================

  protected:

    /** This function does the actual filtering of the spectrum. */
    virtual void processSpectrum(Complex *spectrum, int spectrumSize);

    /** This function estimates the spectral envelope and stores it in the output array 
    'spectralEnvelope'. */
    virtual void estimateSpectralEnvelope(double *spectralMagnitudes, double *spectralEnvelope, 
      int spectrumSize);

    double          smoothingInterval, sampleRate;
    SlewRateLimiter spectralSmoother;

    bool dBMode;

  };

} // end namespace rosic

#endif // rosic_SpectralEnvelopeProcessor_h
