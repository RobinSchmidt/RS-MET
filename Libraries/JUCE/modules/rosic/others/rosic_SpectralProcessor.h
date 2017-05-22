#ifndef rosic_SpectralProcessor_h
#define rosic_SpectralProcessor_h

//// rosic-indcludes:
//#include "rosic_OverlapAddProcessor.h"
//#include "../transforms/rosic_FourierTransformerRadix2.h"

namespace rosic
{

  /**

  This class serves a baseclass for effects that are based on spectral manipulations.

  */

  class SpectralProcessor : public OverlapAddProcessor
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the maximum blocksize here and may optionally pass a maximum 
    overlap- and zero-padding factor. For spectral processing with a cosine^2 window on the input, 
    an overlap and zero-padding of 2 is usually a good choice. */
    SpectralProcessor(int maxBlockSize, int maxOverlapFactor = 4, int maxPaddingFactor = 4); 

    /** Destructor. */
    ~SpectralProcessor();  

    //=============================================================================================

  protected:

    /** This function is the one, you should override in your subclass to do the actual processing.
    The baseclass implementation does nothing. */
    virtual void processBlock(double *block, int blockSize);

    /** This function is to be overriden by your subclass to process the spectrum in some way. The 
    passed 'spectrumSize' will be N/2 where N is the FFT-size - thus, the passed spectrum contains 
    only the non-redundant (positive) frequencies with bin-indices 0...N/2-1. The purely real 
    coefficients for DC and the Nyquist-frequency will be contained in the real and imaginary parts 
    of spectrum[0] respectively. */
    virtual void processSpectrum(Complex *spectrum, int spectrumSize) {}

    FourierTransformerRadix2 transformer;

  };

} // end namespace rosic

#endif // rosic_SpectralProcessor_h
