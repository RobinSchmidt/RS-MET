#ifndef rosic_FourierTransformerBluestein_h
#define rosic_FourierTransformerBluestein_h

//// rosic-indcludes:
//#include "rosic_FourierTransformerRadix2.h"
//#include "../math/rosic_ComplexFunctions.h"

namespace rosic
{

  /**

  This class works similar to the FourierTransformerRadix2 class, except that it accepts arbitrary
  blockSizes. It uses Bluestein's algorithm to evaluate the DFT of a given block which works for 
  arbitrary block-sizes and preserves the desirable O(N*log(N)) complexity of FFT algorithms 
  (although the implicit constants are larger than for radix-2 algorithms).

  \todo implement normalization

  */

  class FourierTransformerBluestein  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    FourierTransformerBluestein();  

    /** Destructor. */
    ~FourierTransformerBluestein(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** FFT-size, can be an arbitrary integer > 1. Powers of 2 will be most efficient, other 
    numbers increase the CPU-load significantly, but O(N*log(N))-scaling is preserved. */
    void setBlockSize(int newBlockSize);     

    /** Sets the direction of the transform (@see: FourierTransformerRadix2::directions). This 
    will affect the sign of the exponent (or equivalently: theimaginary part) in the twiddling 
    factors and the normalization constant. */
    void setDirection(int newDirection);

    /** Sets the mode for normalization of the output 
    (@see: FourierTransformerRadix2::normalizationModes). */
    void setNormalizationMode(int newNormalizationMode);

    //---------------------------------------------------------------------------------------------
    // signal processing:

    /** Transforms a buffer of complex numbers into its (forward or inverse) fourier transform. 
    The inBuffer will remain intact. Both, inBuffer and outBuffer must be of the size which was 
    specified when setting up the blockSize with setBlockSize(). */
    void transformComplexBuffer(Complex *inBuffer, Complex *outBuffer);   

    /** Does the same thing as transformComplexBuffer but performes in-place calculation 
    (overwrites the input buffer). */
    void transformComplexBufferInPlace(Complex *buffer);         

    //=============================================================================================

  protected:

    /** Generates the chirp signal which is needed modulate the input signal (and to retrieve the
    output-spectrum), as well as the spectrum of the conjugate chirp signal which is needed to 
    multiply the spectrum of the chirp-modulated signal. */
    void generateChirp();

    /** Updates the normalizationFactor member variable acording to a new blockSize, direction or
    normalizationMode. */
    void updateNormalizationFactor();

    Complex *h;             /**< The spectrum of the of the h-values. */
    Complex *c;             /**< The modulating chirp-signal. */
    Complex *y;             /**< Internal buffer of size M. */

    int N;                  /**< The blocksize of the Bluestein-FFT. */
    int M;                  /**< The enlarged blocksize for the embedded radix-2 FFT. M ist the 
                            smallest power of two, such that M >= 2*N-1. */

    int direction;          /**< The direction of the transform (@see: directions). */
    int normalizationMode;  /**< The normalization mode (@see: normalizationModes. */

    double normalizationFactor;

    bool blockSizeIsPowerOfTwo;
    /**< Indicates if we have the special case of a power of two blocksize - in this case, the 
    whole transformation is significantly simplified. */

    FourierTransformerRadix2 transformerRadix2;
    /**< This embedded object is used to perform the radix-2 forward and inverse transformations 
    which occur as part of the Bluestein algorithm. */

  };

} // end namespace rosic

#endif // rosic_FourierTransformerBluestein_h
