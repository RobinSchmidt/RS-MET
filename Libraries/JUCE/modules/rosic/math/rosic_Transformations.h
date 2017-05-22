#ifndef rosic_Transformations_h
#define rosic_Transformations_h

//// third party includes:
//#include "../_third_party/kiss_fft_v1_2_6/kiss_fft.h"

//// rosic includes:
//#include "rosic_Complex.h"
//#include "rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  // better use names similar to MatLab like fft, ifft
  void discreteFourierTransform(Complex *in, Complex *out, int length, bool isInverse); 
  /**< Computes the discrete fourier transform of an input sequence. */

  void radixTwoFastFourierTransform(Complex *in, Complex *out, int length, bool isInverse); 
  /**< Computes the discrete fourier transform of an input sequence which has a length of a power 
  of two (this is not checked for) . */


} // end namespace rosic

#endif // #ifndef rosic_Transformations_h