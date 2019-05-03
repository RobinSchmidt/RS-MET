#ifndef LowpassHighpassStereoSSE2_h
#define LowpassHighpassStereoSSE2_h

#include <intrin.h> 
#include <xmmintrin.h>
#include "MoreMath.h"
using namespace MoreMath;

/**

This class combines a first order lowpass-, and a first order highpass-filter
into a single object. The 2 filters, which are originally in a series 
connection, are implemented as a single biquad-stage.

*/

class LowpassHighpassStereoSSE2
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 LowpassHighpassStereoSSE2();   ///< Constructor.
 ~LowpassHighpassStereoSSE2();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setLpfCutoff (double newLpfCutoff);
 ///< Sets the cutoff frequency of the lowpass-filter.

 void setHpfCutoff (double newHpfCutoff);
 ///< Sets the cutoff frequency of the highpass-filter.

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE __m128d getSampleVector(__m128d inVector);
 ///< Calculates one output sample-frame at a time.

 //---------------------------------------------------------------------------
 // others:

 void resetBuffers();
 ///< Resets the internal buffers of the filter to zero.

 //===========================================================================

protected:

 __m128d x1;  // SSE2-vector of two doubles, representing left and right input
              // delayed by one sample
 __m128d x2;  // ...delayed by two samples
 __m128d y1;  // ...output delayed by one sample
 __m128d y2;  // ...output delayed by two samples

 double b0, b1, b2; // feedforward-coeffs
 double a1, a2;     // feedback-coeffs

 // filter parameters:
 doubleA lpfCutoff;
 doubleA hpfCutoff;

 doubleA sampleRate;
 doubleA sampleRateRec;    // reciprocal of the sampleRate
 
 // private member functions:
 void calcCoeffs();     
  // calculates the filter coefficients from filter parameters. The 
  // design-equations for the first order filters come from "The Scientist and
  // Engineers Guide to Digital Signal Processing" (www.dspguide.com)
};

//-----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):


INLINE __m128d LowpassHighpassStereoSSE2::getSampleVector(__m128d inVector)
{
 /*
 __m128d m1, m2, m3, m4, m5, m6, m7, m8;

 // create vector from double-b0-coeff, multiply with input, store:
 m1 = _mm_mul_pd(_mm_set1_pd(b0), in);
 m2 = _mm_mul_pd(_mm_set1_pd(b1), x1);
 m3 = _mm_mul_pd(_mm_set1_pd(b2), x2);
 m4 = _mm_mul_pd(_mm_set1_pd(a1), y1);
 m5 = _mm_mul_pd(_mm_set1_pd(a2), y2);
 m6 = _mm_add_pd(m2, m3);
 m7 = _mm_add_pd(m4, m5);
 m8 = _mm_add_pd(m6, m7);
 m8 = _mm_add_pd(m8, m1);
 */

 // create vector from double-b0-coeff, multiply with input, store:
 __m128d tmp = _mm_mul_pd(_mm_set1_pd(b0), inVector);

 // create vector from double-b1-coeff, multiply, accumulate:
 tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(b1), x1)); 

 // and again and again and again:
 tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(b2), x2));
 tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(a1), y1));
 tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(a2), y2));

 // this is a long dependence-chain - may it can be optimized....

 // update the buffer-variables:
 x2 = x1;
 x1 = inVector;
 y2 = y1;
 y1 = tmp;

 return tmp;
}

#endif // LowpassHighpassStereoSSE2_h
