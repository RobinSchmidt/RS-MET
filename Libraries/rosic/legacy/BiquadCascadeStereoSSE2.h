#ifndef BiquadCascadeStereoSSE2_h
#define BiquadCascadeStereoSSE2_h

//#include "AudioModule.h"
//#include "Definitions.h"
#include "MoreMath.h"
using namespace MoreMath;

/**

This class implements a cascade of up to 12 (= const short maxNumStages) 
biquad-filter stages. The individual stages have their own set of coefficients
which has to be set from outside this class - this class does not do the 
filter-design. The coefficients can be calculated by means of the related 
BiquadDesigner class or the IirDesigner class, for example.

*/

class BiquadCascadeStereoSSE2
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

BiquadCascadeStereoSSE2();  ///< Constructor.     
~BiquadCascadeStereoSSE2();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Set a new sample-rate.

 void setNumStages(int newNumStages);    
 ///< Lets the user set the number of biquad stages.

 // parameters that are supposed to be updated at double-rate (inlined
 // access-functions):

 INLINE void setCoeffs(double *newB0, double *newB1, double *newB2,
                       double *newA1, double *newA2,
                       double newGain = 1.0);
 /**< Allows the user to set the filter coefficients for the individual
      biquad-stages. The difference-equation of each of the biquad stages is:
      y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
                     - a1*y[n-1] - a2*y[n-2]   */


 //---------------------------------------------------------------------------
 // audio processing:

 INLINE __m128d getSampleVector(__m128d inVector); 
 /**< Calculates a single filtered output-vecotr via a cascade of biquads in 
      Direct-Form 1 */

 //---------------------------------------------------------------------------
 // others:

 void initBiquadCoeffs (); 
 /**< Initializes the biquad coefficients to neutral values. */

 void reset(); 
 /**< Sets the buffers for the previous input and output doubles of all biquad
      stages to zero. */

 void getMagnitudeResponse(double *frequencies, 
                           double *magnitudes, 
                           int     numBins, 
                           bool    inDecibels = false, 
                           bool    accumulate = false);
 /** Calculates the magnitudes of the frequency-response at the frequencies
     given in the array "frequencies" (in Hz) and stores them in the 
     array "magnitudes". Both arrays are assumed to be "numBins" long. 
     "inDecibels" indicates, if the frequency response should be returned in
     decibels. If "accumulate" is true, the magnitude response of this 
     biquad-cascade will be multiplied with (or added to, when "inDecibels" 
     is true) to the magnitudes which are already there in the 
     "magnitudes"-array. This is useful for calculating the magnitude response
     of several biquad-cascades in series. */

 void getMagnitudeResponse(float *frequencies, 
                           float *magnitudes, 
                           int    numBins, 
                           bool   inDecibels = false, 
                           bool   accumulate = false);
 /** The same as above but for single-precision floats. */


 //===========================================================================

protected:

 // maximum number of biquad-stages:
 static const int maxNumStages = 12;

 // buffering:
 __m128d x1[maxNumStages]; // arrays of previous biquad input samples, 
 __m128d x2[maxNumStages]; // array index indicates the biquad-stage
 __m128d y1[maxNumStages]; // arrays of previous output samples
 __m128d y2[maxNumStages];

 // filter coefficients:
 doubleA gain;             // a global gain-factor
 doubleA a0[maxNumStages]; // a0-coefficients for the individual stages
 doubleA a1[maxNumStages]; // a1-coefficients for the individual stages
 doubleA a2[maxNumStages]; // a2-coefficients for the individual stages
 doubleA b0[maxNumStages]; // b0-coefficients for the individual stages
 doubleA b1[maxNumStages]; // b1-coefficients for the individual stages
 doubleA b2[maxNumStages]; // b2-coefficients for the individual stages

 // filter parameters:
 intA numStages;

 //
 doubleA sampleRate;
 doubleA sampleRateRec;
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE void BiquadCascadeStereoSSE2::setCoeffs(double *newB0, 
                                     double *newB1, 
                                     double *newB2,                                                
                                     double *newA1, 
                                     double *newA2,
                                     double  newGain)
{
 static intA i;  // for the loop through the stages

 // copy the passed coefficients into the variables:
 for(i=0; i<numStages; i++)
 {
  b0[i] = newB0[i];
  b1[i] = newB1[i];
  b2[i] = newB2[i];
  a1[i] = newA1[i];
  a2[i] = newA2[i];
 }

 gain = newGain;

 return;
}

INLINE __m128d BiquadCascadeStereoSSE2::getSampleVector(__m128d inVector)
{
 __m128d tmp, tmp2;
 intA i;  //for the loop through the stages

 tmp = inVector;

 // calculate current output-sample (y[n]) of all the BiQuad-stages
 // (the output of one stage is the input for the next stage):
 for (i=0; i<numStages; i++)  
 {
  tmp2 = tmp; // for x1[i]

  tmp = _mm_mul_pd(_mm_set1_pd(b0[i]), tmp);

  // create vector from double-b1-coeff, multiply, accumulate:
  tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(b1[i]), x1[i])); 

  // and again and again and again:
  tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(b2[i]), x2[i])); 
  tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(a1[i]), y1[i])); 
  tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_set1_pd(a2[i]), y2[i])); 

  // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next iteration:
  x2[i] = x1[i];
  x1[i] = tmp2;
  y2[i] = y1[i];
  y1[i] = tmp;  

  /*
  // calculate current output-sample (y[n]) of BiQuad-stage i:
  tmp = b0[i]*(tmp) + b1[i]*x1[i] + b2[i]*x2[i]
                    + a1[i]*y1[i] + a2[i]*y2[i]
                    + TINY; // to avoid denorm problems

  // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next iteration:
  x2[i] = x1[i];
  x1[i] = tmp2;
  y2[i] = y1[i];
  y1[i] = tmp;  
  */
 }

 return tmp;
}

#endif // BiquadCascadeStereoSSE2_h
