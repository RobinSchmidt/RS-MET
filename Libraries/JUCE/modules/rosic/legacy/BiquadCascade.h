#ifndef BiquadCascade_h
#define BiquadCascade_h

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

class BiquadCascade
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

BiquadCascade();  ///< Constructor.     
~BiquadCascade();  ///< Destructor.

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

 INLINE double getSampleDirect1(double in); 
 /**< Calculates a single filtered output-sample via a cascade of biquads in 
      Direct-Form 1 */

 INLINE double getSampleDirect2(double in); 
 /**< Calculates a single filtered output-sample via a cascade of biquads in 
      Direct-Form 2 */


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

 // filter parameters:
 intA numStages;

 // maximum number of biquad-stages:
 static const int maxNumStages = 12;

 // filter coefficients:
 doubleA gain;             // a global gain-factor
 doubleA a0[maxNumStages]; // a0-coefficients for the individual stages
 doubleA a1[maxNumStages]; // a1-coefficients for the individual stages
 doubleA a2[maxNumStages]; // a2-coefficients for the individual stages
 doubleA b0[maxNumStages]; // b0-coefficients for the individual stages
 doubleA b1[maxNumStages]; // b1-coefficients for the individual stages
 doubleA b2[maxNumStages]; // b2-coefficients for the individual stages

 // buffering:
 doubleA x1[maxNumStages]; // arrays of previous biquad input samples, 
 doubleA x2[maxNumStages]; // array index indicates the biquad-stage
 doubleA y1[maxNumStages]; // arrays of previous output samples
 doubleA y2[maxNumStages];

 doubleA g1[maxNumStages]; // arrays of previous intermediate samples for DF2
 doubleA g2[maxNumStages];

 //
 doubleA sampleRate;
 doubleA sampleRateRec;
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE void BiquadCascade::setCoeffs(double *newB0, 
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

INLINE double BiquadCascade::getSampleDirect1(double in)
{
 doubleA tmp, tmp2;
 intA i;  //for the loop through the stages

 tmp = in;

 // calculate current output-sample (y[n]) of all the BiQuad-stages
 // (the output of one stage is the input for the next stage):
 for (i=0; i<numStages; i++)  
 {
  tmp2 = tmp; // for x1[i]

  // calculate current output-sample (y[n]) of BiQuad-stage i:
  tmp = b0[i]*(tmp) + b1[i]*x1[i] + b2[i]*x2[i]
                    - a1[i]*y1[i] - a2[i]*y2[i]
                    + TINY; // to avoid denorm problems

  // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next iteration:
  x2[i] = x1[i];
  x1[i] = tmp2;
  y2[i] = y1[i];
  y1[i] = tmp;  
 }

 return gain*tmp;
}

INLINE double BiquadCascade::getSampleDirect2(double in)
{
 static doubleA x, y, g;
 static intA    i;  // for the loop through the stages

 x = in;

 // calculate current output-sample (y[n]) of all the BiQuad-stages
 // (the output of one stage is the input for the next stage):
 for (i=0; i<numStages; i++)  
 {
  // calculate current output-sample (y[n]) of BiQuad-stage i:
  g = x - a1[i]*g1[i] - a2[i]*g2[i];
  y = b0[i]*g + b1[i]*g1[i] + b2[i]*g2[i];

  // set g[n-1], g[n-2] for the next iteration:
  g2[i] = g1[i];
  g1[i] = g;

  x = y; // output of one stage is input to the next
 }

 return gain*y;
}



#endif // BiquadCascade_h
