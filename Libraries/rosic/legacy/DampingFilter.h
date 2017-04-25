#ifndef DampingFilter_h
#define DampingFilter_h

//#include "AudioModule.h"
#include "Definitions.h"
#include "MoreMath.h"
using namespace MoreMath;

/**

This class implements a filter which is meant be used as a frequency dependent
gain inside a feedback loop, so as to apply a frequency dependent damping 
inside the loop. It applies a global gain factor to a signal as well as a 
first order low-shelving and a first order high-shelving filter. 
In that respect, it is very much like the ToneControl class, but here we use a
different definition of the corner-frequency. In the ToneControl class, the
corner-frequency is defined to be the frequency at which the gain is the
geometric mean between the gain the reference gain (which is unity) - here the 
rereference gain does not need to be unity and the relative gain at the 
corner-frequency can be specified arbitrarily. This facilitates the use of the 
filter inside the feedback-loop of a delay-line - here it may be desirable to 
define the corner-freq in terms of decay-time instead of in terms of gain.

*/

class DampingFilter
{

public:

 //---------------------------------------------------------------------------
 //construction/destruction:

 DampingFilter();   ///< Constructor.
 ~DampingFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setGlobalGainFactor(double newGlobalGainFactor);
 /** Sets up a global gain factor (as factor, not in dB). */

 void setLowGainFactor(double newLowGainFactor);
 /**< Sets the realative gain factor of the low-shelving filter  */

 void setLowCrossoverFreq(double newLowCornerFreq);
 /**< Sets the corner frequency of the low-shelving filter */

 void setLowCrossoverGainFactor(double newLowCrossoverGainFactor);
 /**< Sets the relative gain factor at which the crossover-frequency is 
      measured for the low-shelving filter. See Orfanidis' paper about 
      High Order Equalizers for details. */

 void setHighGainFactor(double newHighGainFactor);
 /**< Sets the relative gain factor of the high-shelving filter. */

 void setHighCrossoverFreq(double newHighCornerFreq);
 /**< Sets the corner frequency of the high-shelving filter. */

 void setHighCrossoverGainFactor(double newHighCrossoverGainFactor);
 /**< Sets the relative gain factor at which the crossover-frequency is 
      measured for the high-shelving filter. */

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 //others:

 void resetBuffers();
 ///< Resets the internal buffers of the filter to zero.

 //===========================================================================

protected:

 // buffering:
 doubleA out;
 doubleA x1, x2;  // past input values
 doubleA y1, y2;  // past output values

 // overall 2-pole filter-coefficients:
 doubleA b0, b1, b2; // feedforward coeffs
 doubleA     a1, a2; // feedback coeffs

 // 1-pole filter-coefficients:
 doubleA b0Ls, b1Ls, a1Ls; // low-shelf coeffs
 doubleA b0Hs, b1Hs, a1Hs; // high-shelf coeffs

 // filter parameters:
 doubleA globalGainFactor;

 doubleA lowCrossoverFreq;
 doubleA lowCrossoverGainFactor;
 doubleA lowGainFactor;

 doubleA highCrossoverFreq;
 doubleA highCrossoverGainFactor;
 doubleA highGainFactor;

 doubleA sampleRate;
 doubleA sampleRateRec;  // reciprocal of the sampleRate
 
 // private member functions:
 void calcLowShelfCoeffs();  // calculates the low-shelving coeffs
 void calcHighShelfCoeffs(); // calculates the high-shelving coeffs
 void calcBiquadCoeffs();    // combines the two first order filters into one
                             // biquad filter

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double DampingFilter::getSample(double in)
{
 // calculate output-sample:
 out = b0*in + b1*x1 + b2*x2
             + a1*y1 + a2*y2;

 // update buffer-variables:
 x2 = x1;
 x1 = in;
 y2 = y1;
 y1 = out;

 return out;
}

#endif // DampingFilter_h
