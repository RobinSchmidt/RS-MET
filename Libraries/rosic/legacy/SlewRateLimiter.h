#ifndef SlewRateLimiter_h
#define SlewRateLimiter_h

/**

This is a slew rate limiter with adjustable attack and release time constants.
It works similar to the "EnvelopeFollower" class, except that it does not work
on rectified or squared signals but on the original bipolar signal. It is
useful to smooth discontinouities in control signals such as an LFO-output 
when the LFO produces a quare wave.

© Braindoc 2005 (www.braindoc.de)

*/

#include "Definitions.h"
#include <math.h>

class SlewRateLimiter  
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

	SlewRateLimiter();
 ~SlewRateLimiter();

 //---------------------------------------------------------------------------
 // parameter settings:

	void setSampleRate (double newSampleRate);
 /**< Sets the sample-rate. */

	void setAttackTime (double newAttackTime); 
 /**< Sets the attack-time in seconds */

	void setReleaseTime(double newReleaseTime);
 /**< Sets the release-time in seconds */

 //---------------------------------------------------------------------------
 // audio processing:

	INLINE double getSample(double in);

 //---------------------------------------------------------------------------
 // others:

	void reset();

protected:

 double y_1;	// previous ouput sample of the unit

	double coeffAttack, coeffRelease; 
	 // the recursion coefficient for the recursive RC-Filter - the smoothing is 
  // performed via the recursive filter equation: 
  // y[n] = (1-coeff)*x[n] + coeff*y[n-1]. 
  // Which coeff of the two will be used is determined by the signal - if it 
  // is rising ( y[n-1] < x[n] ), the attack-coefficient is used, if it is 
  // falling, the release coeeficient will be used


	double sampleRate;
	double attackTime, releaseTime; // time constants in seconds
};

//----------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions
// which are supposed to be called at audio-rate (they can't be put into
// the .cpp file):

double SlewRateLimiter::getSample(double in)
{
	static double temp, coeff;

	temp = in;

	// decide, if signal is rising or falling and choose appropriate filter
 // coefficient:
	if( y_1 < temp )
		coeff = coeffAttack;
	else
		coeff = coeffRelease;

	// perform the envelope detection:
	temp = (1.0-coeff)*temp 
        + coeff*y_1
        + TINY;

	// update the memorized output for next call:
	y_1 = temp;

	return temp;
}
// to be optimized.....


#endif // #ifndef SlewRateLimiter_h
