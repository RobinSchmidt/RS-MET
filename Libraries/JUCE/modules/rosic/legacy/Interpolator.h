#ifndef Interpolator_h
#define Interpolator_h

//#include "AudioModule.h"
#include "MoreMath.h"
using namespace MoreMath;

/**

This class implements some interpolation formulas to reconstruct a continous 
signal from samples - that is: it allows for reading out the signal at 
arbitrary non-integer time instants. Many of the implemented interpolation
formulas are taken from the paper "Polynomial Interpolators for High-Quality
Resampling of Oversampled Audio" by Olli Niemitalo.

*/

class Interpolator
{

public:

 //---------------------------------------------------------------------------
 //construction/destruction:

	         Interpolator();
	virtual ~Interpolator();

 //---------------------------------------------------------------------------
 //audio processing:

 INLINE double getSampleLinear(double x, double *y);
 /**< Calclulates one output sample by linearly interpolating between y[0] and
      y[1] at position x. y should point to an array of at least two values, x
      should be a number between 0...1. */

 INLINE double getSampleParabolic4p2o2x(double x, double *y);
 /**< Calclulates one output sample by interpolating parabolically between 4 
      points. Warning: this interpolator has a zero at Nyquist-frequency which
      makes it not so suitable fo non-oversampled audio. */

 INLINE double getSampleHermite4p3o(double x, double *y);
 /**< Calclulates one output sample by interpolating between 4 points with a 
      3rd order Hermite polynomial at position x. The position x should be
      the fractional part of the actual time-instant, the y[1] should be the
      sample value at the integer part of that time-instant such that the 
      *y-array contains two sample-values before and two samples values after
      the time instant, at which the signal is to be read out. */

 INLINE double getSample2ndOsculating4p5o(double x, double *y);
 /**< Calclulates one output sample by interpolating between 4 points with a 
      2nd-order-osculating 5th order polynomial. */

 INLINE double getSampleOptimal4p3o2x(double x, double *y);
 /**< Calclulates one output sample by interpolating between 4 points with a 
      3rd order polynomial, optimized for 2x oversampling. */

 INLINE double getSampleOptimal6p5o2x(double x, double *y);
 /**< Calclulates one output sample by interpolating between 6 points with a 
      5th order polynomial, optimized for 2x oversampling. The position x 
      should be the fractional part of the actual time-instant, the y[2]
      should be the sample value at the integer part of that time-instant such
      that the *y-array contains 3 sample-values before and 3 samples values
      after the time instant, at which the signal is to be read out. */


 //---------------------------------------------------------------------------
 //others:
 void reset();

 //===========================================================================

protected:

 doubleA previousOutput;    // remembers the previous output sample for the 
                           // allpass-interpolator

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double Interpolator::getSampleLinear(double x, double *y)
{
 return y[0] + x*(y[1]-y[0]);
}

INLINE double Interpolator::getSampleParabolic4p2o2x(double x, double *y)
{
 //doubleA y1mym1, c0, c1, c2;

 // 4-point, 2nd-order parabolic 2x (x-form)
 double y1mym1 = y[2]-y[0];
 double c0 = 1/2.0*y[1] + 1/4.0*(y[0]+y[2]);
 double c1 = 1/2.0*y1mym1;
 double c2 = 1/4.0*(y[3]-y[1]-y1mym1);
 return (c2*x+c1)*x+c0;
}

INLINE double Interpolator::getSampleHermite4p3o(double x, double *y)
{
 static doubleA c0, c1, c2, c3;

 // 4-point, 3rd-order Hermite (x-form)
 c0 = y[1];
 c1 = (1.0/2.0)*(y[2]-y[0]);
 c2 = (y[0] - (5.0/2.0)*y[1]) + (2.0*y[2] - (1.0/2.0)*y[3]);
 c3 = (1.0/2.0)*(y[3]-y[0]) + (3.0/2.0)*(y[1]-y[2]);
 return ((c3*x+c2)*x+c1)*x+c0;
}

INLINE double Interpolator::getSample2ndOsculating4p5o(double x, double *y)
{
 static doubleA z, even1, odd1, even2, odd2, c0, c1, c2, c5;

 // 4-point, 5th-order 2nd-order-osculating (z-form)
 z     = x - 1/2.0;
 even1 = y[0]+y[3];
 odd1  = y[0]-y[3];
 even2 = y[1]+y[2]; 
 odd2  = y[1]-y[2];
 c0    = 9/16.0*even2 - 1/16.0*even1;
 c1    = 3/16.0*odd1 - 25/16.0*odd2;
 c2    = 1/4.0*(even1-even2);
 c5    = odd1 - 3*odd2;
 return (((c5*z*z-c5)*z+c2)*z+c1)*z+c0;
}

INLINE double Interpolator::getSampleOptimal4p3o2x(double x, double *y)
{
 static doubleA z, even1, odd1, even2, odd2, c0, c1, c2, c3;

 // Optimal 2x (4-point, 3rd-order) (z-form)
 z     = x - 1/2.0;
 even1 = y[2]+y[1]; 
 odd1  = y[2]-y[1];
 even2 = y[3]+y[0];
 odd2  = y[3]-y[0];
 c0    = even1*0.45868970870461956 + even2*0.04131401926395584;
 c1    = odd1*0.48068024766578432 + odd2*0.17577925564495955;
 c2    = even1*-0.246185007019907091 + even2*0.24614027139700284;
 c3    = odd1*-0.36030925263849456 + odd2*0.10174985775982505;
 return ((c3*z+c2)*z+c1)*z+c0;
}


INLINE double Interpolator::getSampleOptimal6p5o2x(double x, double *y)
{
 doubleA z, 
         even1, odd1, even2, odd2, even3, odd3, 
         c0, c1, c2, c3, c4, c5;

 // Optimal 2x (6-point, 5th-order) (z-form)
 z     = x - 1/2.0;
 even1 = y[3]+y[2]; 
 odd1  = y[3]-y[2];
 even2 = y[4]+y[1];
 odd2  = y[4]-y[1];
 even3 = y[5]+y[0]; 
 odd3  = y[5]-y[0];
 c0    = even1*0.40513396007145713 + 
         even2*0.09251794438424393 + 
         even3*0.00234806603570670;
 c1    = odd1*0.28342806338906690 + 
         odd2*0.21703277024054901 + 
         odd3*0.01309294748731515;
 c2    = even1*-0.191337682540351941 + 
         even2*0.16187844487943592 + 
         even3*0.02946017143111912;
 c3    = odd1*-0.16471626190554542 + 
         odd2*-0.00154547203542499 + 
         odd3*0.03399271444851909;
 c4    = even1*0.03845798729588149 + 
         even2*-0.05712936104242644 + 
         even3*0.01866750929921070;
 c5    = odd1*0.04317950185225609 + 
         odd2*-0.01802814255926417 + 
         odd3*0.00152170021558204;
 return ((((c5*z+c4)*z+c3)*z+c2)*z+c1)*z+c0;
}

#endif // Interpolator_h
