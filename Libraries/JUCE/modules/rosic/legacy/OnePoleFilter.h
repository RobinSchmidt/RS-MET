#ifndef OnePoleFilter_h
#define OnePoleFilter_h

//#include "AudioModule.h"
#include "Definitions.h"
#include "MoreMath.h"

/**

This is an implementation of a simple one-pole filter unit, which mimics an
analog RC-Filter circuit. It has 2 modes: lowpass and highpass.

*/

class OnePoleFilter
{

public:

 /** This is an enumeration of the available filter modes. */
 enum modes
 {
  BYPASS = 0,
  LOWPASS,
  HIGHPASS
 };

 //---------------------------------------------------------------------------
 // construction/destruction:

 OnePoleFilter();   ///< Constructor.
 ~OnePoleFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setMode(int newMode);
 ///< Chooses the filter mode. See the enumeration above for available modes.


 void setCutoff(double newCutoff);


 void setCoeffs(double newB0, 
                double newB1, 
                double newA1);
 ///< Sets the filter coefficients manually.


 //---------------------------------------------------------------------------
 // audio processing:

 INLINE double getSample(double in);


 //---------------------------------------------------------------------------
 // others:
 void resetBuffers();

 //===========================================================================

protected:

 // buffering:
 double x_1;
 double y_1;

 // filter coefficients:
 double b0; // feedforward coeffs
 double b1;
 double a1; // feedback coeff


 // filter parameters:
 double cutoff;
 int    mode;  

 double sampleRate; 
 double sampleRateRec;  // reciprocal of the sampleRate

 // internal functions:
 void calcCoeffs();  // calculates filter coefficients from filter parameters

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double OnePoleFilter::getSample(double in)
{
 // calculate the output sample:
 y_1 =   b0*in 
       + b1*x_1 
       + a1*y_1
       + TINY;

 // update the buffer variables:
 x_1 = in;

 return y_1;
}

#endif // OnePoleFilter_h
