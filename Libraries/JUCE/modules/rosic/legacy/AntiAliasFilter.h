#ifndef AntiAliasFilter_h
#define AntiAliasFilter_h

#include "AudioModule.h"
#include "IirDesigner.h"   // has an embedded IirDesigner object

/** 

This class implements a direct form IIR filter with the following
characteristics:

mode:          lowpass
approximation: Butterworth
order:         4
cutoff:        18 kHz by default, adjustable via setCutoff()

It uses the IiRDesigner class to calculate its coefficients. This filter
should be used before decimating an oversampled signal.

*/

class AntiAliasFilter : public AudioModule
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

	         AntiAliasFilter();  ///< Constructor.
	virtual ~AntiAliasFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);  
 /**< This is the sample-rate at which the filter operates.
      That means: the target sample-rate times the oversampling-factor */

 virtual void setCutoff(double newCutoff);
 ///< The cutoff frequency of the filter, 18kHz by default.

 //others:
 virtual void reset();
 ///< Resets the internal buffers (for past input and output samples).

 //---------------------------------------------------------------------------
 //audio processing:

 virtual INLINE double getSample(double in); 
 /**< Generates one output sample (we are still at the oversampled 
      sample-rate here). */

 //===========================================================================

protected:

 //embedded objects:
 IirDesigner filtDesigner; ///< The embedded IirDesigner object.

 //parameters:
 doubleA cutoff;            ///< The cutoff frequency.

 //filter coefficients:
 doubleA b[5];   ///< The 5 feedforward coefficients b0,...,b4
 doubleA a[5];   ///< The 5 feedback coefficients a0,...,a4 
                /**< (a[0] is actually not used in the filtering process). */

 //buffering:
 doubleA x[5];   ///< The current and past input samples (=x[n]...x[n-4]).
 doubleA y[5];   ///< The current and past output samples (=y[n]...y[n-4])
                /**< where the y[n] sample represents the current output sample. */
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
INLINE double AntiAliasFilter::getSample(double in)
{
  //set up the x[n]...x[n-4] buffer:
  x[4] = x[3];
  x[3] = x[2];
  x[2] = x[1];
  x[1] = x[0];
  x[0] = in;

  //apply the filter:
  y[0] =     b[0] * x[0]
           + b[1] * x[1]
           + b[2] * x[2]
           + b[3] * x[3]
           + b[4] * x[4]
           - a[1] * y[1]
           - a[2] * y[2]
           - a[3] * y[3]
           - a[4] * y[4];

  //set up the y[n-1]...y[n-4] buffer for the next call:
  y[4] = y[3];
  y[3] = y[2];
  y[2] = y[1];
  y[1] = y[0];

  //return the current output-sample:
  return y[0];
}

#endif // AntiAliasFilter_h