#ifndef rosic_ExponentialRamp_h
#define rosic_ExponentialRamp_h

// standard library includes:
#include "math.h"

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

 /** 

 This is a simple ramp envelope with exponential characteristic. It has only
 3 parameters: the start value, the end value and the time in which the
 envelope ramps up or down from the start value to the end value. It is
 implemented with a RC-type lowpass-filter, such that the characteristic of the 
 ramp is exponential (technically speaking: it utilizes the step response
 of a simple one-pole lowpass). This also implies, that the end value will not
 really be reached (in theory) but only approached asymptotically. So the
 time value is not really the time between start and end, but rather a time
 constant tau of the RC unit. The time conastant tau is defined as the time
 until the filter reaches 63.2% of the end value (for an incoming
 step-function). This time constant can be scaled via setTauScale() to 
 re-define the ramp time to other values than 63.2%.

 */

 class ExponentialRamp
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  ExponentialRamp();   ///< Constructor.
  ~ExponentialRamp();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< sets the sets sample-rate.

  void setStart(double newStart);      
  ///< Sets the start value of the ramp.

  double getTime() const;
  /**< Retruns the time between start and end in ms (approximately, see class
       description). */

  void setTime(double newTime);
  /**< Sets the time between start and end in ms (approximately, see class
       description). */

  void setEnd(double newEnd);        
  ///< Sets the end value of the ramp.

  void setTauScale(double newTauScale);   
  /**< Scale the time constant tau if you want to reach other levels 
       than 63.2% in "time" */

  //---------------------------------------------------------------------------
  // event processing:

  void trigger();  
  ///< Triggers the ramp envelope.

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSample();

  //===========================================================================

 protected:

  // filter coefficients:
  doubleA inputCoeff;        // weight for the current input-sample
  doubleA recursionCoeff;	   // weight for the previous output-sample    

  // buffering:
  doubleA out;               // current and previous output value

  // parameters:
  doubleA end, time, start;  // parameters of the envelope
  doubleA tauScale;          // scaler for the time constants
  doubleA sampleRate;

 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):
 INLINE double ExponentialRamp::getSample()
 {
  // calculate output-sample:
  out =   inputCoeff*end 
        + recursionCoeff*out
        + TINY;

  // return the output-sample:
  return out;
 }

} // end namespace rosic

#endif // rosic_ExponentialRamp_h_
