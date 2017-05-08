#ifndef rosic_ExponentialSmooother_h
#define rosic_ExponentialSmooother_h

#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

/** A class to smooth an incoming target value. You set up the target value via setTargetValue() and
per sample retrieve the smoothed value via getSample(). 

\todo: maybe rename to OnePoleSmoother, maybe let the time-constant be passed in milliseconds

WARNING: class has not yet been tested */

class ExponentialSmoother
{

public:

  /** Sets the time constant (in seconds) which is the time it takes to reach around 
  1 - 1/e = 63% of the target value when starting from zero. You must also pass the samplerate
  at which the smoother should operate here. */
  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    coeff = exp(-1.0 / (sampleRate * timeConstant));
  }

  /** Sets the target value that should be approached. */
  void setTargetValue(double newTargetValue)
  {
    target = newTargetValue;
  }

  /** Sometimes, you may want to set the internal state of the current value manually from outside. 
  With this function, you can do this. */
  void setCurrentValue(double newCurrentValue)
  {
    current = newCurrentValue;
  }

  /** Returns a smoothed output sample. */
  INLINE double getSample()
  {
    current = target + coeff * (current - target);
    return current;
  }

protected:

  // member variables:
  doubleA current = 0;
  doubleA target  = 0;
  doubleA coeff   = 0;

};

}


#endif