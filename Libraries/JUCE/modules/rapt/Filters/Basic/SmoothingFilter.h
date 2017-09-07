#ifndef RAPT_SMOOTHINGFILTER_H_INCLUDED
#define RAPT_SMOOTHINGFILTER_H_INCLUDED


/** A filter for smoothing signals in the time domain. It consists of a series connection of first
order lowpass filters (aka leaky integrators). The user can set the order, i.e. the number of 
lowpasses. 

(todo: The filter will scale the time constants according to the order so as to maintain
comparable transition times, so the transition time will not depend on the order. The transition
doesn't get more and more sluggish, when increasing the order (which would happen without the 
automatic downscaling of the time constant) */

// todo: templatize the class
class rsSmoothingFilter
{

public:

  /** Sets the time constant (in seconds) which is the time it takes to reach around 
  1 - 1/e = 63% of the target value when starting from zero. You must also pass the samplerate
  at which the smoother should operate here. */
  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    coeff = exp(-1.0 / (sampleRate * timeConstant));
  }


  /** Returns a smoothed output sample. */
  inline double getSample(double in)
  {
    return y1 = in + coeff*(y1-in); // from rosic::LeakyIntegrator
  }

protected:

  // member variables:
  double y1    = 0;  // y[n-1]
  double coeff = 0;

};

#endif