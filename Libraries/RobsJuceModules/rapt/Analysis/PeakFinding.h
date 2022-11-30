#ifndef RAPT_PEAKFINDING_H
#define RAPT_PEAKFINDING_H

/** A class that detects peaks and then decays exponentially. This can be seen as a simplified
version of an envelope follower with zero attack time. So, in effect, it responds immediately to 
any peaks and then drags an exponentially decaying trail from that peak. If additional smaller 
peaks occur under the umbrella of that trail, they will not be seen separately, they will be 
subsumed/shadowed by the larger peak. This class can be useful for distinguishing major, relevant
peaks from the minor, irrelevant ones. */

template<class T>
class rsPeakTrailDragger  // rename to rsPeakShadower - move to PeakFinding.h/cpp
{

public:

  /** Sets the number of samples it takes, to decay to a specified value after having seen a unit
  impulse. If r is this ratio and d is the number of samples, our multiplier needs to be the
  d-th root of r such that after d successive multiplications, we reach r, if we init with 1. The 
  default value for the ratio is 0.5, in which case d would give the time to decay to 1/2. Another 
  common value is 1/e ...maybe we should use that as default? */
  void setDecaySamples(T numSamples, T targetRatio = 0.5)
  { c = pow(targetRatio, T(1)/numSamples); }


  /** Computes one output sample at a time. */
  T getSample(T x)
  {
    y *= c;
    if(x > y)
      y = x;
    return y;
  }
  // maybe have TSig/TPar template parameters and use rsMax instead of the "if" which should take
  // the element-wise maximum in case of SIMD vector types for the signal
  // maybe implement a getSample with additional dt parameter for use with non-equidistant input

  // function for non-equidistant data - not yet tested:
  T getSample(T x, T dt)
  {
    y *= pow(c, dt);  // is this correct? verify!
    if(x > y)
      y = x;
    return y;
  }


  /** Applies the process running forward through the signal x of length N with time-stamps given 
  in t and writes the result into y. Can be used in place - the buffers x,y may point to the same 
  memory location. */
  void applyForward(const T* t, const T* x, T* y, int N)
  {
    reset();
    this->y = y[0] = x[0]; // first output sample is equal to first input sample
    for(int n = 1; n < N; n++)
      y[n] = getSample(x[n], t[n]-t[n-1]);
  }
  // not yet tested


  void applyBackward(const T* t, const T* x, T* y, int N)
  {
    reset();
    this->y = y[N-1] = x[N-1]; // first output sample is equal to first input sample
    for(int n = N-2; n >= 0; n--)
      y[n] = getSample(x[n], t[n+1]-t[n]);
  }
  // not yet tested

  /** Resets the internal state. */
  void reset() { y = T(0); }

protected:

  T c = T(0);
  T y = T(0);

};
// -maybe make a version with hold - could be useful for limiters
// -maybe make a version with linear instead of exponential decay - maybe the slope should be 
//  adapted in the "if" according to x (be proportional) ...hmm...maybe not
// -i think, it works correctly only if the input signal is always nonnegative - this is the case 
//  for envelopes but maybe it's useful to have a class that works also when the input may go
//  negative

// ToDo:
// -There's some peak-finding/picking related code in MiscUnfinished.h/cpp. That should be moved 
//  here someday (when the API and algos have reasonably stabilized)

#endif