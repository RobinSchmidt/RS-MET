#ifndef RAPT_PEAKFINDING_H
#define RAPT_PEAKFINDING_H

/** A class that detects peaks and then decays exponentially. This can be seen as a simplified
version of an envelope follower with zero attack time. So, in effect, it responds immediately to 
any peaks and then drags an exponentially decaying trail from that peak. If additional smaller 
peaks occur under the umbrella of that trail, they will not be seen separately, they will be 
subsumed/shadowed/masked by the larger peak. This class can be useful for distinguishing major, 
relevant peaks from the minor, irrelevant ones that often sit on the flanks of the major 
mountains. */

template<class T>
class rsPeakMasker
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the number of samples it takes, to decay to a specified value after having seen a unit
  impulse. If r is this ratio and d is the number of samples, our multiplier needs to be the
  d-th root of r such that after d successive multiplications, we reach r, if we init with 1. The 
  default value for the ratio is 0.5, in which case d would give the time to decay to 1/2. Another 
  common value is 1/e ...maybe we should use that as default? */
  void setDecaySamples(T numSamples, T targetRatio = 0.5)
  { c = pow(targetRatio, T(1)/numSamples); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes one output sample at a time. */
  T getSample(T x)
  {
    return y = rsMax(y * c, x);
  }

  /** Per-sample computation function for non-equidistant data. */
  T getSample(T x, T dt)
  {
    return y = rsMax(y * rsPow(c, dt), x);
  }

  /** Applies the process running forward through the signal yIn of length N and writes the result 
  into yOut. Can be used in place: the buffers yIn, yOut may point to the same memory location. */
  void applyForward(const T* yIn, T* yOut, int N)
  {
    reset();
    this->y = yOut[0] = yIn[0]; // first output sample is equal to first input sample
    for(int n = 1; n < N; n++)
      yOut[n] = getSample(yIn[n]);
  }

  /** Like applyForward(const T* yIn, T* yOut, int N) but with an explicitly given x-axis which
  does not necessarily have to be equidistant. */
  void applyForward(const T* x, const T* yIn, T* yOut, int N)
  {
    reset();
    this->y = yOut[0] = yIn[0];
    for(int n = 1; n < N; n++)
      yOut[n] = getSample(yIn[n], x[n]-x[n-1]);
  }

  /** Like applyForward but does a backward pass through the data instead. */
  void applyBackward(const T* yIn, T* yOut, int N)
  {
    reset();
    this->y = yOut[N-1] = yIn[N-1];
    for(int n = N-2; n >= 0; n--)
      yOut[n] = getSample(yIn[n]);
  }

  /** Like applyForward but does a backward pass through the data instead. */
  void applyBackward(const T* x, const T* yIn, T* yOut, int N)
  {
    reset();
    this->y = yOut[N-1] = yIn[N-1];
    for(int n = N-2; n >= 0; n--)
      yOut[n] = getSample(yIn[n], x[n+1]-x[n]);
  }

  /** Resets the internal state. */
  void reset() { y = T(0); }

protected:

  T c = T(0);
  T y = T(0);

};

// ToDo:
// -Maybe make a version with hold - could be useful for limiters
// -Maybe make a version with linear instead of exponential decay - maybe the slope should be 
//  adapted in the "if" according to x (be proportional) ...hmm...maybe not
// -I think, it works correctly only if the input signal is always nonnegative - this is the case 
//  for envelopes but maybe it's useful to have a class that works also when the input may go
//  negative. -> Implement unit tests that cover both cases and make them pass
// -There's some peak-finding/picking related code in MiscUnfinished.h/cpp. That should be moved 
//  here someday (when the API and algos have reasonably stabilized)

#endif