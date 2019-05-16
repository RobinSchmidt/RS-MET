#pragma once

// maybe merge with the EnvelopeFollower.h/cpp files in Unfinished - drag files over and add the 
// new class there

/** An advanced envelope follower implementing the following proceesing chain:

-Butterworth lowpass - gets rid of Gibbs ripples, if any (maybe try Bessel - avoid overshoot)
-full wave rectifier - takes absolute value
-slew rate limiter   - preliminary/raw/simple envelope extraction
-min/max smoother    - extracts average of min and max value in some time window
-Bessel lowpass      - gets rid of jaggies

maybe try to apply the anti-Gibbs-ripple lowpass after taking the absolute value

*/

template<class T> // maybe have separate signal and parameter types - or maybe not - min/max smoothing may need them to be the same...
class rsEnvelopeFollower2
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  // int getDelay()



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


protected:



  rsEngineersFilter<T, T> preFilter;
  rsSlewRateLimiterWithHold<T, T> slewLimiter; // maybe hold is not needed bcs of min/max smoothing
  rsMinMaxSmoother<T> minMaxSmoother;
  rsEngineersFilter<T, T> postFilter;

};