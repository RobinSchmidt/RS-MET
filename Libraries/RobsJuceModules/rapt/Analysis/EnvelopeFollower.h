#ifndef RAPT_ENVELOPEFOLLOWER_H
#define RAPT_ENVELOPEFOLLOWER_H


/** A class that detects peaks and then decays exponentially. This can be seen as a simplified
version of an envelope follower with zero attack time. So, in effect, it responds immediately to 
any peaks and then drags an exponentially decaying trail from that peak. If additional smaller 
peaks occur under teh umbrella of that trail, they will not be seen separately, they will be 
subsumed/shadowed by the larger peak. This class can be useful for distinguishing major/relevant
peaks for the minor irrelevant ones. */

template<class T>
class rsPeakTrailDragger  // maybe rename to rsPeakFollower/rsPeakMeter
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

  /** Resets the internal state. */
  void reset() { y = T(0); }

protected:

  T c = T(0);
  T y = T(0);

};


//=================================================================================================

/** This is an envelope follower with user adjustable attack and release time constants. It squares
the incoming signal or takes the absolute value, smoothes this signal via an attack-release (AR)
averager and (possibly) extracts the square-root. So the output value represents the
instantaneous mean-abs, mean-squared or RMS-value of the incoming signal. The averaging time is
determined by the attack and release time constants which are set up in milliseconds.

References:
(1) DAFX, p. 84

\ToDo: maybe instead of hard-switching between the attack and decay times ta, td, make a
continuous transition based on the relation of the input and output. For example, the time
constant t could be a weighted geometric average between decay and attack time-constant:
t = td^k * ta^(1-k) where k = 0.5 + 0.5*(out-in)/(out+in)
with this formula: if
in  == 0:  k = 1   -> t = td (we set (out-0)/(out+0)= 1, if out==0)
out == 0:  k = 0   -> t = ta (we set (0-in) /(0+in) =-1, if in ==0)
out == in: k = 0.5 -> t = sqrt(td*ta) (geometric mean)
actually, we just have to check, if in==out, set t = sqrt(td*ta) in this case, otherwise use the
formula. As in and out are both nonegative, in != out implies a positive nonzero denominator.

here's an article about optimizing root-mean calculations:
http://www.embedded.com/design/configurable-systems/4006520/Improve-your-root-mean-calculations */

template<class TSig, class TPar>
class rsEnvelopeFollower : public rsSlewRateLimiter<TSig, TPar>
{

public:

  enum detectorModes
  {
    MEAN_ABS,
    MEAN_SQUARE,
    ROOT_MEAN_SQUARE,

    NUM_MODES
  };


  /** \name Construction/Destruction */

  /** Constructor. */
  //rsEnvelopeFollower();

  /** Destructor. */
  //~rsEnvelopeFollower();


  /** \name Setup */

  /** Chooses the mode for the envelope detector. @see detectorModes */
  void setMode(int Mode)
  {
    if( Mode >= MEAN_ABS && Mode < NUM_MODES )
      mode = Mode;
  }


  /** \name Audio Processing */

  /** Smoothes the input value vith the AR-averager. */
  RS_INLINE TSig applySmoothing(TSig in);

  /** Estimates the signal envelope via one of the functions getSampleMeanAbsolute(),
  getSampleMeanSquare() or getSampleRootMeanSquare() depending on the chosen mode. */
  RS_INLINE TSig getSample(TSig in);

  /** Estimates the signal envelope via AR-averaging the mean absolute value. */
  RS_INLINE TSig getSampleMeanAbsolute(TSig in);

  /** Estimates the signal envelope via AR-averaging the mean squared value. */
  RS_INLINE TSig getSampleMeanSquare(TSig in);

  /** Estimates the signal envelope via AR-averaging the mean squared value and extracting the
  square root. */
  RS_INLINE TSig getSampleRootMeanSquare(TSig in);

protected:

  /** \name Data */

  int mode = MEAN_ABS;  /** @see detectorModes */

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsEnvelopeFollower<TSig, TPar>::applySmoothing(TSig in)
{
  return rsSlewRateLimiter<TSig, TPar>::getSample(in);
}

template<class TSig, class TPar>
RS_INLINE TSig rsEnvelopeFollower<TSig, TPar>::getSample(TSig in)
{
  switch(mode)
  {
  case MEAN_ABS:         return getSampleMeanAbsolute(in);
  case MEAN_SQUARE:      return getSampleMeanSquare(in);
  case ROOT_MEAN_SQUARE: return getSampleRootMeanSquare(in);
  default:               return getSampleMeanAbsolute(in);
  }
}

template<class TSig, class TPar>
RS_INLINE TSig rsEnvelopeFollower<TSig, TPar>::getSampleMeanAbsolute(TSig in)
{
  return(applySmoothing(fabs(in)));
}

template<class TSig, class TPar>
RS_INLINE TSig rsEnvelopeFollower<TSig, TPar>::getSampleMeanSquare(TSig in)
{
  return(applySmoothing(in*in));
}

template<class TSig, class TPar>
RS_INLINE TSig rsEnvelopeFollower<TSig, TPar>::getSampleRootMeanSquare(TSig in)
{
  return rsSqrt(getSampleMeanSquare(in));
}

//=================================================================================================


/** An advanced envelope follower implementing the following proceesing chain:

-Butterworth lowpass - gets rid of Gibbs ripples, if any (maybe try Bessel - avoid overshoot)
-full wave rectifier - takes absolute value
-slew rate limiter   - preliminary/raw/simple envelope extraction
-min/max smoother    - extracts average of min and max value in some time window
-Bessel lowpass      - gets rid of jaggies

todo:
-maybe try to apply the anti-Gibbs-ripple lowpass after taking the absolute value
-maybe have separate signal and parameter types - or maybe not - min/max smoothing may need them to 
 be the same... */

template<class T> 
class rsEnvelopeFollower2
{

public:


  rsEnvelopeFollower2();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the sample rate at which this module operates. */
  void setSampleRate(T newSampleRate);

  /** Sets the (expected) frequency of the input signal. This will set up some internal smoothing
  parameters that are most suitable for that particular input frequency. If your signal is not 
  monophonic, you should probably pass the lowest frequency that you expect in the signal. 

  Warning: this parameter does not yet support realtime modulation because the minMaxSmoother has
  to be resetted when this value changes  */
  void setInputFrequency(T newFreq);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  //int getDelay() const 
  //{
  //  return T(3./2.) * minMaxSmoother.getLength();   // factor 3/2 ad hoc from inspection
  //}


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Produces one envelope sample for one input sample at a time. */
  inline T getSample(T in)
  {
    T tmp = getSamplePreFilteredAbs(in); 
    tmp   = getSampleSlewLimited(tmp);
    tmp   = getSampleMinMaxSmoothed(tmp);
    tmp   = getSamplePostFiltered(tmp);
    return tmp;
    //return T(1.23) * tmp; // factor 1.23 compensates gain loss due to lowpass (from inspection)
  }
  // todo: split the function into 
  // getSamplePreFiltered/getSampleSlewLimited/getSampleMinMaxSmoothed/getSamplePostFiltered
  // so we can inspect the signal at all points in the processing chain for making plots

  // todo: the arbitrary factors 3/2 for the delay and 1.23 for the gain should probably be related
  // to the attack/release times of the slewrate limiter - or maybe these times should also be set 
  // up according to the input frequency - some experimentation may be needed - try also very fast
  // envelopes

  /** Resets the internal state variables to their initial values. */
  void reset();


  // These functions are just for making it possible to pick up the signal at various points within
  // the internal processing chain - for experiments, plots, tests, etc. Client code normally just 
  // needs to call the high-level getSample() function:
  //inline T getSamplePreFilteredAbs(T in) { return rsAbs(in); }  // test - no pre-filter
  inline T getSamplePreFilteredAbs(T in) { return rsAbs(preFilter.getSample(in)); } 
  //inline T getSamplePreFilteredAbs(T in) { return rsAbs( T(1.23) * preFilter.getSample(in)); } // try taking abs before the pre-filter
  //inline T getSamplePreFilteredAbs(T in) { return T(1.23)*preFilter.getSample(rsAbs(in)); } // ...nope! much worse!
  inline T getSampleSlewLimited(T in)    { return slewLimiter.getSample(in);      }
  inline T getSampleMinMaxSmoothed(T in) { return minMaxSmoother.getSample(in);   }
  inline T getSamplePostFiltered(T in)   { return postFilter.getSample(in);       }
  // somehere in this chain, we need to apply the gain of 1.23...
  // ...figure out, if this factor (obtained by eyeballing) just applies to the particular 
  // test-signal and/or settings or if it's generally suitable - in the first case, remove it 
  // ...try different waveforms, frequencies (with and without gibbs ripple), envelope 
  // characteristics, etc.



protected:

  /** Updates the min/max smoother and the post-filter according to input frequency and 
  sample-rate. */
  void updateSmoothingFilters();


  rsEngineersFilter<T, T> preFilter;
  rsSlewRateLimiterWithHold<T, T> slewLimiter; // maybe hold is not needed bcs of min/max smoothing
  rsMinMaxSmoother<T> minMaxSmoother;
  rsEngineersFilter<T, T> postFilter;

  // some arbitrary but reasonable initial values:
  T inputFreq  = T(100); 
  T sampleRate = T(44100);


};



#endif