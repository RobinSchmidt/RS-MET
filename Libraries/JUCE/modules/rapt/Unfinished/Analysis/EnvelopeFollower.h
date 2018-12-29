#ifndef RAPT_ENVELOPEFOLLOWER_H
#define RAPT_ENVELOPEFOLLOWER_H

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
class rsEnvelopeFollower : public rsSlewRateLimiter
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
  return rsSlewRateLimiter::getSample(in);
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

#endif
