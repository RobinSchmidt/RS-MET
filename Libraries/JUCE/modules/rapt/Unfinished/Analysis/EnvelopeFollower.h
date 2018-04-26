#ifndef RS_ENVELOPEFOLLOWER_H
#define RS_ENVELOPEFOLLOWER_H

namespace RSLib
{

  /**

  This is an envelope follower with user adjustable attack and release time constants. It squares 
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
  http://www.embedded.com/design/configurable-systems/4006520/Improve-your-root-mean-calculations

  */

  class RSLib_API rsEnvelopeFollower : public rsSlewRateLimiter
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
    rsEnvelopeFollower();  

    /** Destructor. */
    ~rsEnvelopeFollower();  


    /** \name Setup */

    /** Chooses the mode for the envelope detector. @see detectorModes */
    void setMode(int newMode);       


    /** \name Audio Processing */

    /** Smoothes the input value vith the AR-averager. */
    RS_INLINE double applySmoothing(double in);

    /** Estimates the signal envelope via one of the functions getSampleMeanAbsolute(), 
    getSampleMeanSquare() or getSampleRootMeanSquare() depending on the chosen mode. */
    RS_INLINE double getSample(double in);

    /** Estimates the signal envelope via AR-averaging the mean absolute value. */
    RS_INLINE double getSampleMeanAbsolute(double in);

    /** Estimates the signal envelope via AR-averaging the mean squared value. */
    RS_INLINE double getSampleMeanSquare(double in);

    /** Estimates the signal envelope via AR-averaging the mean squared value and extracting the 
    square root. */
    RS_INLINE double getSampleRootMeanSquare(double in);

  protected:

    /** \name Data */

    int mode;  /** @see detectorModes */

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsEnvelopeFollower::applySmoothing(double in)
  {
    return rsSlewRateLimiter::getSample(in);
  }

  RS_INLINE double rsEnvelopeFollower::getSample(double in)
  {
    switch( mode )
    {
    case MEAN_ABS:         return getSampleMeanAbsolute(in);
    case MEAN_SQUARE:      return getSampleMeanSquare(in);
    case ROOT_MEAN_SQUARE: return getSampleRootMeanSquare(in);
    default:               return getSampleMeanAbsolute(in);
    }
  }

  RS_INLINE double rsEnvelopeFollower::getSampleMeanAbsolute(double in)
  {
    return( applySmoothing(fabs(in)) );
  }

  RS_INLINE double rsEnvelopeFollower::getSampleMeanSquare(double in)
  {
    return( applySmoothing(in*in) );
  }

  RS_INLINE double rsEnvelopeFollower::getSampleRootMeanSquare(double in)
  {
    return rsSqrt( getSampleMeanSquare(in) );
  }

}

#endif
