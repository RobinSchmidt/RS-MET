#ifndef rosic_EnvelopeFollower_h
#define rosic_EnvelopeFollower_h

//// rosic-indcludes:
//#include "../others/rosic_SlewRateLimiter.h"

namespace rosic
{

  /**

  This is an envelope follower with user adjustable attack and release time constants. It squares 
  the incoming signal or takes the absolute value, smoothes this signal via an attack-release (AR) 
  averager and (possibly) extracts the square-root. So the output value represents the 
  instantaneous mean-abs, mean-squared or RMS-value of the incoming signal. The averaging time is 
  determined by the attack and release time constants which are set up in milliseconds.

  References: Udo Zoelzer - DAFX (page 84)

  */

  class EnvelopeFollower : public SlewRateLimiter
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
    EnvelopeFollower();  

    /** Destructor. */
    ~EnvelopeFollower();  


    /** \name Setup */

    /** Chooses the mode for the envelope detector. @see detectorModes */
    void setMode(int newMode);       


    /** \name Audio Processing */

    /** Smoothes the input value vith the AR-averager. */
    INLINE double applySmoothing(double in);

    /** Estimates the signal envelope via one of the functions getSampleMeanAbsolute(), 
    getSampleMeanSquare() or getSampleRootMeanSquare() depending on the chosen mode. */
    INLINE double getSample(double in);

    /** Estimates the signal envelope via AR-averaging the mean absolute value. */
    INLINE double getSampleMeanAbsolute(double in);

    /** Estimates the signal envelope via AR-averaging the mean squared value. */
    INLINE double getSampleMeanSquare(double in);

    /** Estimates the signal envelope via AR-averaging the mean squared value and extracting the 
    square root. */
    INLINE double getSampleRootMeanSquare(double in);

  protected:


    /** \name Data */

    int mode;  /** @see detectorModes */

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double EnvelopeFollower::applySmoothing(double in)
  {
    return SlewRateLimiter::getSample(in);
  }

  INLINE double EnvelopeFollower::getSample(double in)
  {
    switch( mode )
    {
    case MEAN_ABS:         return getSampleMeanAbsolute(in);
    case MEAN_SQUARE:      return getSampleMeanSquare(in);
    case ROOT_MEAN_SQUARE: return getSampleRootMeanSquare(in);
    default:               return getSampleMeanAbsolute(in);
    }
  }

  INLINE double EnvelopeFollower::getSampleMeanAbsolute(double in)
  {
    return( applySmoothing(fabs(in)) );
  }

  INLINE double EnvelopeFollower::getSampleMeanSquare(double in)
  {
    return( applySmoothing(in*in) );
  }

  INLINE double EnvelopeFollower::getSampleRootMeanSquare(double in)
  {
    return sqrt( getSampleMeanSquare(in) );
  }

}

#endif
