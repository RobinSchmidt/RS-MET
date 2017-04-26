#ifndef rosic_LevelDetector_h
#define rosic_LevelDetector_h

// rosic-indcludes:
#include "../filters/rosic_LeakyIntegrator.h"
#include "rosic_InstantaneousEnvelopeDetector.h"
#include "../others/rosic_SlewRateLimiter.h"

namespace rosic
{

  /**

  This is a level detector with user adjustable averaging time constant. It squares 
  the incoming signal or takes the absolute value, smoothes this signal via an averager and 
  (possibly) extracts the square-root. So the output value represents the instantaneous 
  mean-abs, mean-squared or RMS-value of the incoming signal. 

  */

  class LevelDetector 
  {

    friend class NoiseGate;  // needs direct access to the embedded SlewRateLimiter
    friend class Limiter;    // ditto

  public:

    enum detectorModes
    {
      MEAN_ABS,
      MEAN_SQUARE,
      ROOT_MEAN_SQUARE,
      INSTANTANEOUS_ENVELOPE,

      NUM_MODES
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    LevelDetector();  

    /** Destructor. */
    ~LevelDetector();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the attack time constant (in ms) for the attack/release follower. */
    void setAttackTime(double newAttackTime);

    /** Sets the release time constant (in ms) for the attack/release follower. */
    void setReleaseTime(double newReleaseTime);

    /** Chooses the mode for the envelope detector. @see detectorModes */
    void setMode(int  newMode);       

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Estimates the level of the incoming signal (in decibels). */
    INLINE double getLevel(double in);

    /** Estimates the amplitude envelope of the incoming signal (as raw amplitude). */
    INLINE double getAmplitudeEnvelope(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    int mode;   // absolute/squared/RMS
    InstantaneousEnvelopeDetector instEnvDetector;
    LeakyIntegrator               preSmoother;
    SlewRateLimiter               attackReleaseFollower;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double LevelDetector::getLevel(double in)
  {
    const double lowAmplitude = 0.0000000000001; // to avoid log-of-zero in amp->dB conversion
    return amp2dBWithCheck(getAmplitudeEnvelope(in), lowAmplitude);
  }

  INLINE double LevelDetector::getAmplitudeEnvelope(double in)
  {
    double tmp;
    switch( mode )
    {
    case MEAN_ABS:               tmp = fabs(in);                                     break;
    case MEAN_SQUARE:            tmp = in*in;                                        break;
    case ROOT_MEAN_SQUARE:       tmp = in*in;                                        break;
    case INSTANTANEOUS_ENVELOPE: tmp = instEnvDetector.getInstantaneousEnvelope(in); break;
    default:                     tmp = fabs(in);
    }
    tmp = preSmoother.getSample(tmp);
    tmp = attackReleaseFollower.getSample(tmp);
    if( mode == ROOT_MEAN_SQUARE )
      tmp = sqrt(tmp);

    return tmp;
  }

} // end namespace rosic

#endif // rosic_LevelDetector
