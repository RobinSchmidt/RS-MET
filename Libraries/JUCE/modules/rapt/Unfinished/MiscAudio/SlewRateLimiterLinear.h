#ifndef RS_SLEWRATELIMITERLINEAR_H
#define RS_SLEWRATELIMITERLINEAR_H

namespace RSLib
{

  /**

  This is a slewrate limiter with user adjustable attack and release time constants. ....

  */

  class RSLib_API rsSlewRateLimiterLinear  
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsSlewRateLimiterLinear();  

    /** Destructor. */
    ~rsSlewRateLimiterLinear();  


    /** \name Setup */

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);    

    /** Sets the attack-time (in milliseconds) - this time which it takes to rise from 0 to 1
    when the input signal makes an upward step from 0 to 1. */
    void setAttackTime(double newAttackTime); 

    /** Sets the release-time (in milliseconds). ... */
    void setReleaseTime(double newReleaseTime);


    /** \name Inquiry */

    /** Returns the attack-time (in milliseconds). */
    double getAttackTime() const { return attackTime; }

    /** Returns the release-time. */
    double getReleaseTime() const { return releaseTime; }


    /** \name Audio Processing */

    /** Returns a smoothed input value. */
    RS_INLINE double getSample(double in);


    /** \name Misc */

    void reset();

  protected:
    
    /** \name Internal */

    /** Calculates the maximum difference between current and past output sample, when the current
    input is greater than the past output. */
    void calculateUpwardLimit();

    /** Calculates the maximum difference between current and past output sample, when the current
    input is smaller than the past output. */
    void calculateDownwardLimit();


    double calculateStepLimit(double unitStepTime);

    double upwardLimit, downwardLimit;  
    double y1;                          // previous output sample
    double sampleRate;                  // the samplerate
    double attackTime, releaseTime;     // in milliseconds

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsSlewRateLimiterLinear::getSample(double in)
  {
    y1 += rsClipToRange(in-y1, -downwardLimit, upwardLimit);
    return y1;
  }

}

#endif
