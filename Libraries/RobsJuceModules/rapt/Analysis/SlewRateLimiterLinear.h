#ifndef RAPT_SLEWRATELIMITERLINEAR_H
#define RAPT_SLEWRATELIMITERLINEAR_H

/** This is a slewrate limiter with user adjustable attack and release time constants.... */

template<class TSig, class TPar>
class rsSlewRateLimiterLinear
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsSlewRateLimiterLinear();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the attack-time (in milliseconds) - this time which it takes to rise from 0 to 1
  when the input signal makes an upward step from 0 to 1. */
  void setAttackTime(TPar newAttackTime);

  /** Sets the release-time (in milliseconds). ... */
  void setReleaseTime(TPar newReleaseTime);


  /** \name Inquiry */

  /** Returns the attack-time (in milliseconds). */
  TPar getAttackTime() const { return attackTime; }

  /** Returns the release-time. */
  TPar getReleaseTime() const { return releaseTime; }


  /** \name Audio Processing */

  /** Returns a smoothed input value. */
  RS_INLINE TSig getSample(TSig in);


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


  TPar calculateStepLimit(TPar unitStepTime);

  TPar upwardLimit, downwardLimit;
  TSig y1;                          // previous output sample
  TPar sampleRate;                  // the samplerate
  TPar attackTime, releaseTime;     // in milliseconds

private:

  // make assignment operator and copy constructor unavailable
  rsSlewRateLimiterLinear& operator=(const rsSlewRateLimiterLinear& /*other*/) { return *this; }
  rsSlewRateLimiterLinear(const rsSlewRateLimiterLinear& /*other*/) { }

};

// ToDo: 
// -Factor out a class that doesn't store the sampleRate. It should keep attack and release in 
//  samples
// -Avoid the unary negation in getSample by storing -downwardLimit instead of downwardLimit

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsSlewRateLimiterLinear<TSig, TPar>::getSample(TSig in)
{
  y1 += rsClip(in-y1, -downwardLimit, upwardLimit);
  return y1;
}

#endif
