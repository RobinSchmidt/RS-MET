#ifndef RAPT_SLEWRATELIMITER_H
#define RAPT_SLEWRATELIMITER_H

/** This is a slewrate limiter with user adjustable attack and release time constants. It smoothes
the incoming signal via an attack-release (AR) averager. The averaging time is determined by the
attack and release time constants which are set up in milliseconds. */

template<class TSig, class TPar>
class rsSlewRateLimiter
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsSlewRateLimiter();

  /** Destructor. */
  //~rsSlewRateLimiter();


  /** \name Setup */

  /** Sets the sample-rate for this envelope detector. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the attack-time (in milliseconds) - this time which it takes to rise to
  1-1/e = 0.63 when the input signal makes an upward step from 0 to 1. */
  void setAttackTime(TPar newAttackTime);

  /** Sets the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
  when the input signal makes a downward step from 1 to 0. */
  void setReleaseTime(TPar newReleaseTime);


  /** \name Inquiry */

  /** Returns the attack-time (in milliseconds) - this time which it takes to rise to
  1-1/e = 0.63 when the input signal makes an upward step from 0 to 1. */
  TPar getAttackTime() const { return attackTime; }

  /** Returns the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
  when the input signal makes a downward step from 1 to 0. */
  TPar getReleaseTime() const { return releaseTime; }


  /** \name Audio Processing */

  /** Smoothes the input value vith the AR-averager. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  void reset();

protected:

  /** \name Internal */

  /** Calculates the attack coefficient. */
  void calculateAttackCoefficient();

  /** Calculates the release coefficient. */
  void calculateReleaseCoefficient();


  /** \name Data */

  TSig y1;                          // previous output sample
  TPar coeffAttack, coeffRelease;   // attack and release filter coefficients
  TPar sampleRate;                  // the samplerate
  TPar attackTime, releaseTime;     // in milliseconds

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsSlewRateLimiter<TSig, TPar>::getSample(TSig in)
{
  if(y1 < in)
    y1 = in + coeffAttack  * (y1-in);
  else
    y1 = in + coeffRelease * (y1-in);
  return y1;
}


//=================================================================================================

template<class TSig, class TPar>
class rsSlewRateLimiterWithHold : public rsSlewRateLimiter<TSig, TPar>
{

public:

  void setSampleRate(TPar newSampleRate) //override ..is only for runtime polymorphism
  {
    rsSlewRateLimiter<TSig, TPar>::setSampleRate(newSampleRate);
    holdSamples = rsRoundToInt(holdTime * sampleRate);
  }

  void setHoldTime(TPar newHoldTime)
  {
    holdTime = newHoldTime;
    holdSamples = rsRoundToInt(holdTime * sampleRate);
  }

  RS_INLINE TSig getSample(TSig in); // override;

  void reset() // override
  {
    rsSlewRateLimiter<TSig, TPar>::reset();
    sampleCounter = 0;
  }

protected:

  TPar holdTime = 0;
  int holdSamples = 0;
  int sampleCounter = 0;

};

template<class TSig, class TPar>
RS_INLINE TSig rsSlewRateLimiterWithHold<TSig, TPar>::getSample(TSig in)
{
  if(y1 > in) {
    if(sampleCounter >= holdSamples)
      y1 = in + coeffRelease * (y1-in);
    else
      sampleCounter++;
  }
  else {
    y1 = in + coeffAttack  * (y1-in);
    sampleCounter = 0;
  }
  return y1;
}


#endif
