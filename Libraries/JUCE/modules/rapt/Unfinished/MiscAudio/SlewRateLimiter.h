#ifndef RAPT_SLEWRATELIMITER_H
#define RAPT_SLEWRATELIMITER_H

/** This is a slewrate limiter with user adjustable attack and release time constants. It smoothes
the incoming signal via an attack-release (AR) averager. The averaging time is determined by the
attack and release time constants which are set up in milliseconds. */

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
  void setSampleRate(double newSampleRate);

  /** Sets the attack-time (in milliseconds) - this time which it takes to rise to
  1-1/e = 0.63 when the input signal makes an upward step from 0 to 1. */
  void setAttackTime(double newAttackTime);

  /** Sets the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
  when the input signal makes a downward step from 1 to 0. */
  void setReleaseTime(double newReleaseTime);


  /** \name Inquiry */

  /** Returns the attack-time (in milliseconds) - this time which it takes to rise to
  1-1/e = 0.63 when the input signal makes an upward step from 0 to 1. */
  double getAttackTime() const;

  /** Returns the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
  when the input signal makes a downward step from 1 to 0. */
  double getReleaseTime() const;


  /** \name Audio Processing */

  /** Smoothes the input value vith the AR-averager. */
  double getSample(double in);


  /** \name Misc */

  void reset();

protected:

  /** \name Internal */

  /** Calculates the attack coefficient. */
  void calculateAttackCoefficient();

  /** Calculates the release coefficient. */
  void calculateReleaseCoefficient();


  /** \name Data */

  double y1;                          // previous output sample
  double coeffAttack, coeffRelease;   // attack and release filter coefficients
  double sampleRate;                  // the samplerate
  double attackTime, releaseTime;     // in milliseconds

};

//-----------------------------------------------------------------------------------------------
// inlined functions:




//=================================================================================================

class rsSlewRateLimiterWithHold : public rsSlewRateLimiter
{

public:

  void setSampleRate(double newSampleRate);

  void setHoldTime(double newHoldTime);

  double getSample(double in); // override;

  void reset();

protected:

  double holdTime = 0;
  int holdSamples = 0;
  int sampleCounter = 0;

};




#endif
