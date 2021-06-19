#ifndef rosic_SlewRateLimiterLinear_h
#define rosic_SlewRateLimiterLinear_h

namespace rosic
{

/** This is a slewrate limiter with user adjustable attack and release time constants. 

ToDo:
-templatize and move to rapt, thereby inline a lot of stuff
-implement a branchless version of rsClip and use that */

class SlewRateLimiterLinear
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SlewRateLimiterLinear();

  /** Destructor. */
  ~SlewRateLimiterLinear();

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the attack-time (in milliseconds) - this time which it takes to rise from 0 to 1
  when the input signal makes an upward step from 0 to 1. */
  void setAttackTime(double newAttackTime);

  /** Sets the release-time (in milliseconds). ... */
  void setReleaseTime(double newReleaseTime);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the attack-time (in milliseconds). */
  double getAttackTime() const { return attackTime; }

  /** Returns the release-time. */
  double getReleaseTime() const { return releaseTime; }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Returns a smoothed input value. */
  double getSample(double in);

  //---------------------------------------------------------------------------------------------
  // others:

  void reset();

  //=============================================================================================

protected:

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






} // end namespace rosic

#endif // rosic_SlewRateLimiterLinear
