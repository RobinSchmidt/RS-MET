#ifndef RS_SLEWRATELIMITER_H
#define RS_SLEWRATELIMITER_H

namespace RSLib
{

  /**

  This is a slewrate limiter with user adjustable attack and release time constants. It smoothes 
  the incoming signal via an attack-release (AR) averager. The averaging time is determined by the 
  attack and release time constants which are set up in milliseconds.

  */

  class RSLib_API rsSlewRateLimiter  
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsSlewRateLimiter();  

    /** Destructor. */
    ~rsSlewRateLimiter();  


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
    double getAttackTime() const { return attackTime; }

    /** Returns the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
    when the input signal makes a downward step from 1 to 0. */
    double getReleaseTime() const { return releaseTime; }


    /** \name Audio Processing */

    /** Smoothes the input value vith the AR-averager. */
    RS_INLINE double getSample(double in);


    /** \name Misc */

    void reset();

  protected:

    /** \name Internal */

    /** Calculates the attack coefficient. */
    void calculateAttackCoefficient();

    /** Calculates the release coefficient. */
    void calculateReleaseCoefficient();


    /** \name Data */

    double coeffAttack, coeffRelease;   // attack and release filter coefficients
    double y1;                          // previous output sample
    double sampleRate;                  // the samplerate
    double attackTime, releaseTime;     // in milliseconds

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsSlewRateLimiter::getSample(double in)
  {
    if( y1 < in )
      y1 = in + coeffAttack  * (y1-in);
    else
      y1 = in + coeffRelease * (y1-in);
    return y1;
  }



  //-----------------------------------------------------------------------------------------------
  // a class to facilitate the use of stereo-slewrate limiters (a better solution would be to use
  // a template parameter SampleFrame in the class above - if a double is passed as 
  // template-parameter, we would have a mono version, if some kind of "DoublePair" is passed, we 
  // would have stereo and we could even pass multichannel arrays)

  class rsSlewRateLimiterStereo
  {

  public:

    void setSampleRate(double newSampleRate)
    { left.setSampleRate(newSampleRate); right.setSampleRate(newSampleRate); }

    void setAttackTime(double newAttackTime)
    { left.setAttackTime(newAttackTime); right.setAttackTime(newAttackTime); }

    void setReleaseTime(double newReleaseTime)
    { left.setReleaseTime(newReleaseTime); right.setReleaseTime(newReleaseTime); }

    void getSampleFrameStereo(double *inOutL, double *inOutR)
    {
      *inOutL = left.getSample( *inOutL);
      *inOutR = right.getSample(*inOutR);
    }

    void reset() { left.reset(); right.reset(); }

  protected:

    rsSlewRateLimiter left, right;

  };

}

#endif
