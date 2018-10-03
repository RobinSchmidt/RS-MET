#ifndef rosic_Modulator_h
#define rosic_Modulator_h

// rosic-indcludes:
#include "rosic_BreakpointModulator.h"
#include "rosic_ExponentialRamp.h"
#include "rosic_LowFrequencyOscillator.h"
#include "rosic_SampleModulator.h"
#include "../others/rosic_SlewRateLimiter.h"

namespace rosic
{

  /** 
  This class augments the BreakpointModulator-class by a SlewRateLimiter with 
  adjustable attack and release time-constants, a fade-in feature via an 
  ExponentialRamp (which is useful for LFO-like usage), a timeScale, 
  timeScaleByKey and timeScaleByVel feature, a depth, depthByKey, depthByVel
  adjustment and an option to switch between additive and multiplicative
  apply-modes.

  */

  class Modulator
  {

  public:

    enum modulationSources
    {
      NO_SOURCE = 0,
      BREAKPOINT_MODULATOR,
      SAMPLE_MODULATOR
    };

    enum applyModes
    {
      ADDITIVE = 0,
      MULTIPLICATIVE,
      MULTIPLY_AND_MIX_WITH_ORIGINAL
    };

    //---------------------------------------------------------------------------
    // construction/destruction/operators:

    Modulator();   
    /**< Standard constructor. Constructs a new Modulator object with the default settings. */

    //Modulator(const Modulator& modulatorToCopy);
    /**< Copy constructor. Constructs a new modulator object based on a existing one. */

    ~Modulator();  
    ///< Destructor.

    void copyDataFrom(const Modulator& source);
    /**< Copies the data (the content of all member variables) from the passed source into this
         instance of Modulator. */

    //---------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setModulationSource(int newModulationSource);
    /**< Chooses the source for the modulating signal. 
         @see modulationSources. */

    void setScaleFactor(double newScaleFactor);
    /**< Sets a factor to scale the overall output-signal value. */

    void setOffset(double newOffset);
    /**< Sets an offset to shift the overall output-signal up or down. */

    void setBeatsPerMinute(double newBpm);
    /**< Sets the tempo to which we sync in beat per minute. */

    void setTimeScale(double newTimeScale);
    /**< Sets the time scale factor. */

    void setTimeScaleByKey(double newTimeScaleByKey);
    /**< Sets the key-dependence of the time scale factor. */

    void setTimeScaleByVel(double newTimeScaleByVel);
    /**< Sets the velocity-dependence of the time scale factor. */

    void setDepth(double newDepth);
    /**< Sets the modulation depth. */

    void setDepthByKey(double newDepthByKey);
    /**< Sets the key-dependence of the modulation depth. */

    void setDepthByVel(double newDepthByVel);
    /**< Sets the velocity-dependence of the modulation depth. */

    void setUpwardSlewRateLimit(double newLimit);
    /**< Sets the slew-rate limit for attack like phases. */

    void setUpwardSlewRateLimitByKey(double newLimitByKey);
    /**< Sets the key dependence of the slew-rate limit for attack like 
    phases. */

    void setUpwardSlewRateLimitByVel(double newLimitByVel);
    /**< Sets the velocity dependence of the slew-rate limit for attack like 
    phases. */

    void setDownwardSlewRateLimit(double newLimit);
    /**< Sets the slew-rate limit for release like phases. */

    void setDownwardSlewRateLimitByKey(double newLimitByKey);
    /**< Sets the key dependence of the slew-rate limit for release like 
         phases. */

    void setDownwardSlewRateLimitByVel(double newLimitByVel);
    /**< Sets the velocity dependence of the slew-rate limit for release like 
         phases. */

    void setFadeInTime(double newFadeInTime);
    /**< Sets the time to rise form zero depth to full depth after noteOn. */

    void setFadeInTimeByKey(double newFadeInTimeByKey);
    /**< Sets the key dependence of the fade-in time. */

    void setFadeInTimeByVel(double newFadeInTimeByVel);
    /**< Sets the velocity dependence of the fade-in time. */

    void switchFadeInOnOff(bool fadeInShouldBeOn);
    /**< ´Switches the fade-in on or off. */

    void setFadeOutTime(double newFadeOutTime);
    /**< Sets the time to decay from full depth to zero depth after noteOff. */

    void setFadeOutTimeByKey(double newFadeOutTimeByKey);
    /**< Sets the key dependence of the fade-out time. */

    void setFadeOutTimeByVel(double newFadeOutTimeByVel);
    /**< Sets the velocity dependence of the fade-out time. */

    void switchFadeOutOnOff(bool fadeOutShouldBeOn);
    /**< Switches the fade-out on or off. */

    void setApplyMode(int newApplyMode);
    /**< Sets the mode in which the modulation is applied (see applyModes). */

    void setLoopMode(int newLoopMode);
    /**< Sets the loop mode. */

    void setLoopStart(int newLoopStart);
    /**< Sets the start index (sample- or breakpoint-index depending on the chosen 
         modulationSource) for the loop. */

    void setLoopEnd(int newLoopEnd);
    /**< Sets the end index (sample- or breakpoint-index depending on the chosen 
         modulationSource) for the loop. */

    void setNumCyclesInLoop(int newtNumberOfCyclesInLoop);
    /**< Sets the number number of cycles which the loop represents. This is needed for the 
         calculation of the readout-frequency. */

    void setSyncMode(bool shouldBeSynced);
    /**< Turns the tempo- synchronization on or off. */

    void switchReleaseOnOff(bool shouldBeOn);
    /**< Choose, whether the release phase will be enetred on note-off or not. */

    //---------------------------------------------------------------------------
    // inquiry:

    double getSampleRate();
    /**< Returns the current sampleRate. */

    int getModulationSource() const;
    /**< Returns the currently chosen source for the modulating signal. 
         @see modulationSources. */

    double getScaleFactor() const;
    /**< Returns the signal scale factor. */

    double getOffset() const;
    /**< Returns the signal offset. */

    double getTimeScale() const;
    /**< Returns the time scale factor. */

    double getTimeScaleByKey() const;
    /**< Returns the key-dependence of the time scale factor. */

    double getTimeScaleByVel() const;
    /**< Returns the velocity-dependence of the time scale factor. */

    double getDepth() const;
    /**< Returns the modulation depth. */

    double getDepthByKey() const;
    /**< Returns the key-dependence of the modulation depth. */

    double getDepthByVel() const;
    /**< Returns the velocity-dependence of the modulation depth. */

    double getUpwardSlewRateLimit() const;
    /**< Returns the slew-rate limit for attack like phases. */

    double getUpwardSlewRateLimitByKey() const;
    /**< Returns the key-dependence of the slew-rate limit for attack 
    like phases. */

    double getUpwardSlewRateLimitByVel() const;
    /**< Returns the velocity-dependence of the slew-rate limit for attack 
    like phases. */

    double getDownwardSlewRateLimit() const;
    /**< Returns the slew-rate limit for release like phases. */

    double getDownwardSlewRateLimitByKey() const;
    /**< Returns the key-dependence of the slew-rate limit for release 
    like phases. */

    double getDownwardSlewRateLimitByVel() const;
    /**< Returns the velocity-dependence of the slew-rate limit for release 
    like phases. */

    double getFadeInTime() const;
    /**< Returns the time to rise from zero depth to full depth after noteOn. */

    double getFadeInTimeByKey() const;
    /**< Returns the key-dependence of the fade-in time. */

    double getFadeInTimeByVel() const;
    /**< Returns the velocity-dependence of the fade-in time. */

    bool isFadeInOn();
    /**< Returns true, if fade-in is on, false otherwise. */

    double getFadeOutTime() const;
    /**< Returns the time to decay from full depth to zero depth after noteOff.*/

    double getFadeOutTimeByKey() const;
    /**< Returns the key-dependence of the fade-out time. */

    double getFadeOutTimeByVel() const;
    /**< Returns the velocity-dependence of the fade-out time. */

    bool isFadeOutOn();
    /**< Returns true, if fade-out is on, false otherwise. */

    int getApplyMode();
    /**< Returns the mode in which the modulation is applied (see applyModes).*/

    int getLoopMode();
    /**< Returns the loop mode. */

    int getLoopStart();
    /**< Returns the start index (sample- or breakpoint-index depending on the chosen 
         modulationSource) for the loop. */

    int getLoopEnd();
    /**< Returns the end index (sample- or breakpoint-index depending on the chosen 
         modulationSource) for the loop. */

    int getNumCyclesInLoop();
    /**< Returns the number of cycles contained in the loop. */

    bool isInSyncMode();
    /**< Informs, whether the modulator is in tempo-sync mode or not. */

    bool isReleaseOn();
    /**< Informs, whether the release phase will be enetred off note-on or not.*/

    bool endIsReached();
    /**< Informs wheter the end of the modulation signal has been reached. */

    //---------------------------------------------------------------------------
    // event handling:

    void noteOn(bool startFromCurrentLevel = false, 
      int newKey = 64, int newVelocity = 64);  
    /**< Causes the envelope to start with its attack-phase. When the parameter
    "startFromCurrentLevel" is true, the internal state will not be reset 
    to startAmp, such that the curve begins at the level, where the envelope
    currently is. */  

    void noteOff(bool startFromCurrentLevel = true);  
    /**< Causes the envelope to start with its release-phase. */

    //---------------------------------------------------------------------------
    // audio processing:

    INLINE double getSample();    
    /**< Calculates one output sample at a time. */

    INLINE double applyModulation(double in, bool applyMultiplicatively = false);    
    /**< Calculates one output sample at a time and applies it to the input 
    sample where "apply" means either "add to" (the default) or 
    "multiply by" (set the applyMultiplicatively-flag to true for this). 
    Note that multiplicative apply is much more expensive than additive - 
    not because of the multiplication itself, but because the 
    depth-regulation needs to call a pow()-function in this case.    */

    INLINE void getSampleFrameStereo(double *inL, double *inR, 
      double *outL, double *outR, 
      double *modValue);
    /**< Calculates one output sample-frame at a time and additionally returns 
    the value of the modulator. */

    //---------------------------------------------------------------------------------------------
    // embedded public modules:

    BreakpointModulator breakpointModulator;
    SampleModulator     sampleModulator;

  protected:

    // embedded objects:
    SlewRateLimiter upDownSlewRateLimiter;
    SlewRateLimiter fadeInOutSlewRateLimiter;
    ///ExponentialRamp riseRamp;

    int currentKey; /**< The currently active note-key. */
    int currentVel; /**< The current velocity. */


    double scaleFactor;
    double offset;

    double timeScaleNominal; /**< Time scaling without key- and vel-scaling. */
    double timeScaleByKey;   /**< Key dependence of the time scaling. */
    double timeScaleByVel;   /**< Velocity dependence of the time scaling. */

    double depthNominal;  /**< Modulation depth without key- and vel-scaling. */
    double depth;         /**< Modulation depth with key and vel-scaling. */
    double depthByKey;    /**< Key dependence of the modulation depth. */
    double depthByVel;    /**< Velocity dependence of the modulation depth. */

    double slewUp;           /**< Upward slew-rate limit. */
    double slewUpByKey; 
    double slewUpByVel;

    double slewDown;         /**< Downward slew-rate limit. */
    double slewDownByKey;
    double slewDownByVel;

    double fadeInTime;       /**< Fade-in time for modulation depth. */
    double fadeInTimeByKey; 
    double fadeInTimeByVel;

    double fadeOutTime;      /**< Fade-out time for modulation depth. */
    double fadeOutTimeByKey; 
    double fadeOutTimeByVel;

    int    applyMode;

    bool   fadeInIsOn,
      fadeOutIsOn;

    bool   noRelease; 
    /**< when this is true, the modulator will continue running the loop even 
    after note-off. */

    int    modulationSource;

    Modulator(const Modulator& modulatorToCopy);
    /**< Make a copy-constructor unavailable because this class needs deep copying (because of the 
         pointers in the MutexLocks). If you need to create a new object based on an existing one,
         simply create a new object with the default constructor and then use copyDataFrom(). */
    
  }; // end of class Modulator : public BreakpointModulator

  //----------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE double Modulator::getSample()
  {
    double tmp;
    if( modulationSource == BREAKPOINT_MODULATOR )
      tmp = breakpointModulator.getSample();
    else
      tmp = sampleModulator.getSample();


    tmp        = upDownSlewRateLimiter.getSample(tmp);

    if( currentKey != -1 ) // note is on
      tmp *= fadeInOutSlewRateLimiter.getSample(depth);
    else  // note is off
      tmp *= fadeInOutSlewRateLimiter.getSample(0.0);

    return tmp;
  }

  INLINE double Modulator::applyModulation(double in, bool applyMultiplicatively)
  {
    double tmp = breakpointModulator.getSample();
    tmp        = upDownSlewRateLimiter.getSample(tmp);

    double instantaneousDepth = 0.0;
    if( currentKey != -1 ) // note is on
      instantaneousDepth = fadeInOutSlewRateLimiter.getSample(depth);
    else  // note is off
      instantaneousDepth = fadeInOutSlewRateLimiter.getSample(0.0);

    if( applyMultiplicatively == false )
      return in + instantaneousDepth * depth * tmp;
    else
      return in * RAPT::rsPowBipolar(tmp, instantaneousDepth*depth);
  }

  INLINE void Modulator::getSampleFrameStereo(double *inL, double *inR, 
    double *outL, double *outR, 
    double *modValue)
  {
    double tmp = breakpointModulator.getSample();

    tmp        = scaleFactor * tmp + offset;
    tmp        = upDownSlewRateLimiter.getSample(tmp);

    double instantaneousDepth = 0.0;

    if( fadeInIsOn && currentKey != -1 ) 
      // fade in the depth:
      instantaneousDepth = fadeInOutSlewRateLimiter.getSample(depth);
    else if( fadeOutIsOn && currentKey == -1 )
      // fade out the depth:
      instantaneousDepth = fadeInOutSlewRateLimiter.getSample(0.0);
    else 
      instantaneousDepth = depth;

    /*
    if( currentKey != -1 ) // note is on
    instantaneousDepth = fadeInOutSlewRateLimiter.getSample(depth);
    else if( fadeOutIsOn ) // note is off and fadeOut-switch is on
    instantaneousDepth = fadeInOutSlewRateLimiter.getSample(0.0);
    else
    */

    switch( applyMode )
    {
    case ADDITIVE:
      {
        *modValue = instantaneousDepth * tmp;
        *outL     = *inL + *modValue;
        *outR     = *inR + *modValue;
      }
      break;
    case MULTIPLICATIVE:
      {
        *modValue = RAPT::rsPowBipolar(tmp, depth);
        *outL = (*inL) * (*modValue);
        *outR = (*inR) * (*modValue);

        //*outL = (*inL) * powBipolar(tmp, instantaneousDepth);
        //*outR = (*inR) * powBipolar(tmp, instantaneousDepth);
      }
      break;
    case MULTIPLY_AND_MIX_WITH_ORIGINAL:
      {
        double wet = instantaneousDepth;
        double dry = 1.0-instantaneousDepth;
        *modValue  = tmp;
        //*modValue  = (wet*instantaneousDepth) + dry;

        *outL = dry * (*inL) + tmp*wet*(*inL);
        *outR = dry * (*inL) + tmp*wet*(*inR);
      }
      break;
    } // end of switch( applyMode )
  }

}  // end of namespace rosic

#endif // rosic_Modulator_h
