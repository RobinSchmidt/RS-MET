#ifndef rosic_DynamicsProcessorBase_h
#define rosic_DynamicsProcessorBase_h

// rosic-indcludes:
#include "../infrastructure/rosic_MutexLock.h"
#include "../others/rosic_SlewRateLimiter.h"
#include "../analysis/rosic_LevelDetector.h"
#include <new> 

#define DEFAULT_LOOKAHEAD_SIZE 8192  // default size for the lookahead buffer

namespace rosic
{

  /**

  This class serves as base-class for various dynamics processors, thereby consolidating their
  common functionality such as the envelope follower, lookahead, etc.

  */

  class DynamicsProcessorBase 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    DynamicsProcessorBase(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~DynamicsProcessorBase();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Chooses the mode for the envelope detector. @see EnvelopeFollower::detectorModes */
    //void setEnvelopeDetectionMode(int newMode) { envelopeFollower.setMode(newMode); }

    /** Sets the attack time of the envelope-follower in milliseconds. */
    void setAttackTime(double newAttackTime);

    /** Sets the release time of the envelope-follower in milliseconds. */
    void setReleaseTime(double newReleaseTime);

    /** Sets the lookahead-time in milliseconds. */
    void setLookAheadTime(double newLookAheadTime);

    /** Sets up a gain factor for the input signal (in dB) .*/
    void setInputGain(double newInputGain) { inputGainFactor = dB2amp(newInputGain); }

    /** Sets up a gain factor for the output signal (in dB) .*/
    void setOutputGain(double newOutputGain) { outputGainFactor = dB2amp(newOutputGain); }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { wet = newDryWet; dry = 1.0-wet; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the lookahead-time in milliseconds. */
    double getLookAheadTime() const { return lookAheadTime; }

    /** Returns the ratio between dry and wet between 0...1. */
    double getDryWetRatio() const { return wet; }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the content of the lookAheadBuffer to zeros and resets the envelope-followers 
    state. */
    void reset();

    //=============================================================================================

  protected:

    /** Sets up the lookhead buffer. */
    void setupLookAhead();

    /** Returns in the output slots the delayed input signal. */
    INLINE void applyLookAheadDelay(double *inL, double *inR, double *outL, double *outR);

    double dry, wet;         // gain factors for dry and wet signal
    double inputGainFactor, outputGainFactor;
    int    tapIn, tapOut;
    int    lookAheadBufferLength;           
    double *lookAheadBufferL, *lookAheadBufferR;
    double lookAheadTime;    // in milliseconds 
    double sampleRate;

    LevelDetector   levelDetector;
    SlewRateLimiter attackReleaseEnveloper;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void DynamicsProcessorBase::applyLookAheadDelay(double *inL, double *inR, 
    double *outL, double *outR)
  {
    // write the incoming signal into the lookahead-buffer:
    lookAheadBufferL[tapIn] = *inL;
    lookAheadBufferR[tapIn] = *inR;

    // read out the possibly delayed signals from the lookahead-buffer:
    *outL = lookAheadBufferL[tapOut];
    *outR = lookAheadBufferR[tapOut];

    // increment tap-pointers:
    tapIn  = wrapAround(tapIn+1,  lookAheadBufferLength);
    tapOut = wrapAround(tapOut+1, lookAheadBufferLength);
  }

} // end namespace rosic

#endif // #ifndef rosic_DynamicsProcessorBase_h
