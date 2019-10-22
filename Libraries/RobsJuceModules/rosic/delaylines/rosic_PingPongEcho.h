#ifndef rosic_PingPongEcho_h
#define rosic_PingPongEcho_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_MutexLock.h"
//#include "../filters/rosic_Equalizer.h"
//#include "../filters/rosic_LowpassHighpass.h"
////#include "../effects/rosic_StereoPan.h"
//#include <stdio.h>

namespace rosic
{

  /**

  This class implements a stereo delay/echo that bounces back and forth between left and right 
  channel.

  \todo: 
  -introduce diffusor filter (with time and amount parameters)
  -idea: cross-diffusor (kind of 2x2 MIMO allpass)

  */

  class PingPongEcho 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
    PingPongEcho(int maximumDelayInSamples = 1048576);

    /** Destructor */
    ~PingPongEcho();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the delay-time in seconds or beats (depending on whether sync is active). */
    void setDelayTime(double newDelayTime);

    /** Switches the tempo-sync on or off. */
    void setSyncMode(bool shouldTempoSync);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets a global gain value for the wet signal (as raw amplitude factor). */
    //void setGlobalGainFactor(double newFactor) { g = newFactor; }

    /** Sets a global gain value for the wet signal (in decibels). */
    void setWetLevel(double newLevel) { g = RAPT::rsDbToAmp(newLevel); calculateGainFactors(); }
     //

    /** Sets the input gain for the signal that enters the delaylines as raw amplitude. */
    void setDelayInputGain(double newGain) { wetInGain = newGain; }

    /** Sets the feedback factor (as raw amplitude factor). */
    void setFeedbackFactor(double newFactor) { feedback = newFactor; }

    /** Sets the feedback amount in percent. */
    void setFeedbackInPercent(double newPercentage) { feedback = 0.01*newPercentage; }

    /** Sets the panorama position between -1...+1. */
    void setPan(double newPan){ pan = newPan; calculateGainFactors(); }
  
    /** Sets the cutoff frequency of the embedded lowpass filter. */
    void setHighDamp(double newCutoff) 
    { filter1.setLowpassCutoff(newCutoff); filter2.setLowpassCutoff(newCutoff); }

    /** Sets the cutoff frequency of the embedded highpass filter. */
    void setLowDamp(double newCutoff) 
    { filter1.setHighpassCutoff(newCutoff); filter2.setHighpassCutoff(newCutoff); }

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) 
    { dryWetRatio = newDryWet; RAPT::rsEqualPowerGainFactors(dryWetRatio, &dry, &wet, 0.0, 1.0); }

    /** Switches the delayline into/out of ping-pong mode. */
    void setPingPongMode(bool shouldBePingPong) { pingPongMode = shouldBePingPong; }

    /** Switches 'true stereo' mode on/off. */
    void setTrueStereoMode(bool shouldBeTrueStereo) 
    { trueStereoMode = shouldBeTrueStereo; calculateGainFactors(); }

    /** Mutes the output of the delayline. */
    void setMute(bool shouldBeMuted) { mute = shouldBeMuted; }

    /** Bypasses the effect entirely. */
    void setBypass(bool shouldBeBypassed) { bypass = shouldBeBypassed; }

    /** Opposite of setBypass. \todo get rid of one of them */
    void setActive(bool shouldBeActive) { setBypass(!shouldBeActive); }

    /** Switches channel swapping for the wet signal on (relevant only in true stereo mode). */
    void setStereoSwap(bool shouldBeSwapped) { stereoSwap = shouldBeSwapped; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the delay-time in seconds or beats (depending on whether sync is active). */
    double getDelayTime() const { return delayTime; }

    /** Returns true when tempo-sync is active, false otherwise. */
    int isInSyncMode() const { return tempoSync; }

    /** Returns the global gain as value raw amplitude factor. */
    double getGlobalGainFactor() const { return g; }

    /** Returns the feedback factor (as raw amplitude factor). */
    double getFeedbackFactor() const { return feedback; }

    /** Returns the amount of feedback factor in percent. */
    double getFeedbackInPercent() const { return 100.0*feedback; }

    /** Returns the panorama position between -1...+1.. */
    double getPan() const { return pan; }

    /** Returns the ratio between dry and wet between 0...1. */
    double getDryWetRatio() const { return dryWetRatio; }

    /** Returns the cutoff frequency of the embedded lowpass filter. */
    double getHighDamp() const { return filter1.getLowpassCutoff(); }

    /** Returns the cutoff frequency of the embedded highpass filter. */
    double getLowDamp() const { return filter1.getHighpassCutoff(); }

    /** Returns true when the delayline is in ping-pong mode. */
    bool isInPingPongMode() const { return pingPongMode; }

    /** Returns true when the delayline is muted. */
    bool isMuted() const { return mute; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the content of the delayline to zeros and resets the filter's state. */
    void reset();

    //---------------------------------------------------------------------------------------------
    // embedded public modules:

    LowpassHighpass filter1, filter2;

    //=============================================================================================

  protected:

    /** Sets up the delay-time in samples according to the chosen delayTime, sync-mode and 
    sample-rate user parameters. */
    void setupDelayInSamples();

    /** Calculates gLL, gLR, gRL and gRR from the pan-setting (and pingPongMode- and trueStereoMode
    settings). */
    void calculateGainFactors();

    double gLL, gRR, gLR, gRL;            // gain factors for left and right for pan and cross-mix
    double dry, wet;                      // gain factors for dry and wet signal
    int    tapIn, tapOut;                 // write- and read pointers in the delayline
    int    length;                        // delayline length
    double *buffer1, *buffer2;            // delayline buffers
    double delayTime;                     // in seconds or beats
    double sampleRate;                    // the sample-rate
    double bpm;                           // tempo in beats per minute
    double feedback;                      // feedback factor
    double g;                             // global gain factor - renameto wetOutGain or something
    double wetInGain = 1.0;               // input gain factor for the delayline
    double pan;                           // pan user parameter (-1...+1)
    double dryWetRatio;                   // dry/wet ratio 0....1
    bool   tempoSync;                     // flag to indicate tempo-synchronization
    bool   pingPongMode;                  // flag to indicate ping-pong mode
    bool   trueStereoMode;                // flag to indicate true-stereo mode
    bool   mute;                          // flag to indicate that this delayline is muted
    bool   bypass;                        // flag to indicate that the effect is bypassed
    bool   stereoSwap;                    // flag to indicate swapping of the two stereo channels

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void PingPongEcho::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    if( bypass == true )
      return;
      
    double inM, y1, y2, yL, yR, tmp1, tmp2;
    if( trueStereoMode == false )
    {
      inM = wetInGain * SQRT2_INV * (*inOutL + *inOutR);
      if( pingPongMode == false )
      {
        // mono/normal:
        y1             = buffer1[tapOut];
        tmp1           = inM + RAPT::rsClip(feedback*y1, -2.0, 2.0);
        tmp1           = filter1.getSample(tmp1);
        buffer1[tapIn] = tmp1;
        yL             = gLL*y1;
        yR             = gRR*y1;
      }
      else  
      {
        // mono/pingpong:
        y1             = buffer1[tapOut];
        y2             = buffer2[tapOut];
        tmp1           = inM + RAPT::rsClip(feedback*y2, -2.0, 2.0);
        tmp1           = filter1.getSample(tmp1);
        buffer1[tapIn] = tmp1;
        buffer2[tapIn] = feedback*y1;
        yL             = gLL*y1 + gRR*y2;
        yR             = gLL*y2 + gRR*y1;
      }
    } 
    else 
    {
      tmp1 = wetInGain * (gLL*(*inOutL) + gRL*(*inOutR));
      tmp2 = wetInGain * (gLR*(*inOutL) + gRR*(*inOutR));
      if( stereoSwap == true )
        RAPT::rsSwap(tmp1, tmp2);
      yL   = buffer1[tapOut];
      yR   = buffer2[tapOut];

      if( pingPongMode == false )
      {
        tmp1 += RAPT::rsClip(feedback*yL, -2.0, 2.0);
        tmp2 += RAPT::rsClip(feedback*yR, -2.0, 2.0);
      }
      else 
      {
        tmp1 += RAPT::rsClip(feedback*yR, -2.0, 2.0);
        tmp2 += RAPT::rsClip(feedback*yL, -2.0, 2.0);
      }

      tmp1           = filter1.getSample(tmp1);
      tmp2           = filter2.getSample(tmp2);        
      buffer1[tapIn] = tmp1;
      buffer2[tapIn] = tmp2;
    }


    tapIn   = RAPT::rsWrapAround(tapIn+1,  length);
    tapOut  = RAPT::rsWrapAround(tapOut+1, length);  // is not used here, however...
    *inOutL = dry * (*inOutL) + wet * yL;
    *inOutR = dry * (*inOutR) + wet * yR;
  }

} // end namespace rosic

#endif // #ifndef rosic_PingPongEcho_h
