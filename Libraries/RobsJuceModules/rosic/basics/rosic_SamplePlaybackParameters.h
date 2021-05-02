#ifndef rosic_SamplePlaybackParameters_h
#define rosic_SamplePlaybackParameters_h

//// standard includes:
//#include <stdlib.h>
//#include <string>
//
//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class can be used to store and manage meta-data for audio samples which will control the way
  in which the sample is played back in sample based instruments.

  */

  class SamplePlaybackParameters  
  {

  public:

    /** This enumeration lists the available sample-categories. */
    /*
    enum sampleCategories
    {
      ONE_SHOT = 0,
      SINGLE_CYCLE,
      MULTI_CYCLE
    };
    */

    /** This enumeration lists the available loop-modes. */
    enum loopModes
    {
      NO_LOOP = 0,
      FORWARD_LOOP,
      FORWARD_BACKWARD_LOOP
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SamplePlaybackParameters();

    /** Destructor */
    ~SamplePlaybackParameters();

    //---------------------------------------------------------------------------------------------
    // functions for setting up the info:

    /** Sets the total number of samples in the sample-file. */
    void setNumSamples(int newNumSamples) { numSamples = newNumSamples; ensureDataValidity(); }

    /** Sets the number of channels in the sample-file. */
    void setNumChannels(int newNumChannels) { numChannels = newNumChannels; ensureDataValidity(); }

    /** Sets the flag to indicate that the sample is muted. */
    void setMute(bool shouldBeMuted) { mute = shouldBeMuted; }

    /** Sets the flag to indicate that the sample should play solo. */
    void setSolo(bool shouldBeSolo) { solo = shouldBeSolo; }

    /** Sets the sample at which playback should be started for this sample when it is 
    triggered. */
    void setPlaybackStart(double newPlaybackStart)
    { playbackStart = newPlaybackStart; ensureLoopValidity(); }

    /** Sets the key dependence  of the palyback start point. */
    void setPlaybackStartByKey(double newPlaybackStartByKey) 
    { playbackStartByKey = newPlaybackStartByKey; }

    /** Sets the velocity dependence  of the palyback start point. */
    void setPlaybackStartByVel(double newPlaybackStartByVel)
    { playbackStartByVel = newPlaybackStartByVel; }

    /** Sets the sample at which playback should end. */
    void setPlaybackEnd(double newPlaybackEnd) 
    { playbackEnd = newPlaybackEnd; ensureLoopValidity(); }

    /** Turns looping on or off. */
    void setLoopMode(int newLoopMode){ loopMode = newLoopMode; ensureLoopValidity(); }

    /** Sets the start-point (in samples) of the loop (if any). */
    void setLoopStart(double newLoopStart) { loopStart = newLoopStart; ensureLoopValidity(); }

    /** Sets the end-point (in samples) of the loop (if any). */
    void setLoopEnd(double newLoopEnd) { loopEnd = newLoopEnd; ensureLoopValidity(); }

    /** Sets the length of the loop (by adjusting the loopEnd). */
    void setLoopLength(double newLength) 
    {  
      if( newLength >= 1.0 )
        loopEnd = loopStart + newLength;
      ensureLoopValidity();
    }

    /** Locks the length of the loop such that changen loopStart will also affect loopEnd 
    and vice versa. */
    void setLoopLengthLock(bool shouldBeLocked) { loopLengthLock = shouldBeLocked; }

    /** Sets the crossfade length for the loop. */
    void setLoopCrossfadeLength(double newLoopCrossfadeLength) 
    { loopCrossfadeLength = newLoopCrossfadeLength; ensureLoopValidity(); }

    /** The number of pitch cycles contained in the loop - this info is useful for determining the 
    fundamental frequency from the loop-length. A zero value indicates that it is unknown or 
    irrelevant. */
    void setNumPitchCyclesInLoop(double newNumPitchCyclesInLoop) 
    { pitchCyclesInLoop = newNumPitchCyclesInLoop; ensureLoopValidity(); }

    /** Sets the root-key of this sample which is supposed to be the note which was played when
    the sample was recorded. */
    void setRootKey(double newRootKey)
    {
      if( newRootKey >= 0.0 && newRootKey <= 127.0 )
        rootKey = newRootKey;
    }

    /** Sets the detuning of the actual fundamental frequency from the root-key in cents. */
    //void setRootDetune(double newRootDetune) { rootDetune = newRootDetune; }

    /** Sets the lowest key for which this sample should be played (as MIDI note). */
    void setLoKey(int newLoKey) { loKey = newLoKey; ensureKeyRangeValidity(); }

    /** Sets the highest key for which this sample should be played (as MIDI note). */
    void setHiKey(int newHiKey) { hiKey = newHiKey; ensureKeyRangeValidity(); }

    /** Sets the key up to which the sample is to be faded in. At this key itself, it will be 
    fully faded in - can be used for key-crossfades. If this is equal to lowKey, no key fade-in 
    will take place. */
    void setKeyFadeIn(int newKeyToFadeInTo) 
    { fadeInToKey = newKeyToFadeInTo; ensureKeyRangeValidity(); }

    /** Sets the key up from which the sample is to be faded out. At this key itself, it will 
    still be fully faded in - can be used for key-crossfades. If this is equal to highKey, no key 
    fade-out will take place. */
    void setKeyFadeOut(int newKeyToFadeOutFrom)
    { fadeOutFromKey = newKeyToFadeOutFrom; ensureKeyRangeValidity(); }

    /** The lowest velocity for which this sample should be played (from 1...127). */
    void setLoVel(int newLoVel) { loVel = newLoVel; ensureVelRangeValidity(); }

    /** The highest velocity for which this sample should be played (from 1...127). */
    void setHiVel(int newHiVel) { hiVel = newHiVel; ensureVelRangeValidity(); }

    /** Sets the velocity up to which the sample is to be faded in. At this vel itself, it will be 
    fully faded in - can be used for vel-crossfades. If this is equal to lowVel, no vel fade-in 
    will take place. */
    void setVelFadeIn(int newVelToFadeInTo) 
    { fadeInToVel = newVelToFadeInTo; ensureVelRangeValidity(); }

    /** Sets the velocity up from which the sample is to be faded out. At this vel itself, it will 
    still be fully faded in - can be used for vel-crossfades. If this is equal to highVel, no vel 
    fade-out will take place. */
    void setVelFadeOut(int newVelToFadeOutFrom)
    { fadeOutFromVel = newVelToFadeOutFrom; ensureVelRangeValidity(); }

    /** Sets the playback tuning in semitones. */
    void setTune(double newTune) { tune = newTune; }

    /** Sets the keytracking for the playback pitch in cents/key - a value of 100.0 means normal 
    chromatic playback. */
    void setTuneByKey(double newTuneByKey) { tuneByKey = newTuneByKey; }

    /** Sets the velocity tracking for the playback pitch in cents/step. */
    void setTuneByVel(double newTuneByVel) { tuneByVel = newTuneByVel; }

    /** Sets an absolute detuning in Hz. */
    void setDetuneHz(double newDetune) { detuneHz = newDetune; }

    /** Sets the playback level in decibels. */
    void setLevel(double newLevel) { level = newLevel; }

    /** Sets the keytracking for the playback level in percent. */
    void setLevelByKey(double newLevelByKey) { levelByKey = newLevelByKey; }

    /** Sets the velocity tracking for the playback level in percent. */
    void setLevelByVel(double newLevelByVel) { levelByVel = newLevelByVel; }

    /** Sets the sample-rate of the original data - this is the sample rate at which the file was 
    recorded and has nothing to do with the playback sample-rate. */
    void setRecordingSampleRate(double newRecordingSampleRate) 
    {
      if( newRecordingSampleRate > 0.01 )
        recordingSampleRate = newRecordingSampleRate;
    }

    /** Switches the lowpass-/highpass filter for this sample on/off. */
    void setFilterOnOff(bool shouldBeOn) { filterOnOff = shouldBeOn; }

    /** Sets the cutoff frequency for the lowpass filter. */
    void setLowpassCutoff(double newCutoff) 
    { 
      if( newCutoff > 0.0 )
        lowpassCutoff = newCutoff; 
    }

    /** Sets the cutoff frequency for the highpass filter. */
    void setHighpassCutoff(double newCutoff) 
    { 
      if( newCutoff > 0.0 )
        highpassCutoff = newCutoff; 
    }

    // todo: key- and vel-tracking of the cutoff frequencies


    /** Sets the fundamental frequency of the sample. It is supposed to correspond to the rootkey
    but doesnt have to. Set it to zero for non-pitched sounds. */
    //void setFundamentalFrequency(double newFundamentalFrequency);

    /** Sets the name (i.e. the relative path) of the current sample) - this functions is only
    for conviniently handling preset management in a plugIn-context */
    void setSampleName(char* newSampleName);

    //---------------------------------------------------------------------------------------------
    // functions for retrieving the info:

    /** @see setNumSamples() */
    int getNumSamples() const { return numSamples; }

    /** @see setNumChannels() */
    int getNumChannels() const { return numChannels; }

    bool isMuted() const { return mute; }

    bool isSolo() const { return solo; }

    /** @see setPlaybackStart() */
    double getPlaybackStart() const { return playbackStart; }

    /** @see setPlaybackStartByKey() */
    double getPlaybackStartByKey() const { return playbackStartByKey; }

    /** @see setPlaybackStartByVel() */
    double getPlaybackStartByVel() const { return playbackStartByVel; }

    /** @see setPlaybackEnd() */
    double getPlaybackEnd() const { return playbackEnd; }

    /** @see setLoopMode() */
    int getLoopMode() const { return loopMode; }

    /** @see setLoopStart() */
    double getLoopStart() const { return loopStart; }

    /** @see setLoopEnd() */
    double getLoopEnd() const { return loopEnd; }

    /** @see setLoopLength() */
    double getLoopLength() const { return loopEnd-loopStart; }

    /** @see setLoopLengthLock() */
    bool isLoopLengthLocked() const { return loopLengthLock; }

    /** @see setLoopCrossfadeLength() */
    double getLoopCrossfadeLength() const { return loopCrossfadeLength; }

    /** @see setCyclesInLoop() */
    double getNumPitchCyclesInLoop() const { return pitchCyclesInLoop; }

    /** @see setRootKey() */
    double getRootKey() const { return rootKey; }

    /** @see setRootDetune() */
    //double getRootDetune() const { return rootDetune; }

    /** @see setLoKey() */
    int getLoKey() const { return loKey; }

    /** @see setHiKey() */
    int getHiKey() const { return hiKey; }

    /** @see setLoVel() */
    int getLoVel() const { return loVel; }

    /** @see setHiVel() */
    int getHiVel() const { return hiVel; }

    /** @see setFadeInToKey() */
    int getFadeInToKey() const { return fadeInToKey; }

    /** @see setFadeOutFromKey() */
    int getFadeOutFromKey() const { return fadeOutFromKey; }

    /** @see setFadeInToVel() */
    int getFadeInToVel() const { return fadeInToVel; }

    /** @see setFadeOutFromVel() */
    int getFadeOutFromVel() const { return fadeOutFromVel; }

    double getTune()       const { return tune; }
    double getTuneByKey()  const { return tuneByKey; }
    double getTuneByVel()  const { return tuneByVel; }
    double getDetuneHz()   const { return detuneHz; }


    double getLevel()       const { return level; }
    double getLevelByKey()  const { return levelByKey; }
    double getLevelByVel()  const { return levelByVel; }

    /** @see setRecordingSampleRate() */
    double getRecordingSampleRate() const { return recordingSampleRate; }

    /** @see setFundamentalFrequency() */
    //double getFundamentalFrequency() const { return pitchToFreq(rootKey+0.01*rootDetune); }
    double getFundamentalFrequency() const { return RAPT::rsPitchToFreq(rootKey); }

    /** @see setSampleName() */
    char* getSampleName() const { return sampleName; }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Initializes all settings. */
    void initSettings();

    /** Initializes the loop settings to the loop the whole sample. */
    void initLoopSettings();

    /** Initializes the key- and velocity-range and related crossfade- and tracking settings. */
    void initKeyAndVelRanges();

    /** Initializes the level settings. */
    void initLevelSettings();

    /** Initializes the pitch settings. */
    void initPitchSettings();

    /** Initializes the pan settings. */
    void initPanSettings();

    /** Initializes the filter settings. */
    void initFilterSettings();

    //=============================================================================================

  protected:

    /** This function makes sure that all parameters are in meaningful ranges. */
    void ensureDataValidity();

    /** This function makes sure that all loop parameters are in valid ranges. */
    void ensureLoopValidity();

    /** This function makes sure that all keymapping parameters are in valid ranges. */
    void ensureKeyRangeValidity();

    /** This function makes sure that all velocity-mapping parameters are in valid ranges. */
    void ensureVelRangeValidity();

    /** \todo: sort declarations for fast access to frequently used values (loopStart, loopEnd, 
        detuneHz, fundamentalFrequency) */

    // size info:
    int numSamples, numChannels;

    bool mute, solo;

    // loop info:
    int loopMode;
    double pitchCyclesInLoop, playbackStart, playbackStartByVel, playbackStartByKey, playbackEnd, 
      loopStart, loopEnd, loopCrossfadeLength, loopCrossfadeShape;
    bool loopLengthLock;

    // key- and velocity mapping info:
    int loKey, hiKey, loVel, hiVel, fadeInToKey, fadeOutFromKey, fadeInToVel, fadeOutFromVel;
    double rootKey; //, rootDetune;

    // level parameters:
    double level, levelByKey, levelByVel; 
    
    // level ramping parameters:
    double levelStart, levelTime, levelRelease;

    // pitch parameters:
    double tune, detuneHz, tuneByKey, tuneByVel;

    // stereo panorama parameters:
    double pan, panByKey, panByVel, panCurve, midSide;

    // filtering parameters:
    bool   filterOnOff;
    double lowpassCutoff, lowpassByKey, lowpassByVel, highpassCutoff, highpassByKey, highpassByVel;

    // some characteristic parameters of the audio data itself:
    double recordingSampleRate; //, fundamentalFrequency;

    // name of the sample:
    char* sampleName;
    //int   sampleNameLength;  // length of the name - do we actually need this?

  };

} // end namespace rosic

#endif // #ifndef rosic_SamplePlaybackParameters_h
