#ifndef rosic_EchoLab_h
#define rosic_EchoLab_h

//// includes from the STL:
//#include <vector>
//using std::vector;
//
//// rosic-indcludes:
//#include "../delaylines/rosic_EchoLabDelayLine.h"
////#include "../delaylines/rosic_ModulatedDelayLine.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This class implements a parallel connection of an arbitrary number of modulated delaylines.

  */

  class EchoLab
  {

  public:

    /** Enumeration of the parameters that exist for each delayline, used internally to handle the thread-safe accesses to the delaylines 
    more compactly. \todo move thread safety out of this class */
    enum delayLineParameters
    {
      DELAY_TIME,
      GAIN_FACTOR,
      FEEDBACK_FACTOR,
      //FEEDFORWARD_FACTOR,
      //BLEND_FACTOR,
      PAN,
      PING_PONG,
      MUTE,
      //DELAY_MOD_CYCLE,
      //DELAY_MOD_DEPTH,
      //DELAY_MOD_PHASE_L,
      //DELAY_MOD_PHASE_R,
      //AMP_MOD_CYCLE,
      //AMP_MOD_DEPTH,
      //AMP_MOD_PHASE_L,
      //AMP_MOD_PHASE_R
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    EchoLab();

    /** Destructor. */
    ~EchoLab();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Copies the data from another instance to bring this instance into the same state with regard to the user parameters. */
    void copyDataFrom(EchoLab *other);

    /** Sets up the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets up the tempo in beats per minute (relevant for tempo sync). */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets up the dry/wet ratio (between 0...1). */
    void setDryWet(double newDryWet);

    /** Sets a level (in dB) for the wet signal. */
    void setWetLevel(double newLevel);

    /** Adds a new delayline into the vector. The integer return-value informs, at which index the new delayline was inserted.  */
    int addDelayLine(double newDelayTime, double newGainFactor);

    /** Removes a delayline from the vector. The return-value informs, if there was actually a delayline removed (if you try to remove a 
    non-existing delayline it will return false */
    bool removeDelayLine(int index);

    /** Removes all delaylines and returns true if there was at least one band to remove. */
    bool removeAllDelayLines();

    /** Sets up a parameter for one of the delaylines with mutex-locking and array bound-checking. The return value informs whether the 
    operation was succesful (i.e. the delayLineIndex was not out of range). @see delayLineParameters */
    bool setDelayLineParameterThreadSafe(int parameterIndex, int delayLineIndex, double newValue);

    /** Changes the delay-time of one of the delaylines (in seconds or beats). */
    bool setDelayTime(int index, double newDelayTime);

    /** Changes the gain of one of the delaylines (as raw amplitude factor). */
    bool setGainFactor(int index, double newGainFactor);

    /** Switches one of the delaylines into or out of solo mode. */
    bool setDelayLineSolo(int index, bool shouldBeSolo);

    /** Switches one of the delaylines into or out of solo mode, pass -1 if none should be solo. */
    bool setDelayLineSolo(int index);

    /** Sets a flag to indicate the the GUI should snap to a time-grid. \todo: move this out of this class. */
    void setSnapToTimeGrid(bool shouldSnap) { snapToTimeGrid = shouldSnap; }

    /** Sets the feedforward factor of one of the delaylines (as raw amplitude factor). */
    //bool setFeedforwardFactor(int index, double newFactor);

    /** Sets the feedback factor of one of the delaylines (as raw amplitude factor). */
    //bool setFeedbackFactor(int index, double newFactor);

    /** Sets the blend factor of one of the delaylines (as raw amplitude factor). */
    //bool setBlendFactor(int index, double newFactor);

    /** Sets the panorama position of one of the delaylines between -1...+1. */
    //bool setPan(int index, double newPan);


    /** Sets the cycle-length/period by which the delay-times are modulated (in seconds or beats) for one of the delaylines. */
    //bool setDelayModulationCycleLength(int index, double newCycleLength);

    /** Sets the depth by which the delay-times are modulated in normlized unit between 0...1 where 1 is the maximum permissible modulation 
    amount (the delay modulates down to zero in this case). \ToDo: write a function that expresses the depth in terms of pitch-deviation in 
    semitones ...has to take frequency into account which means recalculation also for cycle-length and bpm changes
    
    */
    //bool setDelayModulationDepth(int index, double newDepth);

    /** Sets the start phase for the delay-modulation of the left channel (in degrees) for one of the delaylines. */
    //bool setDelayModulationPhaseLeft(int index, double newPhase);

    /** Sets the start phase for the delay-modulation of the right channel (in degrees) for one of the delaylines. */
    //bool setDelayModulationPhaseRight(int index, double newPhase);

    /** Switches tempo-synchronization for the delay-modulation on or off (for all delaylines at once). */
    //void setSyncForDelayModulation(bool shouldBeSynced);

    /** Sets the cycle-length/period by which the amplitudes are modulated (in seconds or beats) for one of the delaylines. */
    //bool setAmplitudeModulationCycleLength(int index, double newCycleLength);

    /** Sets the depth by which the amplitudes are modulated in normlized unit between 0...1 where 1 is the maximum permissible modulation 
    amount (the amplitude modulates down to zero in this case) for one of the delaylines. */
    //bool setAmplitudeModulationDepth(int index, double newDepth);

    /** Sets the start phase for the amp-modulation of the left channel (in degrees) for one of the delaylines. */
    //bool setAmplitudeModulationPhaseLeft(int index, double newPhase);

    /** Sets the start phase for the amp-modulation of the right channel (in degrees) for one of the delaylines. */
    //bool setAmplitudeModulationPhaseRight(int index, double newPhase);

    /** Switches tempo-synchronization for the amplitude-modulation on or off (for all dleaylines at once). */
    //void setSyncForAmplitudeModulation(bool shouldBeSynced);



    /** Switches tempo-synchronization for the delay-times on or off (all at once). */
    void setSyncForDelayTimes(bool shouldBeSynced);

    // pan, mod, feedback, .......

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of delaylines. */
    int getNumDelayLines();

    /** Retruns the dry/wet ratio (between 0...1). */
    double getDryWet() const { return dryWet; }

    /** Retruns the level (in dB) for the wet signal. */
    double getWetLevel() const { return wetLevel; }

    /** Returns, whether or not the delaytimes are tempo-synchronized. */
    bool isDelayTimeSynced() const { return delayTimeSync; }

    /** Requests a parameter of one of the delaylines with mutex-locking and array bound-checking and returns the value.  
    @see delayLineParameters*/
    double getDelayLineParameterThreadSafe(int parameterIndex, int delayLineIndex);

    /** Returns the delaytime of one of the delaylines (0.0 if index is out of range). */
    double getDelayTime(int index);

    /** Returns the gain of one of the delaylines in dB (0.0 if index is out of range). */
    //double getGain(int index);

    /** Returns the gain of one of the delaylines as value raw amplitude factor (0.0 if index is out of range). */
    double getGainFactor(int index);

    /** Returns the index of the delayline which is currently in solo mode (if any), -1 if none. */
    int getSoloedDelayLineIndex() const { return soloedIndex; }

    /** Returns true when the delayline in question is not muted and no other delayline is in solo mode. */
    bool isDelayLineActive(int index);

    bool isSnappingToTimeGrid() const { return snapToTimeGrid; }

    /** Returns the feedforward factor of one of the delaylines as value raw amplitude factor (0.0 if index is out of range). */
    //double getFeedforwardFactor(int index);

    /** Returns the feedback factor of one of the delaylines as value raw amplitude factor (0.0 if index is out of range). */
    //double getFeedbackFactor(int index);

    /** Returns the blend factor of one of the delaylines as value raw amplitude factor (0.0 if index is out of range). */
    //double getBlendFactor(int index);

    /** Returns the panorama position between -1...+1 of one of the delaylines (0.0 if index is out of range). */
    //double getPan(int index);



    /** Returns a pointer to the delayline with the given index. Be careful to invalidate this pointer before you remove the delayline via 
    a call to removeDelayLine or removeAllDelayLines - the object to which the pointer points will be destroyed by such a call. You shall 
    never dereference such a pointer after removing the delayline with that index. Note also that the returned pointer may be NULL if the 
    index is out of range. */
    EchoLabDelayLine* getDelayLine(int index);

    /** Returns a pointer to the equalizer that is embedded in the delayline with the given index. Be careful to invalidate this pointer 
    before you remove the delayline via a call to removeDelayLine or removeAllDelayLines - the object to which the pointer points will be 
    destroyed by such a call. You shall never dereference such a pointer after removing the delayline with that index. Note also that the 
    returned pointer may be NULL if the index is out of range. */
    EqualizerStereo* getFeedbackEqualizer(int index);


    /** Returns a pointer to the equalizer in the input path. @see: getFeedbackEqualizer. */
    EqualizerStereo* getInputEqualizer(int index);


    // pan, mod, feedback, .......

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo output sampleframe. */
    INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);

    /** Calculates an entire block of samples at once. */
    INLINE void processBlock(double* inOutL, double* inOutR, int numSamples);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Aquires the mutex-lock for accessing the vector of the delaylines. */
    void acquireLock() { mutex.lock(); }

    /** Releases the mutex-lock for accessing the vector of the delaylines. */
    void releaseLock() { mutex.unlock(); }

    /** Resets the internal states of all the delaylines filters. */
    void resetDelayLines();

    //=====================================================================================================================================

  protected:

    /** Calculates a stereo output sampleframe without acquiring the mutex-lock and without acquiring the mutex locks for embedded objects. 
    Always wrap calls to that function into acquireLock()/releaseLock() calls and calls that also ensure the mutex for the emebedded 
    objects. */
    INLINE void getSampleFrameStereoWithoutLocks(double *inOutL, double *inOutR);

    std::vector<EchoLabDelayLine*> delayLines;   // vector of the delaylines (we use pointers to 
      // elements because the delayline class does not have a public copy-constructor or 
      // assignement operator - we have to take care about creating and deleting the dlealyline
      // objects ourselves with new and delete
       

    double    sampleRate;                     // sample-rate
    MutexLock mutex;                          // mutex-lock for accessing the array of delaylines

    double bpm, dryWet, wetLevel, dryFactor, wetFactor;
    int    soloedIndex;
    bool   delayTimeSync, delayModulationSync, ampModulationSync, snapToTimeGrid;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void EchoLab::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    mutex.lock();
    getSampleFrameStereoWithoutLocks(inOutL, inOutR);
    mutex.unlock();
  }

  INLINE void EchoLab::getSampleFrameStereoWithoutLocks(double *inOutL, double *inOutR)
  {
    double wetL = 0.0;
    double wetR = 0.0;
    double tmpL = 0.0;
    double tmpR = 0.0;
    double inL  = *inOutL;
    double inR  = *inOutR;

    if( soloedIndex >= 0 && soloedIndex < (int) delayLines.size() )
    {
      //delayLines[soloedIndex]->getSampleFrameStereo(inL, inR, &tmpL, &tmpR, false);
      delayLines[soloedIndex]->getSampleFrameStereoWithoutLocks(&inL, &inR, &tmpL, &tmpR, false);
      wetL = tmpL;
      wetR = tmpR;
    }
    else
    {
      for(unsigned int i=0; i<delayLines.size(); i++)
      {
        //delayLines[i]->getSampleFrameStereo(inL, inR, &tmpL, &tmpR, false);
        delayLines[i]->getSampleFrameStereoWithoutLocks(&inL, &inR, &tmpL, &tmpR, false);
        wetL += tmpL;
        wetR += tmpR;
      }
    }

    *inOutL = inL * dryFactor  + wetL * wetFactor;
    *inOutR = inR * dryFactor  + wetR * wetFactor;
  }

  INLINE void EchoLab::processBlock(double* inOutL, double* inOutR, int numSamples)
  {
    unsigned int i;

    mutex.lock();
    for(i = 0; i < delayLines.size(); i++)
      delayLines[i]->acquireFilterLocks();

    for(int n = 0; n < numSamples; n++)
      getSampleFrameStereoWithoutLocks(&inOutL[n], &inOutR[n]);

    /*
    double tmpL, tmpR;
    for(int n=0; n<numSamples; n++)
    {
      tmpL = (double) inOutL[n];
      tmpR = (double) inOutR[n];
      getSampleFrameStereoWithoutLocks(&tmpL, &tmpR);
      inOutL[n] = (float) tmpL;
      inOutR[n] = (float) tmpR;
    }
    */

    for(i = 0; i < delayLines.size(); i++)
      delayLines[i]->releaseFilterLocks();
    mutex.unlock();
  }

}

#endif 
