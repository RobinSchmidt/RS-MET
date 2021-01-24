#ifndef romos_VoiceAllocator_h
#define romos_VoiceAllocator_h


/** This class can be used for managing voices in a polyphonic instrument. */

class VoiceAllocator
{

public:

  enum voiceStealingModes
  {
    STEAL_OLDEST_VOICE = 0,
    STEAL_NEWEST_VOICE,
    //STEAL_NEAREST_VOICE,
    //STEAL_LOWEST_VOICE,
    //STEAL_HIGHEST_VOICE
    NUM_VOICE_STEALING_MODES
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  VoiceAllocator();

  /** Destructor. */
  ~VoiceAllocator();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the number of voices - should be between 1 and maxNumVoices. */
  void setNumVoices(int newNumVoices);

  /** Chooses one of the voice-stealing modes - should be one of the numbers in the 
  voiceStealingModes enumeration. */
  void setVoiceStealingMode(int newVoiceStealingMode);

  /** Selects the behaviour when a note that is already playing is triggered again. In retriggering
  mode, the same voice will be re-used, in non-retriggering mode, a new voice will be used. */
  void setRetriggerMode(bool shouldRetrigger);

  /** Sets all our trigger-flags back to "false". This should be the first thing to do in the 
  top-level per-sample (or per block) function. */
  INLINE void resetTriggerFlags()
  {
    noteOnTriggerFlags  = 0;
    noteOffTriggerFlags = 0;
  }

  /** Sets the note-on flag for the given voice to "true". Called from noteOn(). */
  void setNoteOnTriggerFlag(int voiceIndex);

  /** Sets the note-off flag for the given voice to "true". Called from noteOff(). */
  void setNoteOffTriggerFlag(int voiceIndex);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns number of available voices as selected by the user, i.e. the currently selected 
  polyphony. This, in turn, determines the amount of memory that is allocated for the voices' 
  signal buffers. The value is at least 1 and at most maxNumVoices. */
  INLINE int getNumVoices() const { return numVoices; }

  /** Returns number of currently playing voices. */
  INLINE int getNumPlayingVoices() const { return numPlayingVoices; }

  /** Returns an array filled with the indices of the currently playing voices. The array contents 
  are valid up to getNumPlayingVoices()-1, above that value, there may be garbage or an 
  access-violation (the latter occurs, when all voices are playing) - in any case, don't access 
  array-indices above getNumPlayingVoices()-1. */
  INLINE const int* const getPlayingVoiceIndices() const { return playingVoiceIndices; }

  /** Returns true when the voice in question is currently playing a note, false otherwsie. */
  INLINE bool isNoteOn(int voiceIndex) const
  {
    return voiceStates[voiceIndex].normalizedVelocity != 0.0;
  }

  /** Returns whether or not the voice with given index is currently playing - this is different 
  from isNoteOn because voices that are in noteOff-state (velocity == 0) due to calling noteOff(), 
  are still considered as playing in their release phase. A voice will only stop playing after you 
  call killVoice for the voice in quiestion. */
  INLINE bool isVoicePlaying(int voiceIndex) const
  {
    return voiceStates[voiceIndex].isPlaying;
  }

  /** Returns the MIDI key of the most recent note that is/was played in the voice with given 
  index. */
  INLINE int getKeyOfVoice(int voiceIndex) const
  {
    return voiceStates[voiceIndex].key;
  }

  /** Returns the frequency of the most recent note that is/was played in the voice with given 
  index. */
  INLINE double getFrequencyOfVoice(int voiceIndex) const
  {
    return voiceStates[voiceIndex].frequency;
  }

  /** Returns the velocity (as value between 0...1) of the note that is currently played in the 
  voice with given index. A zero velocity indicates that no note is playing. Beware that "no note 
  playing" does not necessarily imply that the voice is silent - there can be notes that are still
  in their (possibly long) release phase. */
  INLINE double getNormalizedVelocityOfVoice(int voiceIndex) const
  {
    return voiceStates[voiceIndex].normalizedVelocity;
  }

  /** Returns the flag that indicates that the voice in question has just received a note-on. It 
  will be set in noteOn() and should be reset by client-code using resetTriggerFlags(). */
  INLINE bool getNoteOnTriggerFlag(int voiceIndex) const
  {
    return ((noteOnTriggerFlags << voiceIndex) & maskForFlags) > 0;
  }

  /** Returns the flag that indicates that the voice in question has just received a note-off. It 
  will be set in noteOff() and should be reset by client-code using resetTriggerFlags(). */
  INLINE bool getNoteOffTriggerFlag(int voiceIndex) const
  {
    return ((noteOffTriggerFlags << voiceIndex) & maskForFlags) > 0;
  }

  /** Returns true if any of our note-on/off trigger flags is set to true. */
  INLINE bool isAnyTriggerFlagSet() const
  {
    return noteOnTriggerFlags != 0 || noteOffTriggerFlags != 0;
  }

  //-----------------------------------------------------------------------------------------------
  // event handling:

  /** Handles a note-on by updating numPlayingVoices and the playingVoiceIndices array. It returns 
  the index of the voice that is to be used for the new note. A call to noteOn with velocity == 0 
  will be interpreted as noteOff (according to the MIDI convention), so in this case, noteOn will 
  call noteOff with the given key and return noteOff's return-value. */
  int noteOn(int key, int velocity);

  /** Handles a note-off by updating numPlayingVoices and the playingVoiceIndices array. Returns 
  the index of the voice that receives (or should receive) the note-off */
  int noteOff(int key);

  /** Kills the voice with given index - that is: sets it into "off" state. This is different from 
  noteOff in that noteOff merely sets the voice's velocity to 0 but doesn't remove it from the 
  playingVoices array. The idea is, that a note-off will cause the instrument to enter the 
  release-phase (by turning all NoteGates to zero, etc.) and at some time later, a VoiceKiller 
  module will actually kill the voice (remove it from the array of playing voices and marking it as
  "off"). */
  void killVoice(int voiceIndex);

  /** Resets all voices into default-state. */
  void reset();

  /** Resets the voice with given index into default state. */
  void resetVoice(int voiceIndex);

  //===============================================================================================

protected:

  /** Returns the index of the first voice that is currently free (i.e. not playing). Returns -1 
  if no voice is free. */
  int getFirstFreeVoice();

  /** Returns the index inside the playingVoiceIndices array of the voice that is to be stolen, 
  depending on the setting of the voice-stealing mode. Some voice-stealing modes determine the 
  voice which is to steal by the incoming note's key, so you need to pass it as parameter */
  int getPlayingIndexToSteal(int key);

  /** Returns the index of the voice that has been playing for the longest time. Returns -1, if no 
  voice is currently playing. */
  int getOldestPlayingVoiceIndex();

  /** Returns the index of the voice that has been playing for the shortest time. Returns -1, if no 
  voice is currently playing. */
  int getNewestPlayingVoiceIndex();

  /** Searches in our currently playing voices the first one which is currently playing a note with
  the given key and returns its index in the playingVoiceIndices array. The search can start at an 
  arbitrary index, so you can search for the 2nd, 3rd etc. voice also after you have found the 1st,
  2nd, etc. The last parameter can be used to ignore voices with zero velocity (i.e. voices that 
  have already received a note-off). */
  int findAmongPlayingVoices(int key, int searchStartIndex = 0, 
    bool ignoreVoicesWithZeroVelocity = false);

  /** Searches in our currently playing voices for the voice with given voiceIndex and returns its 
  position inside this array (the "playingIndex"). If the voice in question is not among the 
  currently playing ones, the function returns -1. */
  int getPlayingIndexOfVoice(int voiceIndex);

  /** Removes the voice with given index from our array of playing voice indices (if it is actually 
  in this array, if not, it does nothing).  */
  void removeFromPlayingVoices(int voiceIndex);

  /** Moves the entry of our playingVoiceIndices-array with the given index to the bottom of the 
  array, thereby moving up all entries that were formerly below it. The entry at the bottom of the 
  array must always be the voice that has been triggered most recently, so this function should be 
  called whenever a voice that was already playing is re-used (either due to retriggering or due to
  voice-stealing)  */
  void movePlayingVoiceToBottom(int playingIndex);

  /** Moves the entries of our playingVoiceIndices-array one position up, starting at the index 
  that is passed as parameter (including that index), For example, if you pass 2, the entry 
  playingVoiceIndices[2] will be propagated to playingVoiceIndices[1], playingVoiceIndices[3] 
  becomes playingVoiceIndices[2], etc. Finally playingVoiceIndices[numPlayingVoices-1] will be 
  propagated to playingVoiceIndices[numPlayingVoices-2]. Called from movePlayingVoiceToBottom and 
  killVoice. */
  void moveUpPlayingVoices(int fromWhichIndex);

  //-----------------------------------------------------------------------------------------------
  // data:

  static const unsigned long maskForFlags = 0x80000000;  // 2^31

  static const int maxNumVoices = 32; 
  // should be sizeof(unsigned long) to be compatible with our trigger-flags

  int  numVoices;
  int  stealingMode;
  bool retriggerMode;
  int  numPlayingVoices;
  int  playingVoiceIndices[maxNumVoices];

  unsigned long noteOnTriggerFlags, noteOffTriggerFlags;
  // 32 flags to indicate that the voice was triggered at the current sample, maybe use a 
  // "Flags32" class here later

  // we need a member tuningTable (or tuningFunction, temperament)

  struct VoiceState
  {
    int    key;                 // key as MIDI note value
    double frequency;           // derived from key and the tuning to be used
    double normalizedVelocity;  // velocity normalized to the range 0...1
    bool   isPlaying;           // flag to indicate that the voice is playing
  };

  VoiceState voiceStates[maxNumVoices];

};

extern VoiceAllocator voiceAllocator;

//}

#endif
