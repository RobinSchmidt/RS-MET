#ifndef rosic_PolyphonicInstrumentVoice_h
#define rosic_PolyphonicInstrumentVoice_h

namespace rosic
{


/** This class is the base class for a single voice for polyphonic instruments (which should be 
derived from PolyphonicInstrument). These two classes should be used to handle voice management. A 
voice can keep track of more than one incoming note in order to facilitate glides which 
automatically glide back to their origin. For this reason, a voice keeps a list of notes, the front
of this list always represents the currently playing note. When deriving from this class, make sure 
to wrap accesses to the noteList in mutex locks because such accesses can happen from the 
audio-thread as well as from the GUI-thread (preset-switching may cause calls to reset()). A 
member 'mutex' is already in place for this.

\todo make this last thing with the mutex obsolete by moving thread-sync responsibility to 
rosof - we need to be verrryy careful to not miss to insert ScopedLocks everywhere....
\todo a lot of the data defined here can be shared among the voices - make a struct 
InstrumentData or something and let the voices  maintain pointers to it  */

class PolyphonicInstrumentVoice
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  PolyphonicInstrumentVoice();

  /** Destructor. */
  virtual ~PolyphonicInstrumentVoice();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Purely virtual - is supposed to set the sample-rate for this voice. */
  virtual void setSampleRate(double newSampleRate);

  /** Sets the level of the voice in decibels. */
  virtual void setLevel(double newLevel);

  /** Sets the key dependence of the voice level in dB at 127. */
  virtual void setLevelByKey(double newLevelByKey);

  /** Sets the velocity dependence of the voice level in dB at 127. */
  virtual void setLevelByVel(double newLevelByVel);

  /** Sets the tuning frequency for the A4 reference tone (usually 440 Hz). */
  virtual void setMasterTuneA4(double newTuneA4);

  /** Sets the range for the pitch wheel in semitones. */
  virtual void setPitchWheelRange(double newRange);

  /** Sets the glide time in milliseconds. */
  virtual void setGlideTime(double newGlideTime);

  /** Switches glide on or off. */
  virtual void setGlideMode(bool shouldGlide);

  /** This function is empty here in this base-class, you will need to override it in your subclass 
  if you need sync functionality. */
  virtual void setBeatsPerMinute(double newBpm);

  /** Assigns a rosic::TuningTable object to be used for the note-to-frequency conversions, by 
  default, it will use the standard tuning (equally tempered scale). */
  virtual void setTuningTable(TuningTable* newTable);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the level of the voice in decibels. */
  virtual double getLevel();

  /** Returns the key dependence of the voice level in dB at 127. */
  virtual double getLevelByKey();

  /** Returns the velocity dependence of the voice level in dB at 127. */
  virtual double getLevelByVel();

  /** Returns the tuning frequency for the A4 reference tone (usually 440 Hz). */
  virtual double getMasterTuneA4();

  /** Returns the range for the pitch wheel in semitones. */
  virtual double getPitchWheelRange();

  /** Returns the glide time in milliseconds. */
  virtual double getGlideTime();

  /** Informs, whether glide-mode is on or off. */
  virtual bool isInGlideMode();

  /** Informs, whether the voice is currently in a release phase. */
  ///virtual bool isInReleasePhase();

  /** When the voice is in release phase, this function returns the note number beign released, 
  otherwise -1. */
  //virtual bool getNoteBeingReleased();

  /** Returns the key of the currently playing note. */
  virtual int getCurrentNoteKey();

  /** Returns the velocity of the currently playing note. */
  virtual int getCurrentNoteVelocity();

  /** Checks, whether the note in question is in the list of notes. */
  virtual bool hasNoteInList(int noteKeyToCheckFor);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** This function is supposed to calculate the output-samples for both channels and add them at 
  the adresses of *outL and *outR. The third output tells the current amplitude of the voice 
  (i.e. the output of the amp-envelope.) to enable the PolyphonicInstrument class to apply 
  automatic volume scaling acording to the number of  playing voices and their loudnesses. The 
  implementation of this function should generally be done in the respective subclasses - however,
  this base-class here provides a basic skeleton implementation which takes care of incrementing 
  the 'currentNoteAge' and setting up the 'currentFrequency' member variable. */
  virtual INLINE void getSampleFrameStereo(double* outL, double* outR, double* voiceAmplitude);

  //-----------------------------------------------------------------------------------------------
  // event handling:

  /** Triggers a note on this voice. */
  virtual void noteOn(int newKey, int newVelocity, int newDetune = 0);

  /** Triggers a note off this voice. */
  virtual void noteOff(int newKey);

  /** Triggers a note off for all notes in the list. */
  virtual void allNotesOff();

  /** Sets the value of the pitch-wheel in normalized units between -1...+1. */
  virtual void setPitchBend(double newPitchBend);

  /** Sets the voice into its default state (clears the note-list, etc). */
  virtual void reset();

  //-----------------------------------------------------------------------------------------------
  // public member variables (to be accessed by the outlying PolyphonicInstrument class):

  int  currentNoteAge;        // age of the note in samples (for last note priority voice assignment)
  int  noteBeingReleased;     // note being released, otherwise -1
  int  noteBeingReleasedVel;  // velocity of note being released
  bool isSilent;
  bool isReleasing;

  // make them protected and either provide accessors or declare PolyphonicInstrument as friend

  //===============================================================================================

protected:

  /** Returns the frequency of a note with an integer number according to the assigned 
  rosic::TuningTable or according to a standard conversion formula when no tuning table is 
  assigned. */
  virtual double getNoteFrequency(int noteNumber);

  /** Returns the frequency of a note with a non-integer number according to the assigned 
  rosic::TuningTable or according to a standard conversion formula when no tuning table is 
  assigned. */
  virtual double getNoteFrequency(double noteNumber);

  /** Sets up the variables that are responsible for the amplitude glide which takes place along 
  with the pitch glide (mainly to avoid clicks). When passing true as 3rd argument, the variable 
  currentAmplitude will ramp from its current value to the new target value as determined from the 
  key and velocity that is being glided to. When passing false, it will be set immediately to the 
  target value. */
  virtual void prepareForAmplitudeRamp(double newKey, double newVel, bool shouldRamp);

  /** This is supposed to trigger a note on this voice - typically this will cause the oscillator 
  frequencies to be set up and retrigger the envelope generators. */
  virtual void triggerNote(int newKey, int newVelocity, int newDetune = 0);

  /** Triggers a glide to a new note on this voice - typically this will not retrigger the envelope
  generators. You don't need to implement that function when it's not needed - in that case, it 
  will call triggerNote() by default */
  virtual void glideToNote(int newKey, int newVelocity, int newDetune = 0);

  /** Purely virtual - is supposed to trigger the release phase for this voice - typically this 
  will trigger the release phases of all involved envelope-generators. The secosn parameter is the 
  (note-on) velocity of the note to be released - this is needed for setting up the amp-ramp 
  properly when the user tweaks LevelByKey during note releases - don't confuse it with note-off 
  velocity (which is currently ignored). */
  virtual void triggerRelease(int noteToBeReleased, int noteToBeReleasedVel);

  /** An overall amplitude factor for this voice, determined by the values of 'level', 
  'levelByVel', and the velocity of the currently playing note. */
  //double amplitude;

  /** The voices output level in dB. */
  double level;

  /** The key dependence of the voices output level in dB@127. */
  double levelByKey;

  /** The velocity dependence of the voices output level in dB at velocity == 127. */
  double levelByVel;

  /** Variables to handle click-free amplitude changes when switching to a new key and/or 
  velocity. */
  double targetLevel, targetAmplitude, currentLevel, currentAmplitude, levelIncPerSample, 
    ampFactorPerSample, ampRampTime;
  int remainingAmpRampSamples;

  /** The current nominal frequency as determined by the key of the current note. */
  double targetPitch, targetFrequency;

  /** The current nominal frequency as determined by the key of the current note and possibly by a
  modification due to a glide taking place. */
  double currentPitch, currentFrequency;

  /** Range for the pitch-wheel in semitones. */
  double pitchBendRange;

  /** Value of the pitch-wheel in normalized units between -1...+1. */
  double pitchBend;

  /** Factor to multiply a frequency with, due to pitchbend. */
  double pitchBendFactor;

  /** The current frequency as determined by key, glide and the pitchwheel. */
  double currentPitchWithPitchBend, currentFrequencyWithPitchBend;

  /** The time (in ms), it takes to glide from one note to another, when glide is taking place. */
  double glideTime;

  /** The sample-rate in Hz. */
  double sampleRate;

  /** A factor by which the current frequency is to be mutiplied per sample in order to reach the
  'targetFrequency' after the specified glideTime. */
  double pitchIncPerSample, freqFactorPerSample;

  /** The remaining number of samples until the glide will be finished. */
  int remainingGlideSamples;

  /** The master tuning frequency in Hz. */
  //double masterTuneA4;

  /** A flag to indicate whether glide is active or not. */
  bool glideIsActive;

  /** A list of MIDI note events - we keep track of more than one event in order to implement
  glides which can go back and forth (a typical feature of monophonic synthesizers). */
  std::list<MidiNoteEvent> noteList;
  // maybe a std::vector would be more performant for this - see this talk:
  // https://www.youtube.com/watch?v=LrVi9LHP8Bk
  // up to some 1000s of elements, vector beats list - also, we don't want memory allocations to 
  // happen, or use RAPT::rsDoubleEndedQueue which is based on std::vector

  /** A mutex-lock to avoid threading problems with accesses to the noteList. */
  MutexLock mutex;

  /** Pointer to a TuningTable object - will be initialized with NULL and ignored by default. When
  it is assigned (by an outlying PolyphonicInstrument object) it will be used for the 
  note-to-frequency conversions instead of the standard conversion. */
  TuningTable* tuningTable;

};

//-------------------------------------------------------------------------------------------------
// inlined functions :

INLINE void PolyphonicInstrumentVoice::getSampleFrameStereo(
  double* /*outL*/, double* /*outR*/, double* /*voiceAmplitude*/)
{
  // setup the 'currentFrequency' member variable - this maybe used in subclasses to implement 
  // glide:
  if(remainingGlideSamples <= 0 || glideIsActive == false)
  {
    currentPitch     = targetPitch;
    currentFrequency = targetFrequency;
  }
  else
  {
    currentPitch     += pitchIncPerSample;
    currentFrequency *= freqFactorPerSample;
    remainingGlideSamples--;
  }
  currentPitchWithPitchBend     = currentPitch + pitchBend * pitchBendRange;
  currentFrequencyWithPitchBend = currentFrequency * pitchBendFactor;

  // setup the 'currentAmplitude' member variable - this maybe used in subclasses to implement a 
  // smooth amplitdue ramping for the case that the amplitude changes in between notes played by 
  // this voice (due to key/vel dependence of the amplitude):
  if(remainingAmpRampSamples <= 0)
  {
    currentLevel     = targetLevel;
    currentAmplitude = targetAmplitude;
  }
  else
  {
    currentLevel     += levelIncPerSample;
    currentAmplitude *= ampFactorPerSample;
    remainingAmpRampSamples--;
  }

  // increment the note age:
  currentNoteAge++;
}

} // end namespace rosic

#endif //  rosic_PolyphonicInstrumentVoice_h
