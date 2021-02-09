#ifndef jura_VoiceManager_h
#define jura_VoiceManager_h


//=================================================================================================

/** A class to dispatch various kinds of MIDI messages to specific handler functions that can be 
overriden in a subclass.  */

class JUCE_API rsMidiMessageDispatcher
{

public:

  rsMidiMessageDispatcher() {}

  virtual ~rsMidiMessageDispatcher() {}





  //-----------------------------------------------------------------------------------------------
  // \name Event processing:

  /** Handles a generic MidiMessage. This dispatches to a call of the appropriate specific handler 
  function. */
  virtual void handleMidiMessage(MidiMessage message);
  // Last time i checked sizeof(MidiMessage) was 24 bytes = 192 bits so i think it makes sense
  // to pass by reference. A pointer/reference is just 8 bytes = 64 bits.
  // https://stackoverflow.com/questions/40185665/performance-cost-of-passing-by-value-vs-by-reference-or-by-pointer
  // https://softwareengineering.stackexchange.com/questions/372105/is-passing-arguments-as-const-references-premature-optimization
  // https://www.cplusplus.com/articles/z6vU7k9E/

  // should return some information, how the message was handled - in particular, which voice was 
  // used for note-on events

  /** Triggered by a note-on event. */
  virtual void noteOn(int noteNumber, int velocity) {}
  // return an int for the used voice

  /** Triggered by a note-off event. */
  virtual void noteOff(int noteNumber) {}
  // todo: support note-off velocity

  /** Triggered by an all-notes-off event. */
  virtual void allNotesOff() {}

  /** Overrides setMidiController which is inherited from both base-classes - and we simply we pass
  through the function call to both of them here. */
  virtual void setMidiController(int controllerNumber, float controllerValue) {}

  /** Triggered by a pitch-bend event. */
  virtual void setPitchBend(int pitchBendValue) {}

  /** Triggered by an aftertouch event. */
  virtual void setAfterTouch(int afterTouchValue) {}

  /** Triggered by a channel pressure event. */
  virtual void setChannelPressure(int channelPressureValue) {}





  // new code for supporting events with voice-info:

  /** This function can be overriden alternatively to the regular handleMidiMessage callback in 
  cases when the callee needs to know the voice to which this event applies as additional 
  information. The default implementation will just call the regular handleMidiMessage callback in 
  case of non-note messages and dispatch to noteOn/OffForVoice in case of note messages. These
  two functions in turn have default implementations that also just call their regular, voiceless
  counterparts noteOn/Off but subclasses can override them too. For example, look at the polyphonic
  modulator classes such as AttackDecayEnvelopeModulePoly. It overrides noteOnForVoice and triggers 
  the envelope in the underlying DSP object for the given voice. It can also be used to reset the 
  phase in polyphonic oscillators, reset the state in polyphonic filters, etc. The general idea 
  is that many modules want to receive the same message but only the first one allocates the voice 
  and the other ones need to have the information which voice that was. In ToolChain, the message 
  is first passed to the voiceManager via handleMidiMessageReturnVoice which then allocates the 
  voice and returns the info, which voice it has selected (because rsVoiceManager overrides that
  method accordingly). Then, the same event is also passed to the child modules by calling 
  handleMidiMessageWithVoiceInfo, so they can retrigger the right voice, if necessary. A subclass
  that overrides it should also handle the case when voice = 1 which is used as code for 
  "unknown"..or should it?...maybe we sho. */
  virtual void handleMidiMessageForVoice(MidiMessage msg, int voice)
  {
    if( msg.isNoteOn() )
      noteOnForVoice(msg.getNoteNumber(), msg.getVelocity(), voice);
    else if( msg.isNoteOff() )
      noteOffForVoice(msg.getNoteNumber(), voice);
    else
      handleMidiMessage(msg);
  }
  // todo: maybe have an addtional parameter for the midi channel (defaulting to 1)..or no: that
  // info is already contained in the MidiMessage
  // hmm...maybe this is is not a good idea - just override the regular noteOn/Off and infer the
  // voice from the voiceManager - or, in AudioModulePoly override this instead of noteOn


  /** This is supposed to be used when the caller needs to know which voice was assigned or 
  released for noteOn and noteOff events. But the default implementation will just return -1 as
  code for "unknown". Callers should gracefully handle such an "unknown" return value. Subclasses 
  can override this function to provide that information. */
  virtual int handleMidiMessageReturnVoice(MidiMessage msg)
  {
    int returnInfo = -1;  // code for unknown
    if( msg.isNoteOn() )
      returnInfo = noteOnReturnVoice(msg.getNoteNumber(), msg.getVelocity());
    else if( msg.isNoteOff() )
      returnInfo = noteOffReturnVoice(msg.getNoteNumber());
    else
      handleMidiMessage(msg);
    return returnInfo;
  }

  // We use some negative numbers as codes to be used as return values when a valid voice index is
  // expected to be returned but we can't, either because we don't have the information or no voice
  // could be found/allocated, or more than one was found (in case of noteOff, for example):
  static const int unknownVoice = -1;
  static const int noVoice      = -2;
  static const int manyVoices   = -3;

  /** Can be overriden by subclasses to return the voice index that was used to play the note. For
  example, the subclass rsVoiceManager does this. The baseclass implementation returns -1 as code 
  for "unknown". */
  virtual int noteOnReturnVoice(int key, int vel) { noteOn(key, vel); return unknownVoice; }

  /** Can be overriden by subclasses to return the voice index that was put into release state, 
  i.e. the voice which was currently holding the given note before receiving the event. The 
  baseclass implementation returns -1 as code for "unknown".  
  ...to consider: what if a note-off is received and there was actually no voice playing that 
  note? should a subclass then also return -1 or whould we use another code for "none" - which 
  really is a condition different from "unknown". ...also: could it be that more than voice was 
  playing the same note? ...and the note-off could trigger multiple voices to enter release? what 
  should a subclass return in such a case? maybe another code for "many"? for the time being, we 
  just assume that such cases are not a thing. */
  virtual int noteOffReturnVoice(int key) { noteOff(key); return unknownVoice; }
  


  virtual void noteOnForVoice(int key, int vel, int voice) { noteOn(key, vel); }

  virtual void noteOffForVoice(int key, int voice) { noteOff(key); }


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsMidiMessageDispatcher)
};
// -maybe move into extra file
// -when we make this a subclass of rsMidiMessageHandler, we get crashes - why? ...update: this 
//  baseclass does not exits for that very reason - but it's strange anyway and it deserves some 
//  more research why this happened...

//=================================================================================================

/** A class for managing polyphony in instrument modules. 

This class is still under construction - not all features that the API provides are implemented 
yet... */

class JUCE_API rsVoiceManager : public rsMidiMessageDispatcher
{

public:

  rsVoiceManager()
  {
    //setMaxNumVoices(16);   // seems a good default
    setMaxNumVoices(4);      // for debug
  }


  virtual ~rsVoiceManager() {}

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Modes for voice stealing. */
  enum class StealMode
  {
    oldest,
    newest,
    noSteal
    // todo: oldestInRelease, nearest, farthest, lowest, highest, quietest, 
  };
  void setStealMode(StealMode newMode) { stealMode = newMode; }

  enum class KillMode
  {
    immediately,
    afterSilence
  };
  /** Sets up how voices are killed, i.e. determines the criteria by which it is decided, if a 
  voice can be deactivated. Voices are typically not deactivated immediately after receiving a 
  noteOff for the key the voice is playing. Instead, that voice enters the release phase and there
  is no other event that can tell us, when the release phase is over, so we need to figure it out
  ourselves. In the mode KillMode::immediately, we actually do indeed kill the voice immediately 
  but that mode is mostly for testing purposes. In real-world applications, we would rather apply 
  the following heuristic: If the absolute value of the voice's output remains below a given 
  threshold for some given amount of time, the voice is considered to have finished its release and
  can be turned off. That's what KillMode::afterSilence does. The threshold and time are set up via
  setKillThreshold and setKillTime functions.  */
  void setKillMode(KillMode newMode) { killMode = newMode; }

  void setKillThreshold(double newThreshold) { killThreshold = newThreshold; }

  void setKillTime(double newTime)
  {
    killTimeSeconds = newTime;
    killTimeSamples = (int) ceil(sampleRate * killTimeSeconds);
  }


  enum class RetriggerMode
  {
    useNewVoice,
    reuseOldVoice
  };
  /** Sets up the behavior when a note is triggered for which we have an existing active voice 
  (presumably in release). In such cases, it may make sense to not allocate a fresh voice but 
  rather to revive the existing voice and retrigger its envelopes (but not reset oscillator 
  phases and filter states). This can be done by choosing RetriggerMode::reuseOldVoice and it
  is the default behavior. If the old voice should keep releasing and a new voice should be used 
  for the new note, choose RetriggerMode::useNewVoice. */
  void setRetriggerMode(RetriggerMode newMode) { retriggerMode = newMode; }


  /** Sets the maximum number of voices that should be supported. The function is supposed to be 
  called once shortly after construction and then that setting should remain fixed for the lifetime
  of the plugin. If we later want to allow the user to change that setting at runtime, we will need 
  facilities to trigger a re-allocation of all the required resources (buffer, DSP-objects, 
  etc.). We don't have such facilities yet and maybe that feature is not worth the increased 
  complexity anyway...we'll see... */
  void setMaxNumVoices(int newNumber);

  /** Sets the number of voices that should be available for playing. */
  void setNumVoices(int newNumber)
  {
    numVoices = RAPT::rsMin(newNumber, maxNumVoices);

    //numActiveVoices = RAPT::rsMin(numActiveVoices, numVoices); 
    // hmmm...might be a bit harsh to kill off voices immediately - maybe we should wait for them 
    // to be released naturally - but this may mess with our invariant (see comment at idleVoices)
    // todo: figure out what happens when numVoices is being changed during playing and handle it
    // gracefully
  }

  /** Sets the buffer into which the outlying driver object writes its voice outputs. We need 
  access to the outputs of the individual voices here in order to be able to figure out, when
  a voice may be killed. It is currently assumed that the length is 2*maxNumVoices where the 
  factor two comes from the two channels for stereo signals - so we assume a driver that 
  produces stereo output, as ToolChain does. 
  ToDo: make that more flexible to allow any number of channels  */
  void setVoiceSignalBuffer(double* buffer) { voicesBuffer = buffer; }

  void setSampleRate(double newRate) 
  { 
    sampleRate = newRate; 
    killTimeSamples = (int) ceil(sampleRate * killTimeSeconds);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getMaxNumVoices() const { return maxNumVoices; }

  int getNumVoices() const { return numVoices; }

  int getNumActiveVoices() const { return numActiveVoices; } 

  int getNumIdleVoices() const { return numVoices - numActiveVoices; }

  /** Returns the index of the voice that was triggered most recently. This facilitates 
  compatibility with monophonic modules because they always only look at what the newest voice
  does. Initially (or after reset), when no voice has yet been triggered at all, it will return 
  zero. It will still continue to return the most recently triggered voice index even when all
  voices are already dead. It doesn't care at all about note-off. */
  int getNewestVoice() const { return newestVoice;  }
  // maybe it should return -1 when no voice is active. it doesn't seem to make sense to report
  // the index of the last voice ever used - using that information for anything will produce bad 
  // behavior such as dependency of the sound of a patch on the last note played (if note-based 
  // modulators are connected), even if the note event was hours ago

  int getNewestActiveVoice() const 
  { 
    if(numActiveVoices == 0)
      return -1;
    else
      return activeVoices[numActiveVoices-1];  
  }


  int getNumReleasingVoices() const { return numReleasingVoices; }

  //size_t getNumReleasingVoices() const { return releasingVoices.size(); }
  // The inconsistency is a bit ugly but we want to avoid conversion because conversion is not free
  // and this is called per sample. Maybe make it more consistent by using size_t for 
  // getNumActiveVoices, too. We may remove the numActiveVoices member and use 
  // activeVoices.size by adjusting the array size
  //

  int getActiveVoiceIndex(int i) const { return activeVoices[i]; }



  double getPitchForKey(int key) const { return double(key); }
  // preliminary - todo: take into account tuning table and master-tune

  //bool needsPreRenderUpdate() const { return _needsPreRenderUpdate; } 

  //bool needsPostRenderUpdate() const { return _needsPostRenderUpdate; } 

  bool needsVoiceKillCheck() const { return getNumReleasingVoices() > 0; }


  double getVoicePitch(int i) const { return voiceStates[i].pitch; }

  double getVoiceNormalizedVelocity(int i) const { return voiceStates[i].vel01; }
  // rename to getVoiceVelocity and just always assume normalization - midi values in the range
  // 0..127 are historical baggage that is not used in this modern software framework. 
  // juce::MidiMessage represents velocities as float in 0..1, too


  double getPitchBend() const { return pitchBend; }

  /** Helper function to quantize values in the range 0..1 to values representable by the midi
  standard for velocities and controller values. In midi, these are represented as 7-bit integers 
  in the range 0..127, so converting a normalized 0..1 value back and forth introduces a 
  quantization that is replicated by this function. */
  static double quantize7BitUnsigned(double val) { return round(val * 127.0) / 127.0; }
  // maybe rename to quantize7bit
  // can be sued in functions like setNoteVelocity(int voiceIndex, bool quantizeToMidi), if we want
  // to implement such a thing - but setting/readjusting note-velocities after the fact is not 
  // really supposed to be a common occurrence...but it may actually be a useful thing to 
  // "modulate" the note velocity by some sort of controller - this will typically modulate 
  // multiple parameter low-level parameters at once (gain, cutoff, env-excursions, etc.) - we'll 
  // see

  // todo: quantizeSigned8192: 14-bit values used for pitch-wheel - should expect values in -1..+1


  int getKillTimeSamples() const { return killTimeSamples; }


  /** Performs a couple of sanity checks, verifying some invariants that should always be true and
  returns true, iff everything seems same. This is mainly for debugging purposes. */
  bool isInSaneState();

  //-----------------------------------------------------------------------------------------------
  // \name Event handling


  //virtual void handleMidiMessage(const juce::MidiMessage& message, 
  //  rsMidiMessageHandler::MidiHandleInfo* info) override;

  void noteOn(int key, int vel) override 
  { 
    noteOnReturnVoice(key, vel); 
  }

  void noteOff(int key) override 
  { 
    noteOffReturnVoice(key); 
  }


  virtual int noteOnReturnVoice(int key, int vel) override;

  virtual int noteOffReturnVoice(int key) override;

  virtual void setPitchBend(int pitchBendValue) override;


  //virtual int noteOn(int key, int vel, int voice) override;

  //virtual int noteOff(int key, int voice) override;




  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Should be called once per sample by some outlying driver class before the audio signal for 
  the given sample instant is rendered. */
  //void perSampleUpdatePreRender();

  /** Should be called once per sample by some outlying driver class after the audio signal for 
  the given sample instant was rendered. */
  //void perSampleUpdatePostRender();

  void findAndKillFinishedVoices();


  void reset();


protected:

  int activateAndGetLastIdleVoice();

  int stealVoice(int key, int vel);

  void triggerVoice(int voiceIndex, int key, int vel);



  void releaseVoice(int voiceIndex);

  void deactivateVoice(int voiceindex);
  // maybe rename to killVoice

  /** Returns the instantaneous output amplitude of voice i, defined as the absolute maximum of all
  channels (currently we assume 2 channels). */
  inline double getVoiceAmplitude(int i) const
  {
    double *vb = voicesBuffer;;
    return RAPT::rsMax(RAPT::rsAbs(vb[2*i]), RAPT::rsAbs(vb[2*i+1]));
  }




  int maxNumVoices       = 16;  // maximum number of voices
  int numVoices          =  8;  // number of available voices
  int numActiveVoices    =  0;  // number of currently playing voices (holding or releasing)
  int numReleasingVoices =  0;  // number of voices in release phase
  int newestVoice        =  0;  // most recently triggered voice

  /** Indices of the voices that are currently active and therefore must process audio. A voice 
  is active if it's either currently being held or releasing. */
  std::vector<int> activeVoices;

  /** Indices of the voices that are currently idle and therefore available for new notes. The 
  invariant should be that activeVoices and idleVoices at all times combine to the full set 
  0,1,2,3,...,numVoices-1. The activeVoices array is supposed to be filled up to 
  n = numActiveVoices-1 (the rest being set to -1) and the idleVoices array is filled up to 
  numVoices - numActiveVoices = numVoices-1-n  ..verify this... */
  std::vector<int> idleVoices;

  /** Indices of voices that are currently in release mode. It's a subset of the activeVoices. 
  These are ones whose output signals need to be monitored in order to know when a voice can be 
  deactivated. */
  std::vector<int> releasingVoices;


  struct VoiceState
  {
    VoiceState() { reset(); }

    void reset()
    {
      pitch  = 0.0;
      vel01  = 0.0;
      key    = 0;
      //isHeld = false;
      // maybe we should have a state that can be: active/inactive/revived or alive/dead/revived
      // to facilitate the different retrigger modes
    }

    double pitch;    // midi-pitch computed from key and tuning 
    double vel01;    // velocity as continuous value, normalized to range 0..1
    int    key;      // key as MIDI note value (maybe use char)
    //bool   isHeld;   // flag to indicate that the note is being held - may not be needed
  };

  std::vector<VoiceState> voiceStates;


  StealMode     stealMode     = StealMode::oldest;
  KillMode      killMode      = KillMode::afterSilence;
  //KillMode      killMode      = KillMode::immediately;
  //RetriggerMode retriggerMode = RetriggerMode::reuseOldVoice;
  RetriggerMode retriggerMode = RetriggerMode::useNewVoice;


  // stuff to do
  //bool _needsPreRenderUpdate  = false;
  //bool _needsPostRenderUpdate = false;
  // Flags that is set to true whenever there is some process going on that requires a per-sample
  // update of the state, such as gliding from one note to another. The underscore is just for 
  // avoiding confusion with the method that has the same name.
  // maybe for efficiency have just one flag needsPerSampleUpdate and do all updates either
  // pre-render or post render, if possible - it will make a difference of a one sample delay only
  // and things can probably be arranged such that they do the appropriate thing in either of these
  // variants. perhaps pre-render is more appropriate, so a note-off received at one sample will
  // kill the note off *at* that sample (if killmode is immediate) and not one sample later

  //bool sustainPedalHeld   = false;
  //bool softPedalHeld      = false;  // what whould this do?
  //bool sostenutoPedalHeld = false;  // and this?




  double pitchBend = 0;

  // Stuff for the voice killing functionality:
  //static const int numChannels = 2; // maybe uncomment later
  double *voicesBuffer   = nullptr;     // length should be numChannels*maxNumVoices
  double sampleRate      = 44100.0;
  double killThreshold   = 0.0001;       // -60 dB by default
  double killTimeSeconds = 0.5;         //  100 ms
  int    killTimeSamples = (int) ceil(sampleRate * killTimeSeconds);
  std::vector<int> killCounters;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVoiceManager)
};




#endif