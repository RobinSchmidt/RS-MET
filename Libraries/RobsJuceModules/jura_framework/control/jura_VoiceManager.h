#ifndef jura_VoiceManager_h
#define jura_VoiceManager_h


//=================================================================================================

/** A class to dispatch various kinds of MIDI messages to specific handler functions that can be 
overriden in a subclass. */

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

  virtual int handleMidiMessage(MidiMessage msg, int voiceInfo)
  {
    int returnInfo = -1;  // code for unknown
    if( msg.isNoteOn() )
      returnInfo = noteOn(msg.getNoteNumber(), msg.getVelocity(), voiceInfo);
    else if( msg.isNoteOff() )
      returnInfo = noteOff(msg.getNoteNumber(), voiceInfo);
    else
      handleMidiMessage(msg);
    return returnInfo;
  }

  virtual int noteOn(int key, int vel, int voiceInfo) 
  { 
    noteOn(key, vel);
    return -1;        // code for unknown
  }

  virtual int noteOff(int key, int voiceInfo) 
  { 
    noteOff(key);
    return -1;        // code for unknown
  }




  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsMidiMessageDispatcher)
};
// -maybe move into extra file
// -when we make this a subclass of rsMidiMessageHandler, we get crashes - why?

//=================================================================================================

/** A class for managing polyphony in instrument modules. */

class JUCE_API rsVoiceManager : public rsMidiMessageDispatcher
//class JUCE_API rsVoiceManager : public rsMidiMessageHandler
{

public:

  rsVoiceManager()
  {
    setMaxNumVoices(16);
  }


  virtual ~rsVoiceManager() {}

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Modes for voice stealing. */
  enum class StealMode
  {
    oldest,
    newest
    //nearest, lowest, highest, quietest, dontSteal
  };

  enum class KillMode
  {
    immediately,
    afterSilence
  };

  void setKillMode(KillMode newMode) { killMode = newMode; }

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

  void setSampleRate(double newRate) { sampleRate = newRate; }


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
  int getNewestVoice() const { return newestVoice; }

  size_t getNumReleasingVoices() const { return releasingVoices.size(); }
  // the inconsistency is a bit ugly but we want to avoid conversion

  int getActiveVoiceIndex(int i) const { return activeVoices[i]; }



  double getPitchForKey(int key) const { return double(key); }
  // preliminary - todo: take into account tuning table and master-tune

  //bool needsPreRenderUpdate() const { return _needsPreRenderUpdate; } 

  //bool needsPostRenderUpdate() const { return _needsPostRenderUpdate; } 

  bool needsVoiceKillCheck() const { return getNumReleasingVoices() > 0; }


  double getVoicePitch(int i) const { return voiceStates[i].pitch; }

  double getVoiceNormalizedVelocity(int i) const { return voiceStates[i].vel01; }

  double getPitchBend() const { return pitchBend; }


  //-----------------------------------------------------------------------------------------------
  // \name Event handling


  //virtual void handleMidiMessage(const juce::MidiMessage& message, 
  //  rsMidiMessageHandler::MidiHandleInfo* info) override;


  virtual void noteOn(int key, int vel) override;

  virtual void noteOff(int key) override;

  virtual void setPitchBend(int pitchBendValue) override;




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

  void triggerVoice(int voiceIndex, int key, int vel);

  void stealVoice(int key, int vel);

  void releaseVoice(int voiceIndex);

  void deactivateVoice(int voiceindex);
  // maybe rename to killVoice


  int maxNumVoices    = 16;  // maximum number of voices
  int numVoices       =  8;  // number of available voices
  int numActiveVoices =  0;  // number of currently playing voices
  int newestVoice     =  0;  // most recently triggered voice

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
    }

    double pitch;    // midi-pitch computed from key and tuning 
    double vel01;    // velocity as continuous value, normalized to range 0..1
    int    key;      // key as MIDI note value (maybe use char)
    //bool   isHeld;   // flag to indicate that the note is being held - may not be needed
  };

  std::vector<VoiceState> voiceStates;



  StealMode stealMode = StealMode::oldest;
  KillMode  killMode  = KillMode::immediately;


  //bool _needsPreRenderUpdate  = false;
  //bool _needsPostRenderUpdate = false;
  // Flags that is set to true whenever there is some process going on that requires a per-sample
  // update of the state, such as gliding from one note to another. The underscore is just for 
  // avoiding confusion with the method that has the same name.


  double pitchBend = 0;

  // Stuff for the voice killing functionality:
  //static const int numChannels = 2; // maybe uncomment later
  double *voicesBuffer   = nullptr;     // length should be numChannels*maxNumVoices
  double sampleRate      = 44100.0;
  double killThreshold   = 0.00001;     // -100 dB by default
  double killTimeSeconds = 0.1;         //  100 ms
  int    killTimeSamples = (int) ceil(sampleRate * killTimeSeconds);
  std::vector<int> killCounters;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVoiceManager)
};

// ToDo:
// -this class should also provide a mechanism to kill voices: after a noteOff has been received 
//  for a given voice, monitor the output of that voice (how?) and if it remains below a certain
//  threshold for a given amount of time, the voice can be killed and moved back into the pool
//  of available voices
// -maybe use short or even char instead of int - that may save a little memory

// see romos::VoiceAllocator and rosic::PolyphonicInstrument/Voice for ideas, code and formulas


#endif