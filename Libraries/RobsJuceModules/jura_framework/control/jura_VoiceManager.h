#ifndef jura_VoiceManager_h
#define jura_VoiceManager_h



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


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsMidiMessageDispatcher)
};
// maybe move into extra file
// AudioModuleWithMidiIn should then also become a subclass of this



//=================================================================================================

/** A class for managing polyphony in instrument modules. */

class JUCE_API rsVoiceManager : public rsMidiMessageDispatcher
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
    // to be released naturally
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getMaxNumVoices() const { return maxNumVoices; }

  int getNumVoices() const { return numVoices; }

  int getNumActiveVoices() const { return numActiveVoices; } 

  int getNumIdleVoices() const { return numVoices - numActiveVoices; }


  //-----------------------------------------------------------------------------------------------
  // \name Event handling


  virtual void noteOn(int key, int vel) override;

  virtual void noteOff(int key) override;

  virtual void setPitchBend(int pitchBendValue) override;




  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Should be called once per sample by some outlying driver class. */
  void perSampleUpdate();


  void reset();


protected:

  void stealVoice(int key, int vel);


  int maxNumVoices    = 16;  // maximum number of voices
  int numVoices       =  8;  // number of available voices
  int numActiveVoices =  0;  // number of currently playing voices

  /** Indices of the voices that are currently playing and therefore must process audio. A voice 
  is playing if it's currently being held or releasing. */
  std::vector<int> activeVoices;

  /** Indices of the voices that are currently idle and therefore available for new notes. */
  std::vector<int> idleVoices;

  // The invariant should be that playingVoices and idleVoices at all times combine to the full
  // set 0,1,2,3,...,numVoices-1. The playingVoices is supposed to be filled up to 
  // n = numActiveVoices-1 and the idleVoice is filled up to 
  // numVoices - numActiveVoices = numVoices-1-n = 


  //std::vector<int> releasingVoices;


  struct VoiceState
  {
    VoiceState() { reset(); }

    void reset()
    {
      frequency = 0;
      normalizedVelocity = 0;
      key = 0;
      isHeld = false;
    }

    double frequency;           // derived from key and the tuning to be used - maybe pitch is better
    double normalizedVelocity;  // velocity normalized to the range 0...1
    int    key;                 // key as MIDI note value
    bool   isHeld;              // flag to indicate that the note is being held - may not be needed
  };

  std::vector<VoiceState> voiceStates;



  StealMode stealMode = StealMode::oldest;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVoiceManager)
};

// ToDo:
// -this class should also provide a mechanism to kill voices: after a noteOff has been received 
//  for a given voice, monitor the output of that voice (how?) and if it remains below a certain
//  threshold for a given amount of time, the voice can be killed and moved back into the pool
//  of available voices

// see romos::VoiceAllocator and rosic::PolyphonicInstrument/Voice for ideas, code and formulas


#endif