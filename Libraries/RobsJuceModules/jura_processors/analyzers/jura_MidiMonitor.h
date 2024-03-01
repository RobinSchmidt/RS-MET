#ifndef jura_MidiMonitor_h
#define jura_MidiMonitor_h

class MidiMonitorAudioModule : public AudioModuleWithMidiIn, public ChangeBroadcaster
{

  friend class MidiMonitorModuleEditor;

public:

  MidiMonitorAudioModule(CriticalSection *newPlugInLock);

  AudioModuleEditor* createEditor(int type) override;

  void parameterChanged(Parameter* parameterThatHasChanged) override;

  void handleMidiMessage(MidiMessage message) override;

  void handleMidiMessageForVoice(MidiMessage msg, int voice) override
  { handleMidiMessage(msg); }
  // This is getting called instead of handleMidiMessage in ToolChain when the msg is a note event
  // because note events need a voice dispatch. To not miss note events, we need to override 
  // handleMidiMessageForVoice as well and pass the msg through to handleMidiMessage


  virtual void getSampleFrameStereo(double *inOutL, double *inOutR) { }

protected:

  void initializeAutomatableParameters();

  juce::String      midiMessageString;
  CriticalSection   messageStringLock;
  MidiMessageFilter midiFilter;

  juce_UseDebuggingNewOperator;
};

//================================================================================================

class MidiMonitorModuleEditor : public AudioModuleEditor
{

public:

  MidiMonitorModuleEditor(CriticalSection *newPlugInLock, 
    MidiMonitorAudioModule* newMidiMonitorAudioModule);
  virtual ~MidiMonitorModuleEditor();

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();

protected:

  /** Clears the screen and optionally also clears the messageString in the underlying 
  AudioNodule. */
  virtual void clearScreen(bool clearAlsoMessageString = true);

  /** Updates the content of the screen according to the messageString in the underlying 
  AudioModule. */
  virtual void updateScreen();

  MidiMonitorAudioModule *midiMonitorModuleToEdit;

  RTextField *eventFilterLabel;

  RButton *noteButton, *controllerButton, *pitchWheelButton, *programChangeButton, 
    *aftertouchButton, *channelPressureButton, *sysExButton, *metaEventButton, 
    *transportButton, *songPositionButton, *machineControlButton, *activeSenseButton,
    *clockButton, *otherButton, *clearButton;

  RTextEditor *outputDisplay;

  juce_UseDebuggingNewOperator;
};

#endif 
