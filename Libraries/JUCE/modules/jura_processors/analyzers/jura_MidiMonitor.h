#ifndef jura_MidiMonitor_h
#define jura_MidiMonitor_h

class MidiMonitorAudioModule : public AudioModule, public ChangeBroadcaster
{

  friend class MidiMonitorModuleEditor;

public:

  MidiMonitorAudioModule(CriticalSection *newPlugInLock);

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void handleMidiMessage(MidiMessage message);

  virtual void getSampleFrameStereo(double *inOutL, double *inOutR) { }

protected:

  void initializeAutomatableParameters();

  juce::String      midiMessageString;
  CriticalSection   messageStringLock;
  MidiMessageFilter midiFilter;

  juce_UseDebuggingNewOperator;
};

//================================================================================================

#endif 
