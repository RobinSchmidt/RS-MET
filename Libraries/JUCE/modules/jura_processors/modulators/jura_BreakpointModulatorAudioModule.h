#ifndef jura_BreakpointModulatorAudioModule_h
#define jura_BreakpointModulatorAudioModule_h

/** This class wraps a RAPT::BreakpointModulator into a jura::AudioModule to facilitate its use as 
plugIn or sub-module inside a plugIn. */

class JUCE_API BreakpointModulatorAudioModule 
  : /*public AudioModule*/ public AudioModuleWithMidiIn , public ModulationSource
{

  friend class BreakpointModulatorEditor;
  friend class BreakpointModulatorGlobalEditor;
  friend class BreakpointParameterEditor;
  friend class BreakpointModulatorEditorMulti;
  friend class BreakpointModulatorEditorCompact;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  BreakpointModulatorAudioModule(CriticalSection* newPlugInLock, 
    rosic::BreakpointModulator* newBreakpointModulatorToWrap = nullptr);

  /** Destructor. */
  virtual ~BreakpointModulatorAudioModule();

  //---------------------------------------------------------------------------------------------
  // overrides:

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  virtual void setStateToDefaults() override;

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    *inOutR = *inOutL = wrappedBreakpointModulator->getSample();
  }
  // maybe, this function is obsolete now

  virtual void updateModulationValue() override
  {
    modValue = wrappedBreakpointModulator->getSample();
  }

  virtual AudioModuleEditor *createEditor() override;


  // new overrides (added after dragging the old code over - they are currently only dummies):
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override 
  {
    // we should do nothing here, the computation is done in updateModulationValue()
    // actually, the override may be deleted then
  }

  virtual void noteOn(int noteNumber, int velocity) override 
  {
    wrappedBreakpointModulator->noteOn(true, noteNumber, velocity);
  }

  virtual void noteOff(int noteNumber) override
  {
    wrappedBreakpointModulator->noteOff(true);
  }

  virtual void setBeatsPerMinute(double bpm) override
  {
    wrappedBreakpointModulator->setBeatsPerMinute(bpm);
  }

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying RAPT object which is wrapped. */
  //RAPT::rsBreakpointModulator<double> *wrappedBreakpointModulator;
  rosic::BreakpointModulator *wrappedBreakpointModulator;

  bool wrappedBreakpointModulatorIsOwned = false;


  juce_UseDebuggingNewOperator;
};


#endif 
