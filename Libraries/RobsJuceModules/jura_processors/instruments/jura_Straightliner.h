#ifndef jura_Straightliner_h
#define jura_Straightliner_h


/** This class wraps rosic::Straightliner into a rosof::AudioModule to facilitate its use as 
plugIn. */

class StraightlinerAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class StraightlinerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  ///** Constructor. */
  //StraightlinerAudioModule(CriticalSection *newPlugInLock, 
  //  rosic::Straightliner *straightlinerToWrap); 
  // not needed anymore - for full blown isntruments, we don't need facilities to wrap around an
  // existing rosic object. that is relavant onyl for modules that are supposed to be submodules
  // like a filter in a synth
   
  /** Constructor. */
  StraightlinerAudioModule(CriticalSection *newPlugInLock);


  virtual ~StraightlinerAudioModule();

  virtual AudioModuleEditor* createEditor(int type) override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:
  /*
  virtual juce::String getDefaultPresetLocation() override 
  { 
    return getApplicationDirectory() 
      + juce::String("/StraightlinerPresets/000-InitPatchSawtooth.xml");
  }
  */

  virtual void setSampleRate(double newSampleRate) override
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->setSampleRate(newSampleRate);
  }

  /** Loads the user's custom preferences such as the sample-content path. */
  //virtual void loadPreferences();

  //---------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    if(wrappedStraightliner->isSilent())
    {
      RAPT::rsArrayTools::fillWithZeros(inOutBuffer[0], numSamples);
      RAPT::rsArrayTools::fillWithZeros(inOutBuffer[1], numSamples);
    }
    else
    {
      for(int n=0; n<numSamples; n++)
        wrappedStraightliner->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
    }
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    wrappedStraightliner->getSampleFrameStereo(left, right);
  }

  //---------------------------------------------------------------------------------------------
  // others:

  virtual void reset() override
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->resetAllVoices();
  }

  /** We must override this here because we have - for historical reasons - not the FourOscScion
  as child-module but each osc separately, so we must manually tak care of updating the preset
  field of the FourOscSectionEditor here. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the synth's building blocks 
  // here in order to make them automatable:
  FourOscSectionAudioModule      *oscSectionModule;
  MultiModeFilterAudioModule     *filterModule;
  BreakpointModulatorAudioModule *pitchEnvModule, *filterEnvModule, *ampEnvModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::Straightliner *wrappedStraightliner;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class StraightlinerModuleEditor : public PolyphonicInstrumentEditor
{

public:

  enum defaultColourSchemes
  {
    INHERITED_DEFAULTS = 0,
    RED_BLUE_PURPLE_DARK_BRIGHT_DARK,
    RED_BLUE_PURPLE_DARK_DARK_DARK,
    RED_BLUE_PURPLE_BRIGHT_BRIGHT_BRIGHT
  };

  /** Constructor. */
  StraightlinerModuleEditor(CriticalSection *newPlugInLock, 
    StraightlinerAudioModule* newStraightlinerAudioModule);

  // setup:

  //virtual void loadPreferencesFromFile();
  //virtual void setColourSchemeFromXml(const XmlElement& xmlColorScheme);

  /** Sets up the color-schemes for the embedded sub-editors according to the colour-scheme of 
  this editor. */
  virtual void updateSubEditorColourSchemes();

  // callbacks:
  virtual void copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom) override;
  virtual void updateWidgetsAccordingToState() override;
  //virtual void paint(Graphics &g) override; 
  virtual void paintOverChildren(Graphics& g) override;
  virtual void resized() override;

protected:

  /** This is the actual plugin engine which does all the dsp and automation handling. */
  StraightlinerAudioModule *straightlinerAudioModule;

  // the sub-editors:
  FourOscSectionModuleEditor     *oscSectionEditor;
  MultiModeFilterModuleEditor    *filterEditor;
  BreakpointModulatorEditorMulti *envelopeEditor;

  juce_UseDebuggingNewOperator;
};

#endif
