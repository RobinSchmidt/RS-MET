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

  virtual AudioModuleEditor* createEditor() override;

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->setSampleRate(newSampleRate);
  }

  /** Loads the user's custom preferences such as the sample-content path. */
  //virtual void loadPreferences();

  //---------------------------------------------------------------------------------------------
  // audio processing:

  // i think, this code is obsolete now:
  ///** Calculates a stereo-ouput frame. */
  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  wrappedStraightliner->getSampleFrameStereo(inOutL, inOutR);
  //}

  //virtual void processBlockStereo(float *left, float *right, int numSamples)
  //{
  //  if(wrappedStraightliner->isSilent())
  //  {
  //    rosic::fillWithZeros(left, numSamples);
  //    rosic::fillWithZeros(right, numSamples);
  //  }
  //  else
  //  {
  //    for(int n=0; n<numSamples; n++)
  //    {
  //      double dL = (double)left[n];
  //      double dR = (double)right[n];
  //      wrappedStraightliner->getSampleFrameStereo(&dL, &dR);
  //      left[n]  = (float)dL;
  //      right[n] = (float)dR;
  //    }
  //  }
  //}

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    if(wrappedStraightliner->isSilent())
    {
      rosic::fillWithZeros(inOutBuffer[0], numSamples);
      rosic::fillWithZeros(inOutBuffer[1], numSamples);
    }
    else
    {
      for(int n=0; n<numSamples; n++)
        wrappedStraightliner->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
    }
  }

  //---------------------------------------------------------------------------------------------
  // others:

  virtual void reset()
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->resetAllVoices();
  }

  /** We must override this here because we have - for historical reasons - not the FourOscScion
  as child-module but each osc separately, so we must manually tak care of updating the preset
  field of the FourOscSectionEditor here. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

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
  virtual void copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom);
  virtual void updateWidgetsAccordingToState();
  virtual void paint(Graphics &g); 
  virtual void resized();

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
