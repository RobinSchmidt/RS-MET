#ifndef jura_PitchShifter_h
#define jura_PitchShifter_h

class PitchShifterAudioModule : public AudioModule
{

  friend class PitchShifterModuleEditor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  PitchShifterAudioModule(CriticalSection *newPlugInLock, 
    rosic::PitchShifterGrainAdaptive *pitchShifterToWrap = nullptr);

  virtual ~PitchShifterAudioModule();

  AudioModuleEditor* createEditor() override;

  //-----------------------------------------------------------------------------------------------
  // automation and state management:

  /** Creates the static parameters for this module (i.e. parameters that are not created dynamically 
  and are thus always there). */
  virtual void createStaticParameters();

  /** Restores the state of this module from an XmlElement (which was presumably previously created 
  via getStateAsXml). */
  //virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean);

  /** Converts a state which might possibly be from an older version to the current patch-format. */
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState);

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedPitchShifter->setSampleRate(newSampleRate);
  }

  virtual void setBeatsPerMinute(double newBpm)
  {
    wrappedPitchShifter->setBeatsPerMinute(newBpm);
  }

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedPitchShifter->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    double dL, dR;
    for(int n=0; n<numSamples; n++)
    {
      dL = (double)left[n];
      dR = (double)right[n];
      getSampleFrameStereo(&dL, &dR);
      left[n]  = (float)dL;
      right[n] = (float)dR;
    }
  }

  //-----------------------------------------------------------------------------------------------
  // event processing:

  virtual void reset()
  {
    wrappedPitchShifter->reset();
  }

protected:

  rosic::PitchShifterGrainAdaptive *wrappedPitchShifter;
  bool wrappedPitchShifterIsOwned = false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class PitchShifterModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  PitchShifterModuleEditor(CriticalSection *newPlugInLock, 
    PitchShifterAudioModule* newPitchShifterAudioModule);

  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  PitchShifterAudioModule *pitchShifterModuleToEdit;

  RSlider *coarseSlider, *fineSlider, *grainLengthInMillisecondsSlider, *grainLengthInCyclesSlider,
    *grainLengthInBeatsSlider, *feedbackSlider, *dryWetSlider;
  RComboBox *grainLengthUnitComboBox;
  RButton *invertButton, *reverseButton, *antiAliasButton; // *formantPreserveButton, *monoButton;

  juce_UseDebuggingNewOperator;
};

#endif 
