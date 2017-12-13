#ifndef jura_LowFrequencyOscillator_h
#define jura_LowFrequencyOscillator_h


class LowFrequencyOscillatorAudioModule : public AudioModule
{

  friend class LowFrequencyOscillatorEditor;
  friend class LowFrequencyOscillatorEditorCompact;

public:

  LowFrequencyOscillatorAudioModule(CriticalSection *newPlugInLock, 
    rosic::LowFrequencyOscillator *newLowFrequencyOscillatorToWrap);


  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;  // to be deprecated

  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  virtual void setWaveformFromFile(const juce::String &fileToLoadFrom);

  virtual void setBeatsPerMinute(double newBpm)
  {
    wrappedLowFrequencyOscillator->setBeatsPerMinute(newBpm);
  }

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override 
  {
    for(int n = 0; n < numSamples; n++)
      for(int i = 0; i < numChannels; i++)
        inOutBuffer[i][n] = wrappedLowFrequencyOscillator->getSample();
  }

  //virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  //{
  //  *inOutL = *inOutR = wrappedLowFrequencyOscillator->getSample();
  //}

protected:

  /** Fills the array of automatable parameters. */
  virtual void initializeAutomatableParameters();

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::LowFrequencyOscillator *wrappedLowFrequencyOscillator;

  // child modules:
  WaveTableAudioModule *waveTableModule;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** This is a class for a context menu that can be opened to show and edit a comprehensive set of
parameters for a stereo-oscillator. It is used by LowFrequencyOscillatorEditor by clicking on the 
'More' button. */

class LowFrequencyOscillatorEditor : public AudioModuleEditor, public ChangeBroadcaster
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  LowFrequencyOscillatorEditor(CriticalSection *newPlugInLock, 
    LowFrequencyOscillatorAudioModule* newOscillatorModuleToEdit); 

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the bounds of the popup editor relative to the top-left position of the 
  edit-button. */
  virtual void setPopUpEditorBounds(int x, int y, int w, int h)
  { waveTableEditor->setPopUpEditorBounds(x, y, w, h); }

  //---------------------------------------------------------------------------------------------
  // callbacks:
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  // pointers to the edited objects (wrapped and non-wrapped):
  rosic::LowFrequencyOscillator     *oscillatorToEdit;
  LowFrequencyOscillatorAudioModule *oscillatorModuleToEdit;

  WaveTableModuleEditorCompact      *waveTableEditor;

  RSlider *cycleLengthSlider, *depthSlider;
  RButton *tempoSyncButton; 

  juce_UseDebuggingNewOperator;
};



#endif 
