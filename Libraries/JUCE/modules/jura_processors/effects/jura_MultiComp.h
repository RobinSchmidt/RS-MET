#ifndef jura_MultiComp_h
#define jura_MultiComp_h


class JUCE_API MultiCompAudioModule : public jura::ModulatableAudioModule
{

public:

  MultiCompAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);

  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;

protected:

  rosic::rsMultiBandCompressor multiCompCore;

  int maxNumBands  = 0;  // assigned in constructor
  int selectedBand = 0;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompAudioModule)
};

//=================================================================================================

class MultiCompModuleEditor : public AudioModuleEditor
{

public:

  MultiCompModuleEditor(MultiCompAudioModule* multiCompModuleToEdit);

  virtual void resized();

protected:


  virtual void createWidgets();

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  MultiCompAudioModule *multiCompModule;

  // widgets:
  std::vector<RSlider*> splitFreqSliders, thresholdSliders, ratioSliders, attackSliders, 
    releaseSliders;
  std::vector<RButton*> editButtons, soloButtons, muteButtons;
  RComboBox *boxSplitMode;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompModuleEditor)
};

#endif