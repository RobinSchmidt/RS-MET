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

  void selectBand(int bandToSelect) { selectedBand = bandToSelect; }

  int getSelectedBand() { return selectedBand; }

  int getMaxNumBands() { return maxNumBands; }

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;

  // target functions for the per-band parameter callbacks:
  void setSplitFreq(double newFreq)      { multiCompCore.setSplitFrequency(selectedBand, newFreq);      }
  void setThreshold(double newThreshold) { multiCompCore.setThreshold(     selectedBand, newThreshold); }
  void setRatio(    double newRatio)     { multiCompCore.setRatio(         selectedBand, newRatio);     }
  void setAttack(   double newAttack)    { multiCompCore.setAttackTime(    selectedBand, newAttack);    }
  void setRelease(  double newRelease)   { multiCompCore.setReleaseTime(   selectedBand, newRelease);   }

protected:

  rosic::rsMultiBandCompressor multiCompCore;

  int maxNumBands  =  0;  // assigned in constructor
  int selectedBand =  0;  // -1 is code for "None"


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompAudioModule)
};

//=================================================================================================

class MultiCompModuleEditor : public AudioModuleEditor, public RComboBoxObserver
{

public:

  MultiCompModuleEditor(MultiCompAudioModule* multiCompModuleToEdit);

  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;
  virtual void resized() override;

protected:


  virtual void createWidgets();

  /** Makes currently required widgets visible and currently not required widgets invisible. */
  virtual void updateWidgetVisibility();

  MultiCompAudioModule *multiCompModule;

  // widgets:
  RSlider *numBandsSlider;
  std::vector<RSlider*> splitFreqSliders, thresholdSliders, ratioSliders, attackSliders, 
    releaseSliders;
  //std::vector<RButton*> editButtons, soloButtons, muteButtons;
  RComboBox *splitModeBox, *bandSelectBox;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompModuleEditor)
};

#endif