#ifndef jura_MultiComp_h
#define jura_MultiComp_h


/** Baseclass for multiband effects. */

class JUCE_API MultiBandEffect : public jura::ModulatableAudioModule, public ChangeBroadcaster
{

public:

  MultiBandEffect(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  /** Subclasses should call this *once* in their constructor with a pointer to the concrete 
  multiband effect (subclass). This will also create the splitting-related parameters */
  void setEffectCore( rosic::rsMultiBandEffect* effectCore);

  /** Creates the parameters related to the band-splitting. Called from setEffectCore. */
  virtual void createSplittingParameters();

  virtual void parameterChanged(Parameter* p) override;

  void selectBand(int bandToSelect) { selectedBand = bandToSelect; }


  /** Sets the splitting frequency for the band with given index to the given frequency. */
  void setSplitFreq(int bandIndex, double newFreq);

  /** Sets the splitting frequency for the selected band to the given frequency. */
  void setSplitFreq(double newFreq);




  int getSelectedBand() const { return selectedBand; }

  int getMaxNumBands() const { return maxNumBands; }

  int getBandContainingFrequency(double freq);



  /** Returns a pointer to our core DSP object. */
  rosic::rsMultiBandEffect* getCore() { return core; }

  /** Returns true, if the splitting frequencies between the bands are in (strictly) increasing 
  order. */
  bool areBandsInIncreasingOrder(bool strictly = false);

protected:

  rosic::rsMultiBandEffect* core = nullptr;

  int maxNumBands  =  0;  // assigned in constructor
  int selectedBand =  0;  // -1 is code for "None"

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandEffect)
};

//=================================================================================================

/** Plot editor for multiband effects. Allows to select a band by clicking into the rectangular 
area in the frequency response plot that represents that band, change split frequencies by dragging
vertical line, etc. */

class JUCE_API MultiBandPlotEditor : public ColourSchemeComponent, public ChangeListener
{

public:

  MultiBandPlotEditor(jura::MultiBandEffect* multiBandEffectModuleToEdit);
  virtual ~MultiBandPlotEditor();

  //virtual void parameterChanged(Parameter* p) override;
  virtual void changeListenerCallback(ChangeBroadcaster* source) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void paintOverChildren(Graphics& g) override;
  virtual void resized() override;

protected:

  rosic::rsMultiBandEffect* core;
  jura::MultiBandEffect* module; 

  rsFunctionPlot* freqRespPlot;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandPlotEditor)
};

//=================================================================================================

/** Multiband compressor with up to 16 bands. 

todo: 
 -factor out a class MultiBandEffect (in jura (done) and rosic (done) )
 -write the plot-editor for the MultiBandEffect baseclass (done)
 -derive MultiCompPlotEditor from that (done)
 -ensure that the split frequencies are always sorted from low to high
 -restrict ranges for the split-freqs according to the neighbours
  (maybe do these things in rosic::MultiBandEffect)
 -plot frequency responses: rosic::rsMultiBandEffect::getMagnitudeAt(index, freq)

*/

class JUCE_API MultiCompAudioModule : public jura::MultiBandEffect
{

public:

  MultiCompAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);

  /** Creates the per band compression parameters (threshold, ratio, attack, release). */
  virtual void createCompressionParameters();

  /** Returns a pointer to our core DSP object. */
  rosic::rsMultiBandCompressor* getMultiCompCore() { return &multiCompCore; }

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;

  // target functions for the per-band parameter callbacks:
  void setThreshold(double newThreshold) { multiCompCore.setThreshold(     selectedBand, newThreshold); }
  void setRatio(    double newRatio)     { multiCompCore.setRatio(         selectedBand, newRatio);     }
  void setAttack(   double newAttack)    { multiCompCore.setAttackTime(    selectedBand, newAttack);    }
  void setRelease(  double newRelease)   { multiCompCore.setReleaseTime(   selectedBand, newRelease);   }

protected:

  rosic::rsMultiBandCompressor multiCompCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompAudioModule)
};

//=================================================================================================

/** Plot editor for the multiband compressor. */

class JUCE_API MultiCompPlotEditor : public MultiBandPlotEditor
{

public:

  MultiCompPlotEditor(jura::MultiCompAudioModule* multiCompModuleToEdit);
  //virtual ~MultiCompPlotEditor();

protected:

  rosic::rsMultiBandCompressor* multiCompCore;
  jura::MultiCompAudioModule* multiCompModule; 

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompPlotEditor)
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
  MultiCompPlotEditor* plotEditor;
  RSlider *numBandsSlider;
  std::vector<RSlider*> splitFreqSliders, thresholdSliders, ratioSliders, attackSliders, 
    releaseSliders;
  //std::vector<RButton*> editButtons, soloButtons, muteButtons;
  RComboBox *splitModeBox, *bandSelectBox;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompModuleEditor)
};

#endif