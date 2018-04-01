#ifndef jura_MultiComp_h
#define jura_MultiComp_h


class MultiBandEffect;

/** Baseclass for classes that must keep track of the state of a MultiBandEffect object. */

class JUCE_API MultiBandEffectObserver
{

public:

  virtual void bandWasInserted(MultiBandEffect* mbe, int index) = 0;

  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) = 0;

  virtual void bandWasSelected(MultiBandEffect* mbe, int index) = 0;

};

//=================================================================================================

/** Baseclass for multiband effects. 

todo: 

-ensure that the split frequencies are always sorted from low to high
-restrict ranges for the split-freqs according to the neighbours
 (maybe do these things in rosic::MultiBandEffect)
OR: forget about that sorting and:
-let the user add remove bands via right-clcik context menu with items:
 -add: splits the band inside which the click occurred into two, settings are copied into the
  new band
 -remove: deletes the band in which the click occurred - either the left or the right neighbour
  can cover the range where the band was...so maybe have: merge with left / merge with right
 -the NumBands slider is useless then
 -the band selection menu is superfluous - maybe we should just display and ifo such as
  Band: 2/5 when the 2nd of 5 bands is selected
 -we somehow need the ability to insert/remove bands in rosic::rsMultiBandEffect
  void insertBand(int indexOfLeftNeighbour), void removeBand(int index, bool mergeWithRightNeighbour)
  

-plot frequency responses: rosic::rsMultiBandEffect::getMagnitudeAt(index, freq)


*/

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


  virtual void registerMultiBandObserver(MultiBandEffectObserver *obs)
  {
    appendIfNotAlreadyThere(observers, obs);
  }

  virtual void deRegisterMultiBandObserver(MultiBandEffectObserver *obs)
  {
    removeFirstOccurrence(observers, obs);
  }

  virtual void parameterChanged(Parameter* p) override;

  void selectBand(int bandToSelect) { selectedBand = bandToSelect; }


  /** Sets the number of (active) bands. */
  //void setNumBands(int newNumBands);


  void insertBand(int index, double splitFrequency); 

  void removeBand(int index, bool mergeWithRightNeighbour = false);


  /** Sets the splitting frequency for the band with given index to the given frequency. */
  void setSplitFreq(int bandIndex, double newFreq);

  /** Sets the splitting frequency for the selected band to the given frequency. */
  void setSplitFreq(double newFreq);


  /** Returns a pointer to the parameter for the splitting frequency with given index. */
  Parameter* getSplitFreqParam(int bandIndex);


  int getSelectedBand() const { return selectedBand; }
  int getNumBands()     const { return core->getNumberOfBands(); }
  int getMaxNumBands()  const { return maxNumBands; }

  int getBandContainingFrequency(double freq);



  /** Returns a pointer to our core DSP object. */
  rosic::rsMultiBandEffect* getCore() { return core; }

  /** Returns true, if the splitting frequencies between the bands are in (strictly) increasing 
  order. */
  bool areBandsInIncreasingOrder(bool strictly = false);

protected:

  void sendBandInsertNotification(int index);
  void sendBandRemoveNotification(int index);
  void sendBandSelectNotification(int index);



  rosic::rsMultiBandEffect* core = nullptr;

  int maxNumBands  =  0;  // assigned in constructor
  int selectedBand =  0;  // -1 is code for "None"

  std::vector<Parameter*> splitFreqParams;

  std::vector<MultiBandEffectObserver*> observers;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandEffect)
};

//=================================================================================================

/** Plot editor for multiband effects. Allows to select a band by clicking into the rectangular 
area in the frequency response plot that represents that band, change split frequencies by dragging
vertical line, etc. */

class JUCE_API MultiBandPlotEditor : public ColourSchemeComponent, public ChangeListener, 
  public RPopUpMenuObserver
{

public:



  MultiBandPlotEditor(jura::MultiBandEffect* multiBandEffectModuleToEdit);
  virtual ~MultiBandPlotEditor();

  //virtual void parameterChanged(Parameter* p) override;
  virtual void changeListenerCallback(ChangeBroadcaster* source) override;
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void paintOverChildren(Graphics& g) override;
  virtual void resized() override;

protected:

  virtual void openRightClickMenu();

  rosic::rsMultiBandEffect* core;
  jura::MultiBandEffect* module; 

  rsFunctionPlot* freqRespPlot;

  RPopUpMenu *bandPopup = nullptr; // created when needed the first time
  enum popUpIds
  {
    ADD_BAND = 1,
    REMOVE_BAND
  };



  double freqAtMouse = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandPlotEditor)
};

//=================================================================================================

/** Multiband compressor with up to 16 bands. */

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
  //RSlider *numBandsSlider; // remove
  std::vector<RSlider*> splitFreqSliders, thresholdSliders, ratioSliders, attackSliders, 
    releaseSliders;

  //RButton *addBandButton, *removeBandButton;

  //std::vector<RButton*> editButtons, soloButtons, muteButtons;
  RComboBox *splitModeBox, *bandSelectBox;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiCompModuleEditor)
};

#endif