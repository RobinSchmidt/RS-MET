#ifndef jura_MultiComp_h
#define jura_MultiComp_h


class MultiBandEffect;

/** Baseclass for classes that must keep track of the state of a MultiBandEffect object. */

class JUCE_API MultiBandEffectObserver
{

public:

  /** Called after a new band was inserted. */
  virtual void bandWasInserted(MultiBandEffect* mbe, int index) = 0;

  /** Called before a band will be removed. */
  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) = 0;

  /** Called after another band was selected. */
  virtual void bandWasSelected(MultiBandEffect* mbe, int index) = 0;

  /** Called after the state was completely changed, for example due to preset recall. */
  //virtual void totalRefreshNeeded(MultiBandEffect* mbe);

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

  /*
  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;
  */



  /** Subclasses should call this *once* in their constructor with a pointer to the concrete 
  multiband effect (subclass). This will also create the splitting-related parameters */
  void setEffectCore( rosic::rsMultiBandEffect* effectCore);
    // obsolete

  /** Sets the type of the effect that shall be applied to each band. The argument should be one of
  the strings returned by getAvailableEffectTypes(). */
  void setEffectType(const juce::String& typeString); 

  /** Creates the parameters related to the band-splitting. Called from setEffectCore. */
  //virtual void createSplittingParameters();


  virtual void registerMultiBandObserver(MultiBandEffectObserver *obs)
  {
    appendIfNotAlreadyThere(observers, obs);
  }

  virtual void deRegisterMultiBandObserver(MultiBandEffectObserver *obs)
  {
    removeFirstOccurrence(observers, obs);
  }

  virtual void parameterChanged(Parameter* p) override;


  virtual void insertBand(int index, double splitFrequency, bool sendNotification); 
  virtual void removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification);
  // ...must rename split-freq parameters and re-assign their callbacks


  virtual void selectBand(int bandToSelect, bool sendNotification);

  // get rid of these (handle in subclass):
  void createSplitFreqParams();
  void addSplitFreqParam(int index, double freq);
  void removeSplitFreqParam(int index);


  /** Sets the splitting frequency for the band with given index to the given frequency. */
  void setSplitFreq(int bandIndex, double newFreq);

  /** Sets the splitting frequency for the selected band to the given frequency. */
  void setSplitFreq(double newFreq);


  /** Returns a pointer to the parameter for the splitting frequency with given index. */
  Parameter* getSplitFreqParam(int bandIndex);
    // make purely virtual, override in MultiComp

  /** Returns the upper cutoff frequency for the band with given index. */
  double getSplitFreq(int bandIndex) const { return core->getSplitFrequency(bandIndex); }


  int getSelectedBand() const { return selectedBand; }
  int getNumBands()     const { return core->getNumberOfBands(); }
  int getBandContainingFrequency(double freq);


  /** Returns an array of strings of names of the available effect types. The elements of this 
  array can be used as arguments for setEffectType. */
  std::vector<juce::String> getAvailableEffectTypes();

  /** Returns a pointer to our core DSP object. */
  rosic::rsMultiBandEffect* getCore() { return core; }


  // debug functions:

  /** Returns true, if the splitting frequencies between the bands are in (strictly) increasing 
  order. */
  bool areBandsInIncreasingOrder(bool strictly = false);



protected:

  void sendBandInsertNotification(int index);
  void sendBandRemoveNotification(int index);
  void sendBandSelectNotification(int index);

  /** Inserts a new per band effect module to the array at the given index. */
  void insertBandEffect(int index);

  /** Removes a per band effect module from the array at the givne index (and deletes the 
  object). */
  void removeBandEffect(int index);


  int selectedBand =  -1;  // -1 is code for "None"
  std::vector<Parameter*> splitFreqParams; 
  std::vector<MultiBandEffectObserver*> observers;

  // todo: don't use this anymore:
  rosic::rsMultiBandEffect* core = nullptr; // ...or have a direct member and use only the splitting
                                            // functionality

  // ...instead use that:
  AudioModuleFactory perBandModuleFactory;
  std::vector<AudioModule*> perBandModules;
  juce::String effectTypeString;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandEffect)
};

//=================================================================================================

/** Plot editor for multiband effects. Allows to select a band by clicking into the rectangular 
area in the frequency response plot that represents that band, change split frequencies by dragging
vertical line, etc. */

class JUCE_API MultiBandPlotEditor : public ColourSchemeComponent, public ChangeListener, // ChangeListener obsolete?
  public RPopUpMenuObserver,  public MultiBandEffectObserver
{

public:



  MultiBandPlotEditor(jura::MultiBandEffect* multiBandEffectModuleToEdit);
  virtual ~MultiBandPlotEditor();



  //virtual void parameterChanged(Parameter* p) override;
  virtual void bandWasInserted(MultiBandEffect* mbe, int index) override;
  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) override;
  virtual void bandWasSelected(MultiBandEffect* mbe, int index) override;

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

/** Experimental datastructure to hold the parameters of a single band of a multiband 
compressor.*/

class JUCE_API BandCompParameterSet // public BandParameterSet (todo:...has only freq param)
{

public:

  BandCompParameterSet(rosic::rsMultiBandCompressor* core);
  virtual ~BandCompParameterSet();

  virtual void setBandIndex(int newIndex);
   // should update bandIndex, names of the parameters and callbacks

  int getBandIndex() const { return bandIndex; }


  // callback target functions:
  void setSplitFreq(double newFreq)    { compCore->setSplitFrequency(bandIndex, newFreq);    }
  void setThreshold(double newThresh)  { compCore->setThreshold(     bandIndex, newThresh);  }
  void setRatio(    double newRatio)   { compCore->setRatio(         bandIndex, newRatio);   } 
  void setAttack(   double newAttack)  { compCore->setAttackTime(    bandIndex, newAttack);  } 
  void setRelease(  double newRelease) { compCore->setReleaseTime(   bandIndex, newRelease); } 

  //void dummyCallback(double newValue) {} // used when bandIndex is set to -1

protected:

  Parameter *freq, *thresh, *ratio, *att, *rel; 

  int bandIndex = -1;
  rosic::rsMultiBandCompressor* compCore;

};

//=================================================================================================

/** Multiband compressor with up to 16 bands. */

class JUCE_API MultiCompAudioModule : public jura::MultiBandEffect
{

public:

  MultiCompAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);


  virtual void insertBand(int index, double splitFrequency, bool sendNotification) override; 
  virtual void removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification) override;
    // ...must rename band compression parameters and re-assign callbacks


  void createBandParams();
  void addCompressionParams(int index);
  void removeCompressionParams(int index);

  /** Returns a pointer to our core DSP object. */
  rosic::rsMultiBandCompressor* getMultiCompCore() { return &multiCompCore; }

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;

  // target functions for the per-band parameter callbacks (the "if"s are kludgy):
  void setThreshold(double newThreshold) { if(selectedBand >= 0) multiCompCore.setThreshold(     selectedBand, newThreshold); }
  void setRatio(    double newRatio)     { if(selectedBand >= 0) multiCompCore.setRatio(         selectedBand, newRatio);     }
  void setAttack(   double newAttack)    { if(selectedBand >= 0) multiCompCore.setAttackTime(    selectedBand, newAttack);    }
  void setRelease(  double newRelease)   { if(selectedBand >= 0) multiCompCore.setReleaseTime(   selectedBand, newRelease);   }

protected:

  //bool checkParameterNames

  rosic::rsMultiBandCompressor multiCompCore;

  std::vector<BandCompParameterSet*> bandParamSets; // not yet used

  int numCompParamSets = 0;

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

class MultiCompModuleEditor : public AudioModuleEditor, public RComboBoxObserver, 
  public MultiBandEffectObserver
{

public:

  MultiCompModuleEditor(MultiCompAudioModule* multiCompModuleToEdit);
  virtual ~MultiCompModuleEditor();

  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;

  virtual void bandWasInserted(MultiBandEffect* mbe, int index) override;
  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) override;
  virtual void bandWasSelected(MultiBandEffect* mbe, int index) override;



  virtual void addBandWidgets(int index);
  virtual void removeBandWidgets(int index);

  virtual void resized() override;

protected:


  virtual void createWidgets();

  /** Creates the per-band widgets, if necessarry. */
  virtual void createBandWidgets();



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