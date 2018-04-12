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

/** A general multiband effect in which the type of the effect which should be applied to each band 
can be selected. ...tbc

todo: 
-when the rightmost band is selected, it shows the gui for Compressor1
-split-freq sliders are not shown
-check that the required split-freq parameters and their widgets are handled correctly
-implement state recall
-implement switching the type of effect
-plot frequency responses: rosic::rsMultiBandEffect::getMagnitudeAt(index, freq) */

class JUCE_API MultiBandEffect : public jura::ModulatableAudioModule, public ChangeBroadcaster
{

public:

  MultiBandEffect(CriticalSection *lockToUse, MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);
  virtual ~MultiBandEffect();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  AudioModuleEditor* createEditor() override;

  virtual void parameterChanged(Parameter* p) override;

  /** Creates the parameters related to the band-splitting. Called from setEffectCore. */
  //virtual void createSplittingParameters();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  virtual void registerMultiBandObserver(MultiBandEffectObserver *obs)
  { 
    appendIfNotAlreadyThere(observers, obs); 
  }

  virtual void deRegisterMultiBandObserver(MultiBandEffectObserver *obs)
  { 
    removeFirstOccurrence(observers, obs); 
  }

  /** Sets the type of the effect that shall be applied to each band. The argument should be one of
  the strings returned by getAvailableEffectTypes(). */
  void setEffectType(const juce::String& typeString); 


  virtual void insertBand(int index, double splitFrequency, bool sendNotification); 
  virtual void removeBand(int index, bool mergeWithRightNeighbour, bool sendNotification);
  // ...must rename split-freq parameters and re-assign their callbacks
  virtual void selectBand(int bandToSelect, bool sendNotification);

  void createSplitFreqParams();
  void addSplitFreqParam(int index, double freq);
  void removeSplitFreqParam(int index);


  /** Sets the splitting frequency for the band with given index to the given frequency. */
  void setSplitFreq(int bandIndex, double newFreq);

  /** Sets the splitting frequency for the selected band to the given frequency. */
  void setSplitFreq(double newFreq);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns a pointer to the parameter for the splitting frequency with given index. */
  Parameter* getSplitFreqParam(int bandIndex);

  /** Returns the upper cutoff frequency for the band with given index. */
  double getSplitFreq(int bandIndex) const { return core.getSplitFrequency(bandIndex); }


  int getSelectedBand() const { return selectedBand; }
  int getNumBands()     const { return core.getNumberOfBands(); }
  int getBandContainingFrequency(double freq);


  /** Returns an array of strings of names of the available effect types. The elements of this 
  array can be used as arguments for setEffectType. */
  std::vector<juce::String> getAvailableEffectTypes() const;

  /** Returns a pointer to our rosic::rsMultiBandEffect object that is responsible for the 
  bandsplitting. */
  rosic::rsMultiBandEffect* getCore() { return &core; }

  /** Returns a pointer to the effect module for the given band. */
  jura::AudioModule* getBandEffect(int index) const { return perBandModules[index]; }

  /** Returns true, if the splitting frequencies between the bands are in (strictly) increasing 
  order. Used for debugging. */
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

  /** Removes all per-band effects. */
  void clearBandEffects();

  /** Makes sure that the module for each band has the correct moduleName. */
  void updateBandModuleNames();


  rosic::rsMultiBandEffect core;

  int selectedBand =  -1;  // -1 is code for "None"
  std::vector<Parameter*> splitFreqParams; 
  std::vector<MultiBandEffectObserver*> observers;

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

/** Editor for a generic mutiband effect. */

class JUCE_API MultiBandEffectEditor : public AudioModuleEditor, public MultiBandEffectObserver,
  public RComboBoxObserver
{

public:

  MultiBandEffectEditor(MultiBandEffect* effect);
  virtual ~MultiBandEffectEditor();

  virtual void resized() override;
  virtual void rComboBoxChanged(RComboBox* comboBoxThatHasChanged) override;

  virtual void bandWasInserted(MultiBandEffect* mbe, int index) override;
  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) override;
  virtual void bandWasSelected(MultiBandEffect* mbe, int index) override;


protected:

  /** Inserts a new per band effect editor for given index. */
  virtual void insertBandEditor(int index);

  /** Removes the per band effect editor at given index. */
  virtual void removeBandEditor(int index);

  /** Makes currently required sub-editor not required editors invisible. */
  virtual void updateEditorVisibility();

  /** Creates the global (not per-band) widgets. */
  virtual void createWidgets();

  /** Creates the per-band sub editors, if necessarry. */
  virtual void createBandEditors();

  /** Deletes all the per-band sub editors. */
  virtual void clearBandEditors();

  /** Sets up positions of all per band editors. */
  virtual void positionBandEditors();

  /** Sets up position of the band effect editor with given index. */
  virtual void positionBandEditor(int index);

  // widgets:
  MultiBandPlotEditor* plotEditor;
  RComboBox *effectSelectBox, *splitModeBox;
  std::vector<RSlider*> splitFreqSliders;


  MultiBandEffect* effectToEdit;                  // pointer to the edited multiband effect
  std::vector<AudioModuleEditor*> perBandEditors; // sub editor array (one editor for each band)


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandEffectEditor)
};

















// code below is obsolete, but we may want to do a dedicated MultiBandCompressor in which the 
// threshold/ratio parameters can be adjusted from the plot

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