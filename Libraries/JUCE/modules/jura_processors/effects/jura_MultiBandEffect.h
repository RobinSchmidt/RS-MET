#ifndef jura_MultiComp_h
#define jura_MultiComp_h

// todo: rename files to jura_ MultiBandEffect.h/cpp

class MultiBandEffect;

/** Baseclass for classes that must keep track of the state of a MultiBandEffect object. */

class JUCE_API MultiBandEffectObserver
{

public:

  /** Called after a new band was inserted. */
  virtual void bandWasInserted(MultiBandEffect* mbe, int index) = 0;

  /** Called before a band will be removed. */
  virtual void bandWillBeRemoved(MultiBandEffect* mbe, int index) = 0;

  /** Called after a band has been removed. */
  virtual void bandWasRemoved(MultiBandEffect* mbe, int index) = 0;

  /** Called after another band was selected. */
  virtual void bandWasSelected(MultiBandEffect* mbe, int index) = 0;

  /** Called before all of the band effects and splitters will be removed. An observer should 
  invalidate all its references to per-band effects and split-frequency parameters and possibly 
  even delete associated objects (like editors and widgets). The main purpose is to empty the 
  editor before a new preset is being loaded. */
  virtual void allBandsWillBeRemoved(MultiBandEffect* mbe) = 0;

  /** Called after the state was completely changed, for example after a preset recall. */
  virtual void totalRefreshNeeded(MultiBandEffect* mbe) = 0;

};

//=================================================================================================

/** A general multiband effect in which the type of the effect which should be applied to each band 
can be selected. ...tbc

todo: 
-BUG: create 4 bands, remove 1st - there are now 3 bands but the 1st has wrong split freq (even 
 though its split-freq slider is correct)...i think, when removing a band, we also need to update
 the split-frequencies in the audio engine ...seems fixed
-BUG: crash when removing the very last band
-i think, we need to allow 0 bands in the core - this will be interpreted a bypass

-we may need to acquire the lock all the editor member functions
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
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

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
  virtual void selectBand(int bandToSelect, bool sendNotification);


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

  /** Returns index of the band that is currently selected (-1, if none). */
  int getSelectedBand() const { return selectedBand; }

  /** Returns the number of bands. This can be zero in which case there is no effect at all. If it 
  is one, we have have a regular, single broadband effect. */
  int getNumBands() const;

  /** Returns the number of frequency splitters which is max(0, getNumBands()-1). */
  int getNumSplits() const;

  /** Returns the band that contains the given frequency (assumes ordered, non-overlapping 
  bands ...todo: somehow it should work reasonably also with unordered and/or overlapping 
  bands). */
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

  /** Returns true, if the band is given index is removable (for technical reasons, the rightmost
  band is not removable - todo: lift that restriction). */
  bool isBandRemovable(int index);


protected:

  // parameter management:
  void createSplitFreqParams();
  void clearSplitFreqParams();
  void addSplitFreqParam(int index, double freq);
  void removeSplitFreqParam(int index);
  void updateSplitFreqParamNamesAndCallbacks();

  // notification senders:
  void sendBandInsertNotification(int index);
  void sendBandRemovePreNotification(int index);
  void sendBandRemovePostNotification(int index);
  void sendBandSelectNotification(int index);
  void sendClearBandsNotification();
  void sendTotalRefreshNotification();

  /** Inserts a new per band effect module to the array at the given index. */
  void insertBandEffect(int index);

  /** Removes a per band effect module from the array at the givne index (and deletes the 
  object). */
  void removeBandEffect(int index);

  /** Removes all per-band effects. */
  void clearBandEffects(bool callRemoveNotification);

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
vertical line, etc. 

todo:
-plot frequency curves
-plot gain reduction
-maybe it needs to connect to the split-freq parameters
-maybe instead of embedding a rsFunctionPlot object, we should subcall from it, so we may directly
 override paint instead of paintOverChildren
*/

class JUCE_API MultiBandPlotEditor : public ColourSchemeComponent, public ChangeListener, // ChangeListener obsolete?
  public RPopUpMenuObserver,  public MultiBandEffectObserver
{

public:

  MultiBandPlotEditor(jura::MultiBandEffect* multiBandEffectModuleToEdit);
  virtual ~MultiBandPlotEditor();

  //virtual void parameterChanged(Parameter* p) override;

  // MultiBandEffectObserver overrides:
  virtual void bandWasInserted(      MultiBandEffect* mbe, int index) override;
  virtual void bandWillBeRemoved(    MultiBandEffect* mbe, int index) override;
  virtual void bandWasRemoved(       MultiBandEffect* mbe, int index) override;
  virtual void bandWasSelected(      MultiBandEffect* mbe, int index) override;
  virtual void allBandsWillBeRemoved(MultiBandEffect* mbe) override;
  virtual void totalRefreshNeeded(   MultiBandEffect* mbe) override;

  // other callbacks:
  virtual void changeListenerCallback(ChangeBroadcaster* source) override;
  virtual void rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void resized() override;
  virtual void paintOverChildren(Graphics& g) override;

protected:


  /** Paints the shadings of the various bands (currently to highlight the selected band - maybe 
  later more stuff) */
  void paintBandShadings(Graphics& g);

  /** Paints the vertical lines at the split frequencies. */
  void paintSplitLines(Graphics& g);

  /** Opens the right-click context menu for inserting and removing bands. It takes as parameter 
  the index of the band inside of which the mouse-click occured (required because available options 
  may be different for different bands) */
  virtual void openRightClickMenu(int bandIndex);

  /** Refreshes the content of our freqRespPlot member, i.e. the frequency response curves. */
  virtual void refreshFunctionsToPlot();


  rosic::rsMultiBandEffect* core;
  jura::MultiBandEffect* module; 

  rsFunctionPlot* freqRespPlot;

  RPopUpMenu *bandPopup = nullptr; // created when needed the first time
  enum popUpIds
  {
    ADD_BAND = 1,
    REMOVE_BAND
  };

  double freqAtMouse = 0; // why is this a member?

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandPlotEditor)
};

//=================================================================================================

/** Extends MultiBandPlotEditor to provide animated feedback such as drawing gain-reduction into
the plot (useful for multiband dynamics). */

class JUCE_API MultiBandPlotEditorAnimated : public MultiBandPlotEditor/*, public juce::Timer*/
{

public:

  MultiBandPlotEditorAnimated(jura::MultiBandEffect* moduleToEdit);
  virtual ~MultiBandPlotEditorAnimated();

  virtual void paintOverChildren(Graphics& g) override;

protected:

  /** Paints the gain of the output signal with respect to the input signal for each band. 
  Especially useful for dynamics processors. */
  void paintInOutGains(Graphics& g);


  // todo: somehow, we should buffer the static parts in the background into an image to avoid
  // redrawing them for each frame (which is expensive)


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

  // MultiBandEffectObserver overrides:
  virtual void bandWasInserted(      MultiBandEffect* mbe, int index) override;
  virtual void bandWillBeRemoved(    MultiBandEffect* mbe, int index) override;
  virtual void bandWasRemoved(       MultiBandEffect* mbe, int index) override;
  virtual void bandWasSelected(      MultiBandEffect* mbe, int index) override;
  virtual void allBandsWillBeRemoved(MultiBandEffect* mbe) override;
  virtual void totalRefreshNeeded(   MultiBandEffect* mbe) override;


protected:

  /** Inserts a new per band effect editor for given index. */
  virtual void insertBandEditor(int index);

  /** Removes the per band effect editor at given index. */
  virtual void removeBandEditor(int index);

  /** Makes currently required sub-editor not required editors invisible. */
  virtual void updateEditorVisibility();

  /** Updates the names of the editors according to the module names (the names may change at 
  runtime due to addition and removal of bands). */
  virtual void updateEditorNames();

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

  /** Updates the array of our split-frequency sliders, possibly adding or removing some, 
  connecting the to their parameters and setting up their positions. */
  virtual void updateSplitSliders();

  /** Updates the bounds of the sliders for the split frequencies. */
  virtual void updateSplitSliderPositions();

  // widgets:
  //MultiBandPlotEditor* plotEditor;
  MultiBandPlotEditorAnimated* plotEditor;
  RComboBox *effectSelectBox, *splitModeBox;
  std::vector<RSlider*> splitFreqSliders;         // sliders for splitting frequencies
  std::vector<AudioModuleEditor*> perBandEditors; // sub editor array (one editor for each band)
  MultiBandEffect* effectToEdit;                  // pointer to the edited multiband effect

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MultiBandEffectEditor)
};

#endif