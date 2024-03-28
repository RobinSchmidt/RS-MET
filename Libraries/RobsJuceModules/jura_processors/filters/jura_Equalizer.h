#ifndef jura_Equalizer_h
#define jura_Equalizer_h

/** This class wraps rosic::Equalizer into a rosof::AudioModule to facilitate its use as plugIn.

\todo: maybe handle the selection here also - carry over most of the setters/getters from 
EqualizerPlotEditor, handle ALL threading aspects here - get rid of mutexes in the GUI and DSP-core

-maybe write classes RepaintNotifier (as baseclass for EqualizerAudioModule) and 
RepaintNotifyObserver(as baseclass for the editor) in order to not automatically repaint on all 
ChangeBroadcast messages (this repaints unnecessarily often on preset switches) */

class EqualizerAudioModule : public ModulatableAudioModule, public ParameterSetHolder
{

  friend class EqualizerPlotEditor;
  friend class EqualizerModuleEditor;

public:


  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  EqualizerAudioModule(CriticalSection *newPlugInLock, 
    rosic::EqualizerStereo *equalizerStereoToWrap);

  EqualizerAudioModule(CriticalSection *newPlugInLock);

  void init();

  virtual ~EqualizerAudioModule();

  AudioModuleEditor* createEditor(int type) override;


  //-----------------------------------------------------------------------------------------------
  // automation and state management:

  /** Creates the static parameters for this module (i.e. parameters that are not created 
  dynamically and are thus always there). */
  virtual void createStaticParameters();

  // \todo: addParameter, removeParameter, createMetaParameters

  /** Overrides the callback that is called when one of the parameters has been changed. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  // \todo: introduce callback metaParameterChanged

  /** Returns the state of this module as XmlElement. */
  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  /** Restores the state of this module from an XmlElement (which was presumably previously created 
  via getStateAsXml). */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean) override;

  /** Converts a state which might possibly be from an older version to the current patch-format. */
  virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState) override;

  /** Returns the state of the given channel as XmlElement. */
  virtual XmlElement* getChannelStateAsXml(int channelIndex);

  /** Returns the state of the given band in the given channel as XmlElement. */
  virtual XmlElement* getBandStateAsXml(int channelIndex, int bandIndex);

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the sample-rate. */
  virtual void setSampleRate(double newSampleRate) override;

  /** Marks one of the channels as selected. */
  virtual void selectChannel(int channelToSelect);

  /** Marks one of the bands as selected and returns the index of the selected band (in case of 
  success, this will be 'indexToSelect', otherwise -1) */
  virtual int selectBand(int channel, int indexToSelect);

  /** De-selects the currently selected band (if any). */
  virtual void deSelectBand();

  /** Adds an equalizer band. The selectNewBand parameter detremines whether the newly added band 
  will also be selected, the suppressNotification argument suppresses the notification to our 
  ChangeListeners (presumably GUI elements) - these two arguments are useful when adding more than 
  one band in a loop such as in preset recall (in which case one does not want to select and notify 
  for each of the bands). 
  WARNING: when adding a band, memory re-allocations for our Parameter-arrays may occur. If some 
  outlying class keeps a pointer to one of the parameters in these arrays (such as class 
  EqualizerPlotEditor does), this pointer will become invalid - therefore, normally a notification 
  should occur (such that these classes can re-assign their pointers). If you suppress the 
  notification, you should make sure that outlying classes have invalidated their pointers before, 
  for example by calling removeAllBands before (which removes all parameters and then sends out a 
  notification which presumably causes our ParameterSetObservers to invalidate their pointers). */
  virtual void addBand(int channel, int mode, double frequency, double gain, 
    double bandwidth = 2.0*asinh(1.0/sqrt(2.0))/log(2.0), bool selectNewBand = true, 
    bool suppressNotification = false);

  /** Removes one of the bands and returns true when removal was successful. */
  virtual bool removeBand(int channel, int indexToRemove);

  /** Removes all of the bands. */
  virtual void removeAllBands();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of channels that should be drawn in the plot - in case of stereo-linked 
  and mono, this will be one, else two. */
  virtual int getNumChannelsToPlot();

  /** Returns the total number of bands in the underlying equalizer. */
  virtual int getNumBands(int channel);

  /** Returns the index of the currently selected band (or -1, if none is selected). */
  virtual int getSelectedBandIndex() { return selectedIndex; }

  /** Writes the magnitude response at the given frequencies into the passed array. */
  virtual void getMagnitudeResponse(int channel, double *frequencies, double *magnitudes, 
    int numBins);
 
  /** Returns the magnitude response expressed in decibels of the given channel at the given 
  frequency. */  
  double getDecibelsAt(int channel, double frequency);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    for(int n = 0; n < numSamples; n++)
      wrappedEqualizerStereo->getSampleFrameStereo(&inOutBuffer[0][n], &inOutBuffer[1][n]);
  }

  virtual void processStereoFrame(double *left, double *right) override
  {
    wrappedEqualizerStereo->getSampleFrameStereo(left, right);
  }

  //-----------------------------------------------------------------------------------------------
  // others:

  virtual void reset() override;

protected:

  /** Re-assigns the callback-objects inside the dynamically created Parameter objects - this must 
  be called whenever the number of bands was changed since such a change may result in 
  re-allocation of all the objects that represent the bands in the underlying core DSP module. */
  virtual void assignCallbacksForDynamicParameters();

  rosic::EqualizerStereo *wrappedEqualizerStereo;
  bool wrappedEqualizerIsOwned = false;


  int selectedChannel; // currently selected channel (0 for L or M, 1 for R or S)
  int selectedIndex;   // index of currently selected band (-1 if none)

  // arrays for the per-band parameters (later: ModulatableParameter):
  juce::OwnedArray<Parameter> filterModeParameters[2];
  juce::OwnedArray<Parameter> frequencyParameters[2];
  juce::OwnedArray<Parameter> gainParameters[2];
  juce::OwnedArray<Parameter> bandwidthParameters[2];


  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/**

This class plots the frequency responses of a rosic::Equalizer object and allows for editing 
parameters like the center frequencies. 

\todo right-click on empty area shows menu to insert different kinds of bands

*/

class EqualizerPlotEditor	: virtual public rsSpectrumPlot, public ParameterObserver, 
  public ParameterSetObserver 
{

  friend class EqualizerModuleEditor;  // we need this to reach through to the equalizerModuleLock from the editor

  /** Enumeration of the handles that can be grabbed and dragged by the mouse.  */
  enum dragHandles
  {
    NONE = 0,
    FREQUENCY_AND_GAIN,
    BANDWIDTH_AND_GAIN_LEFT,
    BANDWIDTH_AND_GAIN_RIGHT,
    GLOBALGAIN_LINE,
    GAIN, 
    FREQUENCY
  };

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. You must pass a pointer an EqualizerAudioModule object which is to be edited and 
  a pointer to a CriticalSection which will be used to wrap all accesses to this pointer - 
  typically, this CriticalSection should be the plugInLock which is member of the underlying 
  AudioPlugIn object. */  
  EqualizerPlotEditor(CriticalSection *newPlugInLock, 
    EqualizerAudioModule* newEqualizerModuleToEdit);

  /** Destructor. */
  virtual ~EqualizerPlotEditor(); 

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the EqualizerAudioModule object which is to be edited. */
  virtual void setEqualizerModuleToEdit(EqualizerAudioModule* newEqualizerModuleToEdit);

  /** Un-assigns our pointer-members for the parameters and resets them to NULL. */
  virtual void unAssignParameters();

  /** Assigns our pointer-members for the parameters to the parameter that is currently selected 
  in the underlying EqualizerAudioModule. */
  virtual void assignParametersToSelectedBand();


  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the index of the band the is represented by a dot at the given pixel position. It 
  will return -1 if there isn't any band at the position in question. */
  virtual int getBandIndexAtPixelPosition(int x, int y);
  //virtual int getChannelAndIndexAtPixelPosition(int x, int y, int &channel, int& index);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides the changeListenerCallback to update the plot. */
  //virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  // we should get rid of this - we now update ourselves in parameterChanged


  /** Overrides the purely virtual parameterSetChanged() method of the ParameterSetObserver base 
  class. */
  virtual void parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged) override;

  /** Overrides the purely virtual parameterChanged() method of the ParameterObserver base class. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;

  /** Overrides the purely virtual method of the ParameterObserver base class in order to 
  invalidate our pointer-member 'assignedParameter'. */
  virtual void parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted) override;

  /** Overrides mouseMove in order to update the cursor according to what is under the mouse. */
  virtual void mouseMove(const MouseEvent &e) override;

  /** Overrides mouseDown for adjusting the frequency and resonance and lets a context menu pop up 
  when the right button is clicked for MIDI-learn functionality. */
  virtual void mouseDown(const MouseEvent& e) override;

  /** Overrides mouseDrag for adjusting the frequency and resonance. */
  virtual void mouseDrag(const MouseEvent& e) override;

  /** Overrides mouseUp to reset the currentDragHandle to NONE. */
  virtual void mouseUp(const MouseEvent& e) override;

  /** Overrides mouseWheelMove to adjust the bandwidth on wheel moves. */
  virtual void mouseWheelMove (const MouseEvent& ev, const MouseWheelDetails& wheel) override;

  /** Overrides the resized-method. */
  virtual void resized() override;

  /** Updates the frequency response plot. */
  virtual void updatePlot();


protected:

  /** Returns the handle for mouse grab/drag under the specified position (in pixels) as one of 
  the values in enum dragHandles. */
  virtual int getDragHandleAt(int x, int y);

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the handles. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL, 
    XmlElement *targetSVG = NULL) override;

  /** Creates the Popup menu that is used by openRightClickPopupMenu(). */
  //virtual void createRightClickPopupMenu(PopupMenu*& menu);

  /** Handles the result of opening the right-click popoup menu. */
  //virtual void handleRightClickPopupMenuResult(int result, int x, int y);

  /** Opens the PopupMenu that appears on right clicks. */
  //virtual void openRightClickPopupMenu(int x, int y);

  /** Converts a pixel-position to the corresponing frequency and gain values. */
  virtual void xyToFrequencyAndGain(double &x, double &y);

  CriticalSection      *plugInLock;              // mutex to access the edited AudioModule object 
  EqualizerAudioModule *equalizerModuleToEdit;   // the edited AudioModule object

  int currentlyDraggedHandle;  // kind of the handle that is currently being dragged

  // the parameters which wil cause re-plotting and therefore must be listened to:
  Parameter *filterModeParameter, *frequencyParameter, *gainParameter, *bandwidthParameter, 
    *globalGainParameter;

  // magnitude response display stuff:
  int    numBins;


  //double *frequencies, *magnitudes1, *magnitudes2;   // OLD

  std::vector<double> frequencies, magnitudes1, magnitudes2;  // NEW

  double *magnitudes[2];  // Needed for call to setSpectra()

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** New version - not yet used. */

class rsEqualizerNodeEditor : public rsNodeEditor
{

public:

  //rsEqualizerNodeEditor(EqualizerAudioModule* eqModule);
  //rsEqualizerNodeEditor();

protected:

  //EqualizerAudioModule* equalizerModule;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsEqualizerNodeEditor)
};

/** New version of EqualizerPlotEditor based on new plotting facilities and rsNodeEditor.*/

class rsEqualizerPlotEditor : public ColourSchemeComponent, public rsNodeEditorObserver
{

public:

  rsEqualizerPlotEditor(EqualizerAudioModule* eqModule);

  // overrides:
  virtual void nodeWasAdded(rsNodeEditor* editor, int nodeIndex) override;
  virtual void nodeWillBeRemoved(rsNodeEditor* editor, int nodeIndex) override;
  virtual void nodeWasMoved(rsNodeEditor* editor, int nodeIndex) override;

protected:


  EqualizerAudioModule* equalizerModule;
  rsFunctionPlot*       freqRespPlot;
  rsNodeEditor*         nodeEditor;
  //rsEqualizerNodeEditor* eqNodeEditor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsEqualizerPlotEditor)
};

// maybe we need a class rsNodeEditorObserver with callbacks nodeChanged(int index), 
// nodeWasAdded(int index), nodeWillBeRemoved(int index) then derive rsEqualizerPlotEditor from
// rsNodeEditorObserver and have a rsNodeEditor member that is being observed

//=================================================================================================

class EqualizerModuleEditor : public AudioModuleEditor, public RComboBoxObserver, 
  public ParameterSetObserver 
{

public:

  enum layouts
  {
    SLIDERS_RIGHT,
    SLIDERS_BELOW,
    SLIDERS_ABOVE
    //FOR_QUADRIFEX
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  EqualizerModuleEditor(CriticalSection *newPlugInLock, 
    EqualizerAudioModule* newEqualizerAudioModule);

  virtual ~EqualizerModuleEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the EqualizerAudioModule object which is to be edited. */
  virtual void setEqualizerModuleToEdit(EqualizerAudioModule* newEqualizerModuleToEdit);

  virtual void setLayout(int newLayout) { layout = newLayout; }

  virtual void setUseShortSliderNames(bool shouldBeShort);
  virtual void setUseSmallComboBox(   bool shouldBeSmall);

  //virtual void updateDynamicParametersToWidgets();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox  *rComboBoxThatHasChanged);
  //virtual void rSliderValueChanged(RSlider *rSliderThatHasChanged);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged);
  virtual void paint(Graphics &g);
  virtual void resized();

  virtual void updateWidgetsAccordingToState();



protected:

  void createWidgets();
  virtual void updateWidgetVisibility();
  virtual void updateWidgetAppearance(); 
  virtual void updatePlotRange();

  //virtual void showGainScaleSlider(); ...etc.

  EqualizerAudioModule* equalizerModule;

  EqualizerPlotEditor* plotEditor;

  // rectangles for organizing the gui (can be done by baseclass later...)
  juce::Rectangle<int> rightSectionRectangle, bottomSectionRectangle; 

  // widgets:
  RTextField        *bandParametersLabel;
  rsAutomatableButton *bypassButton;
  RButton           *copyButton, *pasteButton, *invertButton;
  RNamedComboBox    *stereoModeComboBox, *gainRangeComboBox, *filterModeComboBox;
  rsModulatableSlider *frequencySlider, *gainSlider, *bandwidthSlider, *globalGainSlider;
  RRadioButton      *channelSelectButton1, *channelSelectButton2;
  RRadioButtonGroup channelSelectRadioGroup;

  int layout;

  bool useShortSliderNames, useSmallComboBox;

  juce_UseDebuggingNewOperator;
};

#endif 
