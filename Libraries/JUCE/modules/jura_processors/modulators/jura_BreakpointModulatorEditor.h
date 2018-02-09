#ifndef jura_BreakpointModulatorEditor_h
#define jura_BreakpointModulatorEditor_h

//=================================================================================================
// class BreakpointModulatorGlobalEditor

/** This class encapsulates the widgets for the global settings that are used in a 
BreakpointModulatorEditor in order to facilitate working with arrays of such sets as needed in 
BreakpointModulatorEditorMulti. */

class JUCE_API BreakpointModulatorGlobalEditor : public AudioModuleEditor
{

  friend class BreakpointModulatorEditor;
  friend class BreakpointModulatorEditorMulti;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. \todo remove the CriticalSection parameter */
  BreakpointModulatorGlobalEditor(CriticalSection *newPlugInLock, 
    BreakpointModulatorAudioModule* newModulatorToEdit);

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::Modulator object which is to be edited. Make sure to 
  call this function again with a NULL-pointer when the object get deleted for some reason. */
  virtual void setModulatorToEdit(BreakpointModulatorAudioModule* newModulatorToEdit);

  /** Selects one of the available widget-layouts. 
  \todo: give them names in an enumeration  */
  virtual void setLayout(int newLayout);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  // pointer to the edited object:
  BreakpointModulatorAudioModule* modulatorToEdit;

  // widgets:
  AutomatableSlider *timeScaleSlider, *timeScaleByKeySlider, *timeScaleByVelSlider,
    *depthSlider, *depthByKeySlider, *depthByVelSlider;
  //AutomatableButton *loopButton, *syncButton; // later - make them automatable
  RButton *loopButton, *syncButton;
  RButton *editButton;

  // data members:
  int layout;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class BreakpointParameterEditor

/** This class encapsulates the widgets for the per-breakpoint settings that are used in a 
BreakpointModulatorEditor in order to facilitate to treat this part of the editor as one entity (to 
change the color-scheme at once, for example). */

class JUCE_API BreakpointParameterEditor : public AudioModuleEditor, public ChangeBroadcaster, 
  public RSliderListener, public RComboBoxObserver
{

  friend class BreakpointModulatorEditorMulti;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  BreakpointParameterEditor(CriticalSection *newPlugInLock);

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::Modulator object which is to be edited. Make sure to 
  call this function again with a NULL-pointer when the object get deleted for some reason. */
  virtual void setModulatorToEdit(BreakpointModulatorAudioModule* newModulatorToEdit);

  /** Selects a breakpoint and updates the widgets accordingly. */
  virtual void selectBreakpoint(int index);

  /** De-selects the currently selected breakpoint (if any) and updates the GUI accordingly (makes 
  some widgets invisible) */
  virtual void deSelectBreakpoint();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();


protected:

  // pointer to the edited object:
  BreakpointModulatorAudioModule* modulatorToEdit;

  // widgets:
  RTextField *indexLabel, *indexValueLabel;
  RSlider    *timeSlider, *levelSlider, *shapeSlider;
  RButton    *shapeToAllButton;
  RNamedComboBox *shapeComboBox;

  // data members:
  int selectedBreakpointIndex;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class BreakpointModulatorEditor

/** This class is a component, intended to serve as base-class for all components that represent 
some kind of sub-editor, for example an envelope-editor inside an editor for a synthesizer. It is 
also a ChangeBroadcaster such that an outlying editor-class can respond to changes inside the 
sub-editor. Sub-classes must override the purely virtual methods setStateFromXml() and 
getStateAsXml().

\todo: use a TimeGridComboBox in place of the generic one */

class JUCE_API BreakpointModulatorEditor : virtual public AudioModuleEditor, 
  public ChangeBroadcaster, public RComboBoxObserver
{

  friend class BreakpointModulatorEditorCompact;

public:

  enum parameters
  {
    NONE = 0,
    SOME_BREAKPOINT,
    TIME_SCALE,

    //...

    ALL // when we load a preset, everything changes at once
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. 
  \todo: get rid of the pluginLock parameter */
  BreakpointModulatorEditor(CriticalSection *newPlugInLock, 
    BreakpointModulatorAudioModule* newModulatorToEdit);

  /** Destructor. */
  virtual ~BreakpointModulatorEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::Modulator object which is to be edited. Make sure to 
  call this function again with a NULL-pointer when the object get deleted for some reason. */
  virtual void setModulatorToEdit(rosic::BreakpointModulator* newModulatorToEdit);

  /** Sets the juce::Label in which the descriptions for the widgets will appear. */
  //virtual void setDescriptionField(RLabel* newDescriptionField);

  /** De-selects the currently selected breakpoint (if any) and updates the GUI accordingly. */
  virtual void deSelectBreakpoint();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void rComboBoxChanged(RComboBox *rComboBoxThatHasChanged);
  //virtual void rLabelTextChanged(RLabel* rLabelThatHasChanged);
  //virtual void rSliderValueChanged(RSlider *sliderThatHasChanged);
  virtual void paint(Graphics &g);
  virtual void resized();

  /** Updates the sliders, buttons, etc. accordning to the state of the rosic::Modulator object 
  which is being edited. It makes sense to de-select any possibly selected after restoring a state,
  so this might be done here optionally too. */
  virtual void updateWidgetsAccordingToState(bool deSelectBreakpoint);

  /** Calls updateWidgetsAccordingToState(bool) with true as argument (we need to implement this 
  purely virtual function here). */
  virtual void updateWidgetsAccordingToState();

  /** Updates the passed ModulatorCurveEditor "plot" and the grid/snap widgets according to the 
  passed BreakpointModulatorAudioModule "m". Used internally in updateWidgetsAccordingToState. */
  void updatePlotAndGridWidgets(BreakpointModulatorAudioModule* m, ModulatorCurveEditor* plot);

protected:

  /** Assigns the grid/snap widgets to the corresponding Parameters in the passed 
  BreakpointModulatorAudioModule. */
  void assignGridAndSnapWidgets(BreakpointModulatorAudioModule* m);

  /** Returns the grid-interval which belongs to a given interval-index. */
  virtual double gridIntervalFromIndex(int index);

  /** Returns the index which belongs to a given grid-interval. */
  virtual int indexFromGridInterval(double interval);

  /** Returns the time-interval which belongs to a given interval-index. */
  virtual double timeIntervalFromIndex(int index);

  /** Returns the index which belongs to a given time-interval. */
  virtual int indexFromTimeInterval(double interval);

  /** Automatically adjusts the x-axis plot-range according to the current content. */
  virtual void autoAdjustPlotRangeX();

  /** Automatically adjusts the y-axis plot-range according to the current content. */
  virtual void autoAdjustPlotRangeY();

  // pointer to the AudioModule object which is being edited:
  jura::BreakpointModulatorAudioModule* modulatorModule;

  // pointer to the actual Modulator object which is being edited:
  rosic::BreakpointModulator* modulatorToEdit;

  // some rectangles to define functional groups:
  juce::Rectangle<int> breakpointGroupRectangle, timeAndDepthGroupRectangle, snapRectangle;

  // sub-editors:
  BreakpointModulatorGlobalEditor* globalEditor;
  BreakpointParameterEditor*       breakpointParameterEditor;

  // widgets:
  RButton   *snapXButton, *snapYButton;
  RComboBox *gridXComboBox, *gridYComboBox;

  // plot and related stuff:
  ModulatorCurveEditor*      breakpointEditor;
  rsPlotZoomer* breakpointZoomer;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** Compact editor for the BreakpointModulatorAudioModule. 
maybe this needs to ba a ChangeListener, too - same as the non-compact version
*/

class JUCE_API BreakpointModulatorEditorCompact : virtual public AudioModuleEditor
{

public:

  enum layouts
  {
    STANDARD,
    COMPACT,

    NUM_LAYOUTS
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  BreakpointModulatorEditorCompact(CriticalSection *newPlugInLock, 
    BreakpointModulatorAudioModule* breakpointModulatorModuleToEdit);
  virtual ~BreakpointModulatorEditorCompact();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Selects one of the predefined layouts for the widgets. @see layouts */
  virtual void setLayout(int newLayout);

  /** Sets up the bounds of the popup editor relative to the top-left position of the 
  edit-button. */
  virtual void setPopUpEditorBounds(int x, int y, int w, int h);

  /** Sets the headline text for both, the compact editor and the embedded popup editor. */
  virtual void setHeadlineText(const juce::String& newHeadlineText);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();
  virtual void updateWidgetsAccordingToState();
  virtual void updatePlot();


protected:

  // pointers to the edited objects (wrapped and non-wrapped):
  rosic::BreakpointModulator *modulatorToEdit;
  BreakpointModulatorAudioModule *modulatorModuleToEdit;

  // big editor that opens when clicking on the 'More' button:
  BreakpointModulatorEditor *popUpEditor;

  // widgets:
  RButton *editButton;

  // plot and related stuff:
  rsDataPlot  *plot;
  double *xValues, *yValues;
  int    numSamplesInPlot;

  // bounds of the big editor relative to the top-left position of the edit-button:
  int popUpEditorX, popUpEditorY, popUpEditorW, popUpEditorH;

  int layout;

  juce_UseDebuggingNewOperator;
};

#endif  