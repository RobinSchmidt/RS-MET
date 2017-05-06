#ifndef jura_ModulatorCurveEditor_h
#define jura_ModulatorCurveEditor_h

/** This class is intended to be used as a graphical editor for a BreakpointModulator-object.
Breakpoints can be created, removed and dragged around.

\todo: implement new locking strategy and parameter handling */

class JUCE_API ModulatorCurveEditor : virtual public CurveFamilyPlotOld,
  virtual public InteractiveCoordinateSystemOld, public ParameterObserver,
  public ChangeBroadcaster
{

public:

  enum actions
  {
    NO_ACTION = 0,
    BREAKPOINT_INSERTED,
    BREAKPOINT_MODIFIED,
    BREAKPOINT_REMOVED,
    LOOP_START_MODIFIED,
    LOOP_END_MODIFIED
  };

  enum mousableObjects
  {
    NO_OBJECT = 0,
    START_LOCATOR,
    LOOP_START_LOCATOR,
    LOOP_END_LOCATOR,
    SOME_BREAKPOINT
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ModulatorCurveEditor(const juce::String& name = juce::String("ModulatorCurveEditor"));

  /** Destructor. */
  virtual ~ModulatorCurveEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a pointer the the actual rosic::Modulator object which is to be edited. Make sure to 
  call this function again with a NULL-pointer when the object get deleted for some reason. */
  virtual void setModulatorToEdit(rosic::BreakpointModulator* newModulatorToEdit);

  /** Switches loop on or off. */
  //virtual void setLoopMode(bool shouldBeLooped);

  /** Switches sync on or off. */
  //virtual void setSyncMode(bool shouldBeSynced);

  /** Selects a new currently active breakpoint. */
  virtual bool setSelectedBreakpointIndex(int indexToActivate);

  /** Changes the time of the currently selected breakpoint. */
  virtual bool setSelectedBreakpointTime(double newTime, bool broadcastMessage = false);

  /** Changes the level of the currently selected breakpoint. */
  virtual bool setSelectedBreakpointLevel(double newLevel, bool broadcastMessage = false);

  /** Changes the shape of the currently selected breakpoint. */
  virtual bool setSelectedBreakpointShape(int newShape, bool broadcastMessage = false);

  /** Changes the shape all breakpoints at once. */
  virtual void setAllBreakpointShapes(int newShape, bool broadcastMessage = false);

  /** Changes the shape-amount of the currently selected breakpoint. */
  virtual bool setSelectedBreakpointShapeAmount(double newShapeAmount, 
    bool broadcastMessage = false);

  /** Changes the shape all breakpoints at once. */
  virtual void setAllBreakpointShapeAmounts(double newShapeAmount, 
    bool broadcastMessage = false);

  /** Overrides CoordinateSystem::setCurrentRangeX to update the plot. */
  virtual void setCurrentRangeX(double newMinX, double newMaxX);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of controlpoints currently in use. */
  virtual int getNumBreakpoints() const;

  ///** Returns an Rectangle wich encloses the curve. Optionally, a margin can be specified in
  //percent. */
  //virtual CoordinateSystemRange
  //  getMaximumMeaningfulRange(double relativeMarginLeft = 10.0, 
  //  double relativeMarginRight  = 10.0, double relativeMarginTop  = 10.0, 
  //  double relativeMarginBottom = 10.0);

  /** Tells, whether loop is on or off. */
  virtual bool getLoopMode();

  /** Tells, whether sync is on or off. */
  //virtual bool getSyncMode();

  /** Returns the time of a breakpoint (-1.0 if the index is out of range). */
  virtual double getBreakpointTime(int index);

  /** Returns the level of a breakpoint (0.0 if the index is out of range). */
  virtual double getBreakpointLevel(int index);

  /** Returns the shape of a breakpoint (-1 if the index is out of range). */
  virtual int getBreakpointShape(int index);

  /** Returns the shapeAmount of a breakpoint (0.0 if the index is out of range). */
  virtual double getBreakpointShapeAmount(int index);

  /** Returns the index of the currently selected breakpoint (which is the one
  which was most recently  modified and which has the corona). */
  virtual int getSelectedBreakpointIndex();

  /** Returns the time of the currently selected breakpoint. */
  virtual double getSelectedBreakpointTime();

  /** Returns the minimum time-value of the currently selected breakpoint which
  does not violate the constraint of being in between it neighbours. */
  virtual double getSelectedBreakpointMinTime();

  /** Returns the maximum time-value of the currently selected breakpoint which
  does not violate the constraint of being in between it neighbours. */
  virtual double getSelectedBreakpointMaxTime();

  /** Returns the level of the currently selected breakpoint. */
  virtual double getSelectedBreakpointLevel();

  /** Returns the shape of the currently selected breakpoint. */
  virtual int getSelectedBreakpointShape();

  /** Returns the shape-amount of the currently selected breakpoint. */
  virtual double getSelectedBreakpointShapeAmount();

  /** Returns the index of the most recently inserted/modified/removed breakpoint. Which of
  the three possible actions it was, can be retrieved with the getMostRecentAction() function.
  If a breakpoint was removed, the index returned will be the one, which the breakpoint had
  before it was removed. */
  virtual int getMostRecentBreakpointIndex();

  /** Returns the action that was most recently perfomed on a breakpoint (see enumeration
  actions for possible values). The index of the breakpoint on which the action was performed
  can be retrived with the getMostRecentBreakpointIndex() function. */
  virtual int getMostRecentBreakpointAction();

  /** Returns an index of the object which is currently under the mouse-cursor.
  @see mousableObjects */
  virtual int whatIsUnderTheMouseCursor(const MouseEvent &e);

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void parameterChanged(Parameter* parameterThatHasChanged);
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted);
  virtual void mouseDown(const MouseEvent &e);
  virtual void mouseDrag(const MouseEvent &e);
  virtual void mouseMove(const MouseEvent &e);
  virtual void mouseUp  (const MouseEvent &e);
  virtual void resized();

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Fills the plotDataX and plotDataY arrays with new data according to the current state of
  the modulator that is passed. */
  virtual void updatePlotCurveData(int curveIndex, rosic::BreakpointModulator* modulator,
    bool updateGUI);

  /** Calls updatePlotCurveData(int, BreakpointModulator*, bool) with arguments
  editedModulatorIndex, modulatorToEdit, true */
  virtual void updatePlotCurveData();

  /** Updates the maximum range of the inherited CoordinateSystem according to the postions of the
  most extreme breakpoints. */
  virtual void updateMaximumRange(bool alsoUpdateCurrentRange = false);


protected:

  /** Overrides CurveFamilyPlot::plotCurveFamily in order to additionally draw the nodes. */
  virtual void plotCurveFamily(Graphics &g, juce::Image *targetImage = NULL,
    XmlElement *targetSVG = NULL);

  /** Plots the breakpoints as dots of some modulator. */
  virtual void plotBreakpoints(Graphics &g, juce::Image *targetImage, 
    rosic::BreakpointModulator* modulator, const Colour& dotColour);

  /** Plots the loop locators of some modulator. If fullHeight is false, it will draw only
  small locators which are more suitable for unfocused curves when there are more curves in one
  plot. */
  virtual void plotLoopLocators(Graphics &g, juce::Image *targetImage, 
    rosic::BreakpointModulator* modulator, const Colour& locatorColour, bool fullHeight = true);


  /** Pointer to the actual Modulator object which is being edited. */
  rosic::BreakpointModulator* modulatorToEdit;

  int   locatorBeingDragged; // index of the locator which is being dragged 
  // (see enum above)

  // array for the previews of the envelopes:
  //static const int numPreviewSamples = 1000;
  //double previewSamplesX[numPreviewSamples];
  //double previewSamplesY[numPreviewSamples];

  /** Allocates the arrays which hold the data for the plots. */
  virtual void allocatePlotBufferArrays();

  /** Deletes the arrays which hold the data for the plots. */
  virtual void freePlotBufferArrays();

  // buffers for the plot (and related stuff):
  int    numModulators;
  int    numSamplesInPlot;
  //double *plotDataX;
  double *plotDataFlatX; // y-values as flat array
  double **plotDataX;    // 
  double *plotDataFlatY; // y-values as flat array
  double **plotDataY;    // 
  int    editedModulatorIndex; // -1, if none

  int selectedBreakpoint;
  int breakpointBeingDragged;

  int mostRecentBreakpointIndex;
  int mostRecentAction;

  juce_UseDebuggingNewOperator;
};

#endif  
