#ifndef jura_InteractiveCoordinateSystemOld_h
#define jura_InteractiveCoordinateSystemOld_h

/** This class is a rsPlot with an array CoordinatesystemControlHandles which - 
in their simplest form - are just points/dots which can be created, dragged around and removed 
by the mouse. It also provides some functionality which is likely to be useful for subclasses such 
as a snapToGrid()-method. */

class JUCE_API InteractiveCoordinateSystemOld : virtual public rsPlot
{

public:

  enum locatorArrowPositions
  {
    NO_ARROW = 0,
    ARROW_AT_TOP,
    ARROW_AT_BOTTOM,
    ARROW_AT_MIDDLE
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  InteractiveCoordinateSystemOld(const juce::String& name);

  /** Destructor. */
  virtual ~InteractiveCoordinateSystemOld();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Switches the snap-to-grid mode on or off for the coarse x-grid. */
  virtual void setSnapToCoarseGridX(bool shouldSnap) { snapToCoarseGridX = shouldSnap; }

  /** Switches the snap-to-grid mode on or off for the coarse y-grid. */
  virtual void setSnapToCoarseGridY(bool shouldSnap) { snapToCoarseGridY = shouldSnap; }

  /** Switches the snap-to-grid mode on or off for the coarse x-grid. */
  virtual void setSnapToFineGridX(bool shouldSnap)   { snapToFineGridX = shouldSnap; }

  /** Switches the snap-to-grid mode on or off for the coarse x-grid. */
  virtual void setSnapToFineGridY(bool shouldSnap)   { snapToFineGridY = shouldSnap; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns whether of not snapping to the coarse grid for the x-axis is active. */
  virtual bool isSnappingToCoarseGridX() const { return snapToCoarseGridX; }

  /** Returns whether of not snapping to the coarse grid for the y-axis is active. */
  virtual bool isSnappingToCoarseGridY() const { return snapToCoarseGridY; }

  /** Returns whether of not snapping to the fine grid for the x-axis is active. */
  virtual bool isSnappingToFineGridX() const { return snapToFineGridX; }

  /** Returns whether of not snapping to the fine grid for the y-axis is active. */
  virtual bool isSnappingToFineGridY() const { return snapToFineGridY; }

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Overrides mouseDown to open the context menu for the magnetic grids on right-click. */
  virtual void mouseDown (const MouseEvent& e);

  /** Overrides mouseEnter to call the respective function in both base classes. */
  //virtual void mouseEnter(const MouseEvent& e);

  /** Overrides mouseExit to call the respective function in both base classes. */
  //virtual void mouseExit(const MouseEvent& e);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Overides the getStateAsXml()-method from the rsPlot base-class. */
  virtual XmlElement* getStateAsXml(
    const juce::String& stateName = juce::String("InteractiveCoordinateSystemOldState")) const;

  /** Overwrites the setStateFromXml()-method from the rsPlot base-class. */
  virtual bool setStateFromXml(const XmlElement &xmlState);


protected:

  /** Overrides rsPlot::openRightClickPopupMenu() to give more options. */
  void openRightClickPopupMenu();

  /** Quantizes a point to the closest coarseGrid-position (either horizontally, vertically, both
  or none - depending on the snapToCoarseGridX and snapToCoarseGridY flags). */
  virtual void snapToCoarseGrid(double &x, double &y);

  /** Quantizes a point to the closest fineGrid-position (either horizontally, vertically, both
  or none - depending on the snapToFineGridX and snapToFineGridY flags). */
  virtual void snapToFineGrid(double &x, double &y);

  /** Just calls snapToFineGrid() and snapToCoarseGrid() (in that order). */
  virtual void snapToGrid(double &x, double &y);

  /** Draws the control-handles on top of the rsPlot which is should
  be to be drawn before (in the paint()-method). */
  //virtual void drawControlHandles(Graphics &g);

  /** Draws a left locator (for marking a start position inside a sample etc.). */
  virtual void drawLeftLocator(Graphics &g, float x, int arrowPosition = ARROW_AT_TOP,
    const Colour& locatorColour = Colours::darkviolet, juce::Image *targetImage = NULL);

  /** Draws a right locator (for marking an end position inside a sample etc.). */
  virtual void drawRightLocator(Graphics &g, float x, int arrowPosition = ARROW_AT_TOP,
    const Colour& locatorColour = Colours::darkviolet, juce::Image *targetImage = NULL);

  /** Draws a position locator (just a vertical line). */
  virtual void drawCurrentPositionLocator(Graphics &g, float x, int arrowPosition = NO_ARROW,
    const Colour& locatorColour = Colours::red, juce::Image *targetImage = NULL);

  bool   snapToCoarseGridX, snapToCoarseGridY, snapToFineGridX, snapToFineGridY;
  int    mouseX, mouseY;

  double dotRadius;
  Colour dotColour;
  // these two members will probably removed after the test-phase

  // color-scheme management:
  Colour loopLocatorColour, matchedLoopConnectorColour, unmatchedLoopConnectorColour;

  // arrays with the possible snap intervals:
  juce::Array<double> xSnapIntervals;
  juce::Array<double> ySnapIntervals;

  juce_UseDebuggingNewOperator;
};

#endif  
