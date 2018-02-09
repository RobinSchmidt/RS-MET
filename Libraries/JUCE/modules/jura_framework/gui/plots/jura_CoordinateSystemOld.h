#ifndef jura_CoordinateSystemOld_h
#define jura_CoordinateSystemOld_h 

/** This class is a component, intended to be used as base-class for all components that need some
underlying coordinate-system, such as function-plots, XY-pads, etc. It takes care of the coordinate
axes, a coarse and a fine grid, conversion between component-coordinates and the coordinates in the
desired coordinate-system (which can be lin- or log-scaled).

\todo:

-rename to rsPlot, rsPlotBase or rsPlotComponent  */

class JUCE_API CoordinateSystemOld : virtual public DescribedComponent
{

  friend class CoordinateSystemZoomer;

public:



  CoordinateSystemOld(const juce::String& newDescription = juce::String("some 2D widget"));
  virtual ~CoordinateSystemOld();

  //-----------------------------------------------------------------------------------------------
  // component-overrides:

  /** Lets a context menu pop up when the right button is clicked to allow export of the content
  as image or svg drawing. */
  virtual void mouseDown(const MouseEvent& e);

  /** Overrides mouseEnter for displaying the inspection Label. */
  virtual void mouseEnter(const MouseEvent& e);

  /** Overrides mouseMove for displaying the inspection Label. */
  virtual void mouseMove(const MouseEvent& e);


  virtual void resized();
  virtual void paint(Graphics &g);

  //-----------------------------------------------------------------------------------------------
  // range-management:

  virtual void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  /**< Sets the maximum for the currently visible range. For logarithmic x- and/or y-axis-scaling,
  make sure that the respective minimum value is greater than zero! */

  virtual void setMaximumRange(rsPlotRange newMaximumRange);
  /**< Sets the maximum for the currently visible range. */

  virtual void setMaximumRangeX(double newMinX, double newMaxX);
  /**< Sets the maximum visible range for the y-axis. */

  virtual void setMaximumRangeY(double newMinY, double newMaxY);
  /**< Sets the maximum visible range for the y-axis. */

  virtual rsPlotRange getMaximumRange() { return plotSettings.maximumRange; }
  /**< Returns the maximum for the currently visible range. */

  virtual void setMaximumRangeMinX(double newMinX);
  /**< Sets the minimum value for the range of x. */

  virtual double getMaximumRangeMinX() { return plotSettings.maximumRange.getMinX(); }
  /**< Returns the minimum value for the range of x. */

  virtual void setMaximumRangeMaxX(double newMaxX);
  /**< Sets the maximum value for the range of x. */

  virtual double getMaximumRangeMaxX() { return plotSettings.maximumRange.getMaxX(); }
  /**< Returns the maximum value for the range of x. */

  virtual void setMaximumRangeMinY(double newMinY);
  /**< Sets the minimum value for the range of y. */

  virtual double getMaximumRangeMinY() { return plotSettings.maximumRange.getMinY(); }
  /**< Returns the minimum value for the range of y. */

  virtual void setMaximumRangeMaxY(double newMaxY);
  /**< Sets the maximum value for the range of y. */

  virtual double getMaximumRangeMaxY() { return plotSettings.maximumRange.getMaxY(); }
  /**< Returns the maximum value for the range of y. */

  virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  /**< Sets the currently visible range. For logarithmic x- and/or y-axis-scaling, make sure that
  the respective minimum value is greater than zero! */

  virtual void setCurrentRange(rsPlotRange newRange);
  /**< Sets the currently visible range. */

  virtual void setCurrentRangeX(double newMinX, double newMaxX);
  /**< Sets the currently visible range for the y-axis. */

  virtual void setCurrentRangeY(double newMinY, double newMaxY);
  /**< Sets the currently visible range for the y-axis. */

  virtual rsPlotRange getCurrentRange() { return plotSettings.currentRange; }
  /**< Returns the currently visible range. */

  virtual void setCurrentRangeMinX(double newMinX);
  /**< Sets the minimum value of x. */

  virtual double getCurrentRangeMinX() { return plotSettings.currentRange.getMinX(); }
  /**< Returns the minimum value of x. */

  virtual void setCurrentRangeMaxX(double newMaxX);
  /**< Sets the maximum value of x. */

  virtual double getCurrentRangeMaxX() { return plotSettings.currentRange.getMaxX(); }
  /**< Returns the maximum value of x. */

  virtual void setCurrentRangeMinY(double newMinY);
  /**< Sets the minimum value of y. */

  virtual double getCurrentRangeMinY() { return plotSettings.currentRange.getMinY(); }
  /**< Returns the minimum value of y. */

  virtual void setCurrentRangeMaxY(double newMaxY);
  /**< Sets the maximum value of y. */

  virtual double getCurrentRangeMaxY() { return plotSettings.currentRange.getMaxY(); }
  /**< Returns the maximum value of y. */

  /** Returns a string that represents the info-line to be shown when mouse is over pixel (x,y). */
  virtual juce::String getInfoLineForPixelPosition(int x, int y);


  virtual PlotColourScheme getPlotColourScheme() const { return plotColourScheme; }

  //-----------------------------------------------------------------------------------------------
  // \name Appearance Setup

  /** Sets up the colour-scheme. */
  virtual void setColourScheme(const PlotColourScheme& newColourScheme);

  /** Sets up the colour-scheme to be used for the right-click popup menu. */
  virtual void setPopUpColourScheme(const WidgetColourScheme& newColourScheme)
  {
    popUpColourScheme = newColourScheme;
  }

  /** Sets up the colour-scheme from an XmlElement. */
  virtual void setColourSchemeFromXml(const XmlElement& xml);

  /** Sets up a caption for the CoordinateSystemOld and the position where it should appear. */
  virtual void setCaption(const juce::String &newCaption, 
    int newPosition = rsPlotSettings::TOP_CENTER);

  /** Sets the position of the x-axis. For possible values see enum positions. */
  virtual void setAxisPositionX(int newAxisPositionX);

  /** Sets the position of the y-axis. For possible values see enum positions. */
  virtual void setAxisPositionY(int newAxisPositionY);

  /** Sets up several x-axis parameters at once */
  virtual void setupAxisX(double newMin, double newMax, bool shouldBeLogScaled, double newLogBase,
    int newAxisPosition, double newCoarseGridInterval, double newFineGridInterval);

  /** Sets up several y-axis parameters at once */
  virtual void setupAxisY(double newMin, double newMax, bool shouldBeLogScaled, double newLogBase,
    int newAxisPosition, double newCoarseGridInterval, double newFineGridInterval);

  /** Sets the visibility of the horizontal coarse grid. */
  virtual void setHorizontalCoarseGridVisible(bool shouldBeVisible);

  /** Sets the interval of the horizontal coarse grid. */
  virtual void setHorizontalCoarseGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the horizontal coarse grid. */
  virtual void setHorizontalCoarseGrid(double newGridInterval, bool shouldBeVisible);

  /** Sets the visibility of the horizontal fine grid. */
  virtual void setHorizontalFineGridVisible(bool shouldBeVisible);

  /** Sets the interval of the horizontal fine grid. */
  virtual void setHorizontalFineGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the horizontal fine grid. */
  virtual void setHorizontalFineGrid(double newGridInterval, bool   shouldBeVisible);

  /** Sets the visibility of the vertical coarse grid. */
  virtual void setVerticalCoarseGridVisible(bool shouldBeVisible);

  /** Sets the interval of the vertical coarse grid. */
  virtual void setVerticalCoarseGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the vertical coarse grid. */
  virtual void setVerticalCoarseGrid(double newGridInterval,
    bool   shouldBeVisible);

  /** Sets the visibility of the vertical fine grid. */
  virtual void setVerticalFineGridVisible(bool shouldBeVisible);

  /** Sets the interval of the vertical fine grid. */
  virtual void setVerticalFineGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the vertical fine grid. */
  virtual void setVerticalFineGrid(double newGridInterval, bool shouldBeVisible);

  /** Sets the visibility of the radial coarse grid. */
  virtual void setRadialCoarseGridVisible(bool shouldBeVisible);

  /** Sets the interval of the radial coarse grid. */
  virtual void setRadialCoarseGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the radial coarse grid. */
  virtual void setRadialCoarseGrid(double newGridInterval, bool shouldBeVisible);

  /** Sets the visibility of the radial fine grid. */
  virtual void setRadialFineGridVisible(bool shouldBeVisible);

  /** Sets the interval of the radial fine grid. */
  virtual void setRadialFineGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the radial fine grid. */
  virtual void setRadialFineGrid(double newGridInterval, bool shouldBeVisible);

  /** Sets the visibility of the angular coarse grid. */
  virtual void setAngularCoarseGridVisible(bool shouldBeVisible);

  /** Sets the interval of the angular coarse grid. */
  virtual void setAngularCoarseGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the angular coarse grid. */
  virtual void setAngularCoarseGrid(double newGridInterval, bool shouldBeVisible);

  /** Sets the visibility of the angular fine grid. */
  virtual void setAngularFineGridVisible(bool shouldBeVisible);

  /** Sets the interval of the angular fine grid. */
  virtual void setAngularFineGridInterval(double newGridInterval);

  /** Sets the interval and visibility of the angular fine grid. */
  virtual void setAngularFineGrid(double newGridInterval, bool shouldBeVisible);

  //-----------------------------------------------------------------------------------------------
  // \name Appearance Inquiry

  virtual bool isHorizontalCoarseGridVisible();
  /**< Informs, if the horizontal coarse grid is visible. */

  virtual bool isHorizontalFineGridVisible();
  /**< Informs, if the horizontal fine grid is visible. */

  virtual bool isVerticalCoarseGridVisible();
  /**< Informs, if the vertical coarse grid is visible. */

  virtual bool isVerticalFineGridVisible();
  /**< Informs, if the vertical fine grid is visible. */

  virtual bool isRadialCoarseGridVisible();
  /**< Informs, if the radial coarse grid is visible. */

  virtual bool isRadialFineGridVisible();
  /**< Informs, if the radial fine grid is visible. */

  virtual bool isAngularCoarseGridVisible();
  /**< Informs, if the angular coarse grid is visible. */

  virtual bool isAngularFineGridVisible();
  /**< Informs, if the angular fine grid is visible. */

  virtual double getHorizontalCoarseGridInterval();
  /**< Returns the interval of the horizontal coarse grid. */

  virtual double getHorizontalFineGridInterval();
  /**< Returns the interval of the horizontal fine grid. */

  virtual double getVerticalCoarseGridInterval();
  /**< Returns the interval of the vertical coarse grid. */

  virtual double getVerticalFineGridInterval();
  /**< Returns the interval of the vertical fine grid. */

  virtual double getRadialCoarseGridInterval();
  /**< Returns the interval of the radial coarse grid. */

  virtual double getRadialFineGridInterval();
  /**< Returns the interval of the radial fine grid. */

  virtual double getAngularCoarseGridInterval();
  /**< Returns the interval of the angular coarse grid. */

  virtual double getAngularFineGridInterval();
  /**< Returns the interval of the angular fine grid. */

  virtual void useLogarithmicScale(bool   shouldBeLogScaledX,
    bool   shouldBeLogScaledY,
    double newLogBaseX = 2.0,
    double newLogBaseY = 2.0);
  /**< Decides if either the x-axis or the y-axis or both should be
  logarithmically scaled and sets up the base for the logarithms. */

  virtual void useLogarithmicScaleX(bool   shouldBeLogScaledX,
    double newLogBaseX = 2.0);
  /**< Decides, if the x-axis should be logarithmically scaled and sets up the
  base for the logarithm. */

  virtual bool isLogScaledX();
  /**< Informs, whether the x-axis is logarithmically scaled or not. */

  virtual void useLogarithmicScaleY(bool   shouldBeLogScaledY,
    double newLogBaseY = 2.0);
  /**< Decides, if the y-axis should be logarithmically scaled and sets up the
  base for the logarithm. */

  virtual bool isLogScaledY();
  /**< Informs, whether the y-axis is logarithmically scaled or not. */

  virtual void setAxisLabels(const juce::String &newLabelX, const juce::String &newLabelY,
    int newLabelPositionX = rsPlotSettings::ABOVE_AXIS, 
    int newLabelPositionY = rsPlotSettings::RIGHT_TO_AXIS);
  /**< Sets the labels for the axes and their position. */

  virtual void setAxisLabelX(const juce::String &newLabelX, 
    int newLabelPositionX = rsPlotSettings::ABOVE_AXIS);
  /**< Sets the label for the x-axis and its position. */

  virtual void setAxisLabelY(const juce::String &newLabelY, 
    int newLabelPositionY = rsPlotSettings::RIGHT_TO_AXIS);
  /**< Sets the label for the y-axis and its position. */

  virtual void setAxisValuesPositionX(int newValuesPositionX);
  /**< Switches x-value annotation between below or above the x-axis
  (or off). */

  virtual void setAxisValuesPositionY(int newValuesPositionY);
  /**< Switches y-value annotation between left to or right to the y-axis
  (or off). */

  /** This function is used to pass a function-pointer with the address of a function which has a 
  double-parameter and a juce::String as return-value. The function will be used to convert the 
  values on the x-axis into corresponding strings for display on the axis. */
  virtual void setStringConversionForAxisX(
    juce::String (*newConversionFunction) (double valueToBeConverted));

  /** Like setStringConversionForAxisX but for y-axis. */
  virtual void setStringConversionForAxisY(
    juce::String (*newConversionFunction) (double valueToBeConverted));

  /** Like setStringConversionForAxisX but not for the axis but for the info-line. */
  virtual void setStringConversionForInfoLineX(
    juce::String (*newConversionFunction) (double valueToBeConverted));

  /** Like setStringConversionForAxisY but not for the axis but for the info-line. */
  virtual void setStringConversionForInfoLineY(
    juce::String (*newConversionFunction) (double valueToBeConverted));

  //-----------------------------------------------------------------------------------------------
  // new getters



  //-----------------------------------------------------------------------------------------------
  // state-management:


  /** Creates an XmlElement from the current state and returns it. */
  virtual XmlElement* getStateAsXml(
    const juce::String& stateName = juce::String("CoordinateSystemState")) const;

  /** Restores a state based on an XmlElement which should have been created
  with the getStateAsXml()-function. */
  virtual bool setStateFromXml(const XmlElement &xmlState);





  virtual XmlElement* getPlotAsSVG(int width, int height);
  /**< Returns the drawing as SVG compliant XmlElement. The caller must take care to delete the
  pointer to the XmlElement when it's not needed anymore. */

  virtual Image* getPlotAsImage(int width, int height);
  /**< Renders the plot to an image object of given width and height. The caller must take care to
  delete the pointer to the image when it's not needed anymore. */

  virtual void openExportDialog(int defaultWidth, int defaultHeight, 
    const juce::String &defaultFormat, const juce::File& defaultTargetFile);
  /**< Opens a dialog window to export the content of the CoordinateSystemOld to a png-image file 
  or svg vector drawing file. */

  virtual void setAutoReRendering(bool shouldAutomaticallyReRender);
  /**< With this function, the automatic re-rendering of the underlying image can be turned on or
  off. If on (default), everytime you change a parameter which will change the appearance of
  the CoordinateSystemOld, it will be re-rendered. However, when you want to change many
  parameters at a time, this can be wasteful in terms of CPU-load. In these cases, it can be
  useful to switch the automatic re-rendering temporarily off. */

  MouseCursor currentMouseCursor;
  /**< We define a member for the mouse-cursor which is to be showed in order to let a
  CoordinateSystemZoomer access this. This is required, because the zoomer object must be above
  the actual CoordinateSystemOld and therefore prevent the CoordinateSystemOld to set it's own
  MouseCursor. Instead we just assign the member mouse-cursor and let the zoomer retrieve it. */

protected:

  /** Opens the PopupMenu that appears on right clicks. */
  void openRightClickPopupMenu();


  virtual void drawCoordinateSystem(Graphics &g);

  //virtual void drawCoordinateSystem(Graphics &g, Image* targetImage = NULL, 
  //  XmlElement* targetSVG = NULL);
  /**< Draws all the stuff either on the internal image which will be displayed as the components
  content or on an arbitrary image (if only the first optional pointer argumnet is nonzero) or on
  an arbitrary image and an SVG compliant XmlElement (if both poiters are nonzero). */


  //virtual void drawCaption(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);
  /**< Draws the caption/headline. Gets called by drawCoordinateSystem(). */



  // also to be obsolete soon:

  virtual void transformToImageCoordinates(double &x, double &y, const Image* theImage);
  /**< Function for converting the x- and y-coordinate values into the corresponding coordinates
  in an image. */

  virtual void transformFromImageCoordinates(double &x, double &y, const Image* theImage);
  /**< Function for converting the x- and y-coordinate values measured in the image's coordinate
  system to the corresponding coordinates of our plot. */

  virtual void transformToComponentsCoordinates(double &x, double &y);
  /**< Function for converting the x- and y-coordinate values into the corresponding coordinates in
  the component (double precision version).*/

  virtual void transformToComponentsCoordinates(float &x, float &y);
  /**< Function for converting the x- and y-coordinate values into the corresponding coordinates in
  the component (single precision version).*/

  virtual void transformFromComponentsCoordinates(double &x, double &y);
  /**< Function for converting the x- and y-coordinate values measured in the components coordinate
  system to the corresponding coordinates of our plot (double precision version). */

  virtual void transformFromComponentsCoordinates(float &x, float &y);
  /**< Function for converting the x- and y-coordinate values measured in the components coordinate
  system to the corresponding coordinates of our plot (single precision version). */



  /** Updates the image object (re-draws it). Will be called, when something about the
  CoordinateSystemOld's appearance-settings was changed. */
  virtual void updateBackgroundImage();


  /** Sets up the output range (i.e. the pixel width and height) in our coordinateMapper. If a 
  non-nullptr is passed for targetImage, the image size will be used, else if a non-nullptr for
  the targetSVG is passed, its size will be used (the xml should already have "width" and "height"
  attributes), else this Component's size will be used. */
  void updateMapperOutputRange(Image* targetImage = nullptr, XmlElement* targetSVG = nullptr);

  /** Updates the input range of our coordinate mapper. */
  void updateMapperInputRange();



  double getPlotHeight(Image *targetImage = NULL);
  /**< Returns either the height of this component or the height of the image (if the pointer is
  non-NULL). */

  double getPlotWidth(Image *targetImage = NULL);
  /**< Returns either the height of this component or the height of the image (if the pointer is
  non-NULL). */

  rsPlotSettings plotSettings;


  // new - to be used soon in the drawing code:
  RAPT::rsCoordinateMapper2D<double> coordinateMapper; 
    // get rid of that - it's handled in rsPlotDrawer now...o...but we need to map form 
    // pixel-coordinates to component coordinates for node-editing..but maybe that can be handled
    // in rsNodeEditor? i think, it has a mapper itself?


  bool showPositionAsDescription;
  bool showPopUpOnRightClick;


  Image*  backgroundImage;
  /**< This image will be used for the appearance of the coodinate system, it will be updated via
  the updateBackgroundImage()-function when something about the coordinate-system (axes, grid,
  etc.) changes - but only if autoReRender is true (which it is by default). */

  bool autoReRenderImage;
  // see above


  juce::String (*stringConversionForInfoLineX) (double valueToConvert);
  juce::String (*stringConversionForInfoLineY) (double valueToConvert);

  // color-scheme management:
  PlotColourScheme   plotColourScheme;
  WidgetColourScheme popUpColourScheme;

  //---------------------------------------------------------------------------------------------
  // make obsolete inherited methods unavailable to client code:

  /*
  virtual void setColours(const Colour newBackgroundColour, const Colour newOutlineColour,
    const Colour newHandleColour, const Colour newTextColour, const Colour newSpecialColour1,
    const Colour newSpecialColour2) {};
    */
  //virtual void setColourScheme(const WidgetColourScheme& newColourScheme) {};
  //virtual WidgetColourScheme getColourScheme() const { return RWidget::getColourScheme(); }

  juce_UseDebuggingNewOperator;
};


#endif
