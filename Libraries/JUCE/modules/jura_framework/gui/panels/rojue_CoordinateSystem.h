#ifndef rojue_CoordinateSystem_h
#define rojue_CoordinateSystem_h

#include "rojue_ThreadedDrawingPanel.h"
#include "../../graphics/rojue_GraphicsTools.h"
#include "../../graphics/rojue_ColourScheme.h"
#include "../../graphics/rojue_GlobalFontInstances.h"
#include "../widgets/rojue_RWidget.h"
#include "../../misc/rojue_FileTools.h"

namespace rojue
{

  /**

  This class is a component, intended to be used as base-class for all components that need some
  underlying coordinate-system, such as function-plots, XY-pads, etc. It takes care of the coordinate
  axes, a coarse and a fine grid, conversion between component-coordinates and the coordinates in the
  desired coordinate-system (which can be lin- or log-scaled).

  */

  class CoordinateSystem : virtual public ThreadedDrawingPanel, virtual public RWidget
  {

    friend class CoordinateSystemZoomer;

  public:

    enum captionPositions
    {
      NO_CAPTION = 0,
      TOP_CENTER,
      CENTER
    };

    enum axisPositions
    {
      INVISIBLE = 0,
      ZERO,
      LEFT,
      RIGHT,
      TOP,
      BOTTOM
    };

    enum axisAnnotationPositions
    {
      NO_ANNOTATION = 0,
      LEFT_TO_AXIS, 
      RIGHT_TO_AXIS,
      ABOVE_AXIS,
      BELOW_AXIS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    CoordinateSystem(const juce::String& newDescription = juce::String(T("some 2D widget")));   

    /** Destructor. */
    virtual ~CoordinateSystem();            

    //---------------------------------------------------------------------------------------------
    // callbacks:

    /** Lets a context menu pop up when the right button is clicked to allow export of the content 
    as image or svg drawing. */
    virtual void mouseDown(const MouseEvent& e);

    /** Overrides mouseEnter for displaying the inspection Label. */
    virtual void mouseEnter(const MouseEvent& e);

    /** Overrides mouseMove for displaying the inspection Label. */
    virtual void mouseMove(const MouseEvent& e);

    /** Overrides mouseExit for displaying the inspection Label. */
    virtual void mouseExit(const MouseEvent& e);

    /**< Overrides mouseDrag to let it be called by a CoordinateSystemZoomer (override it in 
    subclasses, if you want to respond to mouseDrag events). */
    virtual void mouseDrag(const MouseEvent& e);

    /** Overrides mouseUp to let it be called by a CoordinateSystemZoomer (override it in 
    subclasses, if you want to respond to mouseUp events). */
    virtual void mouseUp(const MouseEvent& e);

    /** Overrides mouseDoubleClick to let it be called by a CoordinateSystemZoomer (override it in 
    subclasses, if you want to respond to mouseDoubleClick events). */
    virtual void mouseDoubleClick(const MouseEvent& e);

    /** Overrides mouseWheelMove to let it be called by a CoordinateSystemZoomer (override it in 
    subclasses, if you want to respond to mouseWheelMove events). */
    virtual void mouseWheelMove(const MouseEvent &e, float wheelIncrementX, float wheelIncrementY);



    /** Overrides the resized()-function of the component base-class. */
    virtual void resized();

    /** Overrides the paint-function of the component base-class. */
    virtual void paint(Graphics &g);

    /** Overrides setDirty() to initiate a re-drawing. */
    //virtual void setDirty(bool shouldSetToDirty = true);

    virtual void enablementChanged() { RWidget::enablementChanged(); }


    //-----------------------------------------------------------------------------------------------
    // range-management:


    //-----------------------------------------------------------------------------------------------
    // appearance:

    /** Sets up the colour-scheme from. */
    virtual void setColourScheme(const PlotColourScheme& newColourScheme);

    /** Sets up the colour-scheme from an XmlElement. */
    virtual void setColourSchemeFromXml(const XmlElement& xml);

    /** Changes one of the colours for the graphs if a colour with this index exists (and returns 
    true in this case) - if the index is out of range, it does nothing and returns false. */
    virtual bool changeGraphColour(int index, Colour newColour);

    /** Sets up a caption for the CoordinateSystem and the position where it should appear. */
    virtual void setCaption(const juce::String &newCaption, int newPosition = TOP_CENTER);

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
    virtual void setHorizontalFineGrid(double newGridInterval, bool shouldBeVisible); 

    /** Sets the visibility of the vertical coarse grid. */
    virtual void setVerticalCoarseGridVisible(bool shouldBeVisible); 

    /** Sets the interval of the vertical coarse grid. */
    virtual void setVerticalCoarseGridInterval(double newGridInterval); 

    /** Sets the interval and visibility of the vertical coarse grid. */
    virtual void setVerticalCoarseGrid(double newGridInterval, bool shouldBeVisible); 

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

    /** Sets the unit of the angle (as used by the angular grid) to degrees. If
    false, radiant will be assumed. */
    virtual void setAngleUnitToDegrees(bool shouldBeInDegrees = true);

    /** Informs, if the horizontal coarse grid is visible. */
    virtual bool isHorizontalCoarseGridVisible();

    /** Informs, if the horizontal fine grid is visible. */
    virtual bool isHorizontalFineGridVisible();

    /** Informs, if the vertical coarse grid is visible. */
    virtual bool isVerticalCoarseGridVisible();

    /** Informs, if the vertical fine grid is visible. */
    virtual bool isVerticalFineGridVisible();

    /** Informs, if the radial coarse grid is visible. */
    virtual bool isRadialCoarseGridVisible();

    /** Informs, if the radial fine grid is visible. */
    virtual bool isRadialFineGridVisible();

    /** Informs, if the angular coarse grid is visible. */
    virtual bool isAngularCoarseGridVisible();

    /** Informs, if the angular fine grid is visible. */
    virtual bool isAngularFineGridVisible();

    /** Returns the interval of the horizontal coarse grid. */
    virtual double getHorizontalCoarseGridInterval();

    /** Returns the interval of the horizontal fine grid. */
    virtual double getHorizontalFineGridInterval();

    /** Returns the interval of the vertical coarse grid. */
    virtual double getVerticalCoarseGridInterval();

    /** Returns the interval of the vertical fine grid. */
    virtual double getVerticalFineGridInterval();

    /** Returns the interval of the radial coarse grid. */
    virtual double getRadialCoarseGridInterval();

    /** Returns the interval of the radial fine grid. */
    virtual double getRadialFineGridInterval();

    /** Returns the interval of the angular coarse grid. */
    virtual double getAngularCoarseGridInterval();

    /** Returns the interval of the angular fine grid. */
    virtual double getAngularFineGridInterval();

    /** Decides if either the x-axis or the y-axis or both should be logarithmically scaled and 
    sets up the base for the logarithms. */
    virtual void useLogarithmicScale(bool shouldBeLogScaledX, bool shouldBeLogScaledY, 
      double newLogBaseX = 2.0, double newLogBaseY = 2.0);





    /** Decides, if the x-axis should be logarithmically scaled and sets up the base for the 
    logarithm. */
    virtual void useLogarithmicScaleX(bool shouldBeLogScaledX, double newLogBaseX = 2.0);

    /** Informs, whether the x-axis is logarithmically scaled or not. */
    virtual bool isLogScaledX();

    /** Decides, if the y-axis should be logarithmically scaled and sets up the base for the 
    logarithm. */
    virtual void useLogarithmicScaleY(bool shouldBeLogScaledY, double newLogBaseY = 2.0);

    /** Informs, whether the y-axis is logarithmically scaled or not. */
    virtual bool isLogScaledY();

    /** Sets the labels for the axes and their position. */
    virtual void setAxisLabels(const juce::String &newLabelX, const juce::String &newLabelY, 
      int newLabelPositionX = ABOVE_AXIS, int newLabelPositionY = RIGHT_TO_AXIS);

    /** Sets the label for the x-axis and its position. */
    virtual void setAxisLabelX(const juce::String &newLabelX, 
      int newLabelPositionX = ABOVE_AXIS);

    /** Sets the label for the y-axis and its position. */
    virtual void setAxisLabelY(const juce::String &newLabelY, 
      int newLabelPositionY = RIGHT_TO_AXIS);

    /** Switches x-value annotation between below or above the x-axis (or off). */
    virtual void setAxisValuesPositionX(int newValuesPositionX);

    /** Switches y-value annotation between left to or right to the y-axis (or off). */
    virtual void setAxisValuesPositionY(int newValuesPositionY);

    /** This function is used to pass a function-pointer. This pointer has to be the address of a 
    function which has a double-parameter and a juce::String as return-value. The function will be 
    used to convert the values on the x-axis into corresponding strings for display. */
    virtual void setStringConversionFunctionX(juce::String (*newConversionFunction) 
      (double valueToBeConverted) );

    /** see setStringConversionFunctionX() - same for y-axis.  */
    virtual void setStringConversionFunctionY(juce::String (*newConversionFunction) 
      (double valueToBeConverted) );

    /** Selects whether or not the field with the values should pop up. */
    virtual void setValueFieldPopup(bool shouldPopUp);

    //---------------------------------------------------------------------------------------------
    // state-management:

    /** Creates an XmlElement from the current state and returns it. */
    virtual XmlElement* getStateAsXml(const juce::String& stateName = 
      juce::String(T("CoordinateSystemState"))) const;

    /** Restores a state based on an XmlElement which should have been created with the 
    getStateAsXml()-function. */
    virtual bool setStateFromXml(const XmlElement &xmlState);

    /** Returns the drawing as SVG compliant XmlElement. The caller must take care to delete the 
    pointer to the XmlElement when it's not needed anymore. */
    virtual XmlElement* getPlotAsSVG(int width, int height);

    /** Renders the plot to an image object of given width and height. The caller must take care to 
    delete the pointer to the image when it's not needed anymore. */
    virtual Image* getPlotAsImage(int width, int height);

    /** Opens a dialog window to export the content of the CoordinateSystem to a png-image file or 
    svg vector drawing file. */
    virtual void openExportDialog(int defaultWidth, int defaultHeight, 
      const juce::String &defaultFormat, const File& defaultTargetFile);

    /** With this function, the automatic re-rendering of the underlying image can be turned on or 
    off. If on (default), everytime you change a parameter which will change the appearance of 
    the CoordinateSystem, it will be re-rendered. However, when you want to change many 
    parameters at a time, this can be wasteful in terms of CPU-load. In these cases, it can be
    useful to switch the automatic re-rendering temporarily off. */
    //virtual void setAutoReRendering(bool shouldAutomaticallyReRender);

    /** We define a member for the mouse-cursor which is to be showed in order to let a 
    CoordinateSystemZoomer access this. This is required, because the zoomer object must be above
    the actual CoordinateSystem and therefore prevent the CoordinateSystem to set it's own 
    MouseCursor. Instead we just assign the member mouse-cursor and let the zoomer retrieve it. */
    MouseCursor currentMouseCursor;

    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Overrides drawComponent inherited from ThreadedDrawingComponent in order to do the actual 
    drawing operations. */
    virtual void drawComponent(Image* imageToDrawOnto);

    /** Opens the PopupMenu that appears on right clicks. */
    void openRightClickPopupMenu();

    /** Draws all the stuff either on the internal image which will be displayed as the components 
    content or on an arbitrary image (if only the first optional pointer argumnet is nonzero) or on
    an arbitrary image and an SVG compliant XmlElement (if both poiters are nonzero). */
    virtual void drawCoordinateSystem(Graphics &g, Image* targetImage = NULL, 
      XmlElement* targetSVG = NULL);

    /** Draws a horizontal grid with a given interval in a given colour. Gets called by 
    drawCoordinateSystem(). */
    virtual void drawHorizontalGrid(Graphics &g, double interval, 
      bool exponentialSpacing, Colour gridColour, float lineThickness, 
      Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws a vertical grid with a given interval in a given colour. Gets called by 
    drawCoordinateSystem(). */
    virtual void drawVerticalGrid(Graphics &g, double interval, 
      bool exponentialSpacing, Colour gridColour, float lineThickness, 
      Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws a radial grid with a given interval in a given colour. Gets called by 
    drawCoordinateSystem(). */
    virtual void drawRadialGrid(Graphics &g, double interval, bool exponentialSpacing, 
      Colour gridColour, float lineThickness, Image* targetImage = NULL, 
      XmlElement* targetSVG = NULL);

    /** Draws an angular grid with a given interval in a given colour. Gets called by 
    drawCoordinateSystem(). */
    virtual void drawAngularGrid(Graphics &g, double interval, Colour gridColour, 
      float lineThickness, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the caption/headline. Gets called by drawCoordinateSystem(). */
    virtual void drawCaption(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the x-axis. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisX(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the y-axis. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisY(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the x-axis' label. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisLabelX(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the y-axis' label. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisLabelY(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the numeric values at the x-axis. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisValuesX(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Draws the numeric values at the y-axis. Gets called by drawCoordinateSystem(). */
    virtual void drawAxisValuesY(Graphics &g, Image* targetImage = NULL, XmlElement* targetSVG = NULL);

    /** Function for converting the x- and y-coordinate values into the corresponding coordinates 
    in an image. */
    virtual void transformToImageCoordinates(double &x, double &y, const Image* theImage) const;

    /** Function for converting the x- and y-coordinate values measured in the image's coordinate
    system to the corresponding coordinates of our plot. */
    virtual void transformFromImageCoordinates(double &x, double &y, const Image* theImage) const;

    /** Function for converting the x- and y-coordinate values into the corresponding coordinates in
    the component (double precision version).*/
    virtual void transformToComponentsCoordinates(double &x, double &y) const;

    /** Function for converting the x- and y-coordinate values into the corresponding coordinates in 
    the component (single precision version).*/
    virtual void transformToComponentsCoordinates(float &x, float &y) const;

    /** Function for converting the x- and y-coordinate values measured in the components coordinate
    system to the corresponding coordinates of our plot (double precision version). */
    virtual void transformFromComponentsCoordinates(double &x, double &y) const;

    /** Function for converting the x- and y-coordinate values measured in the components coordinate
    system to the corresponding coordinates of our plot (single precision version). */
    virtual void transformFromComponentsCoordinates(float &x, float &y) const;

    /** Adds a line to an SVG drawing. */
    virtual void addLineToSvgDrawing(XmlElement* theSVG, float x1, float y1, float x2, float y2, 
      float thickness, Colour colour, bool withArrowHead = false);

    /** Adds text-string to an SVG drawing. */
    virtual void addTextToSvgDrawing(XmlElement* theSVG, juce::String theText, float x, float y, 
      Justification justification = Justification::centredLeft );

    /** Draws a text on the Coordinatesystem using a BitmapFont. */
    virtual void drawBitmapText(Graphics &g, const juce::String &text, double x, double y, double w, 
      double h, BitmapFont const* font, Justification justification); 

    /** Updates the image object (re-draws it). Will be called, when something about the 
    CoordinateSystem's appearance-settings was changed. */
    //virtual void updateBackgroundImage();

    /** Updates the scale-factors which are needed when transforming from the CoordinateSystem's 
    coordinates to Component's coordinates and vice versa. Will be called by setBounds(), 
    setRange() and useLogarithmicScale(). */
    virtual void updateScaleFactors();

    /** Returns either the height of this component or the height of the image (if the pointer is 
    non-NULL). */
    double getPlotHeight(Image *targetImage = NULL);

    /** Returns either the height of this component or the height of the image (if the pointer is 
    non-NULL). */
    double getPlotWidth(Image *targetImage = NULL);

    //double scaleX;
    //double scaleY;
    double pixelsPerIntervalX;
    double pixelsPerIntervalY;

    int    axisPositionX;
    int    axisPositionY;
    int    axisLabelPositionX;
    int    axisLabelPositionY;
    int    axisValuesPositionX;
    int    axisValuesPositionY;

    juce::String axisLabelX;
    juce::String axisLabelY;

    bool   horizontalCoarseGridIsVisible;
    bool   horizontalFineGridIsVisible;
    bool   verticalCoarseGridIsVisible;
    bool   verticalFineGridIsVisible;
    bool   radialCoarseGridIsVisible;
    bool   radialFineGridIsVisible;
    bool   angularCoarseGridIsVisible;
    bool   angularFineGridIsVisible;

    //double horizontalCoarseGridInterval;
    //double horizontalFineGridInterval;
    //double verticalCoarseGridInterval;
    //double verticalFineGridInterval;
    double radialCoarseGridInterval;
    double radialFineGridInterval;
    double angularCoarseGridInterval;
    double angularFineGridInterval;

    bool   angleIsInDegrees;

    bool   logScaledX;
    double logBaseX;
    bool   logScaledY;
    double logBaseY;
    bool   logScaledRadius;
    double logBaseRadius;


    bool   valuePopup;
    bool   useBitmapFont;

    int          captionPosition;
    juce::String captionString;
    DrawableText caption;

    //Image*  backgroundImage; 
    /**< This image will be used for the appearance of the coodinate system, it will be updated via
    the updateBackgroundImage()-function when something about the coordinate-system (axes, grid, 
    etc.) changes - but only if autoReRender is true (which it is by default). */

    //bool autoReRenderImage;
    // see above

    TextEditor* inspectionField;
    /**< This is a text-field which shows up with the current coordinates when the user moves the 
    mouse-coursor over the coordinate system. */

    juce::String (*stringConversionFunctionX) (double valueToConvert);
    /**< A pointer to the function which converts a x-value into a juce-string. */

    juce::String (*stringConversionFunctionY) (double valueToConvert);
    /**< A pointer to the function which converts a y-value into a juce-string. */

    // color-scheme management:
    PlotColourScheme colourScheme;
    /*
    Colour topLeftColour, topRightColour, bottomLeftColour, bottomRightColour, 
      axesColour, coarseGridColour, fineGridColour;
    juce::Array<Colour> plotColours;
    */


    //---------------------------------------------------------------------------------------------
    // make obsolete inherited methods unavailable to client code:

    virtual void setColours(const Colour newBackgroundColour, const Colour newOutlineColour, 
      const Colour newHandleColour, const Colour newTextColour, const Colour newSpecialColour1,
      const Colour newSpecialColour2) {}; 
    virtual void setColourScheme(const WidgetColourScheme& newColourScheme) {}

  };

}

#endif  // rojue_CoordinateSystem_h
