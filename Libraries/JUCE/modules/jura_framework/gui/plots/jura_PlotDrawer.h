#ifndef jura_PlotDrawer_h
#define jura_PlotDrawer_h

/** This class is used for drawing a plot. Simply create one on the stack and call drawPlot to do 
the actual drawing */

class JUCE_API rsPlotDrawer
{

public:

  /** Constructor. You need to pass references to an rsPlotSettings and an PlotColourSchme object 
  that will determine the appearance of the plot */
  rsPlotDrawer(const rsPlotSettings& plotSettings, const PlotColourScheme& colorScheme);

  //-----------------------------------------------------------------------------------------------
  // \name Drawing with juce::Graphics

  /** Draws the coordinate system onto the given Graphics object within the given rectangle
  (it's not yet tested, if the rectangle bounds actually work).  */
  virtual void drawPlot(Graphics& g, double x, double y, double width, double height);
    // i think, we need members x,y,w,h


  // todo:

  /** Draws the grids, coordinate axes and axis annotations, i.e. everything that appears behind 
  the actual curves or scattered dots or whatever is going to be plotted to represent the data. */
  virtual void drawPlotBackground(Graphics& g);

  /** Draws the caption which should appear in front of the plotted curves. */
  virtual void drawPlotForeground(Graphics& g);

  /** Plots a curve by connecting the given datapoints with lines. */
  virtual void drawWithLines(Graphics& g, int numValues, float* valuesX, float* valuesY);

  /** Draws the givne datapoints as dots like in a scatterplot. */
  virtual void drawAsDots(Graphics& g, int numValues, float* valuesX, float* valuesY);

  //virtual void drawFunction(std::function<double(double)> function, double increment = 1);
  // goes through the x-range with given increment (in pixels) and evaluates the given function
  // at the (unmapped) x-values and draws the function. very convenient but maybe not very 
  // efficient...but that depends on the situation

  // drawWithBars, drawWithStems, etc.



  //-----------------------------------------------------------------------------------------------
  // \name Drawing on an SVG

  /** Analog to the other drawPlot version but draws onto an svg draing instead of a graphics 
  object. Useful for implementing export of plots to svg files. */
  virtual void drawPlot(XmlElement* svg, double x, double y, double w, double h);
   // maybe move into subclass rsPlotDrawerWithSvg...or something
   // maybe have a subclass that uses OpenGL
   // or have an (possibly abstract) rsPlotDrawer baseclass and rsPlotDrawerNative, 
   // rsPlotDrawerOpenGL, rsPlotDrawerSVG subclasses


protected:

  /** Sets up our mapper member according to the settings member and given corrdinates and 
  extents. */
  void setupMapper(double x, double y, double width, double height);

  /** Returns the x-coordinate for the y-axis in model coordinates. */
  double getVerticalAxisX();

  /** Returns the y-coordinate for the x-axis in model coordinates. */
  double getHorizontalAxisY();


  const rsPlotSettings& settings;
  const PlotColourScheme& colors;

  RAPT::rsCoordinateMapper2D<double> mapper;

};

#endif