#ifndef jura_PlotDrawer_h
#define jura_PlotDrawer_h

/** A class to simplify line-drawing (not yet tested). */

class JUCE_API rsPolyLineDrawer
{

public:

  rsPolyLineDrawer(Graphics& g, const RAPT::rsCoordinateMapper2D<double>& coordinateMapper, 
    float lineThickness) : gfx(g), mapper(coordinateMapper), thickness(lineThickness)
  {}

  /** Starts a new polyline at given x,y (in model coordinates). */
  void startPolyLine(double x, double y)
  {
    xOld = (float)mapper.mapX(x);
    yOld = (float)mapper.mapY(y);
  }

  /** Draws a line from the previous line endpoint (or polyline start point) to the given x,y (in
  model coordinates). */
  void lineTo(double x, double y)
  {
    x = mapper.mapX(x);
    y = mapper.mapY(y);
    gfx.drawLine((float)xOld, (float)yOld, (float)x, (float)y, thickness);
    xOld = (float)x;
    yOld = (float)y;
  }

protected:

  Graphics& gfx;
  const RAPT::rsCoordinateMapper2D<double>& mapper;
  float xOld, yOld;
  float thickness;

};

//=================================================================================================

/** This class is used for drawing a plot. Simply create one on the stack and call drawPlot to do 
the actual drawing */

class JUCE_API rsPlotDrawer
{

public:

  /** Constructor. You need to pass references to an rsPlotSettings and an PlotColourSchme object 
  that will determine the appearance of the plot */
  rsPlotDrawer(const rsPlotSettings& plotSettings, const PlotColourScheme& colorScheme,
    double xLeft, double yTop, double width, double height);

  // void setDrawRectangle(x,y,w,h)

  //-----------------------------------------------------------------------------------------------
  // \name Drawing with juce::Graphics

  /** Draws the coordinate system onto the given Graphics object within the given rectangle
  (it's not yet tested, if the rectangle bounds actually work).  */
  virtual void drawPlot(Graphics& g);

  // todo:

  /** Draws the grids, coordinate axes and axis annotations, i.e. everything that appears behind 
  the actual curves or scattered dots or whatever is going to be plotted to represent the data. */
  virtual void drawPlotBackground(Graphics& g);
    // rename to drawBackground

  /** Draws the caption which should appear in front of the plotted curves. */
  virtual void drawPlotForeground(Graphics& g);

  /** Plots a curve by connecting the given datapoints with lines. */
  template<class T>
  void drawWithLines(Graphics& g, int numValues, T* valuesX, T* valuesY);

  /** Draws the given function by sampling it at x-values corresponding to each pixel's 
  model-coordinate. You can also set a sampling finer or coarser sampling interval by passing a
  (pixel) increment between sampling points other than unity. */
  template<class T>
  void drawWithLines(Graphics& g, std::function<T(T)>& function, double increment = 1.0, 
    float thickness = 2.f);

  /** Like drawWithLines drawWithLines(Graphics&, std::function<T(T)>&, double, float), but with an
  additional array of special x-values at which the function should be evaluated, regardless 
  whether or not one of the regular sampling points falls onto one of them. The idea is that some 
  function may have special key values, such as poles (resonances) or zeros (notches) at which we
  always want a datapoint because otherwise, there will be artifacts in the graphical 
  representation (such as showing a resonance with wrong amplitude). */
  template<class T>
  void drawWithLines(Graphics& g, std::function<T(T)>& function, 
    const std::vector<T>& specialValues, double increment = 1.0, float thickness = 2.f);


  /** Fills the area between the function given by the datapoints and the x-axis. */
  template<class T>
  void fillFunction(Graphics& g, int numValues, T* valuesX, T* valuesY);
    // not yet tested, maybe factor out into rsDrawer

  /** Draws the given datapoints as dots like in a scatterplot. */
  template<class T>
  void drawAsDots(Graphics& g, int numValues, T* valuesX, T* valuesY, float dotSize = 5.f, 
    bool stemsX = false, bool stemsY = false);
    // not yet tested ...maybe factor out the stem-drawing into separate functions drawAsStemsX
    // drawAsStemsY ..more modular

  //virtual void drawFunction(std::function<double(double)> function, double increment = 1);
  // goes through the x-range with given increment (in pixels) and evaluates the given function
  // at the (unmapped) x-values and draws the function. very convenient but maybe not very 
  // efficient...but that depends on the situation

  // drawWithBars, drawWithStems, etc.



  //-----------------------------------------------------------------------------------------------
  // \name Drawing on an SVG

  /** Analog to the other drawPlot version but draws onto an svg draing instead of a graphics 
  object. Useful for implementing export of plots to svg files. */
  virtual void drawPlotBackground(XmlElement* svg);
   // maybe move into subclass rsPlotDrawerWithSvg...or something
   // maybe have a subclass that uses OpenGL
   // or have an (possibly abstract) rsPlotDrawer baseclass and rsPlotDrawerNative, 
   // rsPlotDrawerOpenGL, rsPlotDrawerSVG subclasses

  template<class T>
  void drawWithLines(XmlElement* svg, int numValues, T* valuesX, T* valuesY);


protected:

  /** Sets up our mapper member according to the settings member and x, y, w, h variables. */
  void setupMapper();

  /** Returns the x-coordinate for the y-axis in model coordinates. */
  double getVerticalAxisX();

  /** Returns the y-coordinate for the x-axis in model coordinates. */
  double getHorizontalAxisY();


  const rsPlotSettings& settings;
  const PlotColourScheme& colors;

  RAPT::rsCoordinateMapper2D<double> mapper;

  double x, y, w, h;  // drawing rectangle in pixel coordinates...maybe rename

};

#endif