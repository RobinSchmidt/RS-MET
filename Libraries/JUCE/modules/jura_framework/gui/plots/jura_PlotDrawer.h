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

  /** Draws the coordinate system onto the given Graphics object within the given rectangle
  (it's not yet tested, if the rectangle bounds actually work).  */
  virtual void drawPlot(Graphics& g, double x, double y, double width, double height);

  /** Analog to the other drawPlot version but draws onto an svg draing instead of a graphics 
  object. Useful for implementing export of plots to svg files. */
  virtual void drawPlot(XmlElement* svg, double x, double y, double w, double h);


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