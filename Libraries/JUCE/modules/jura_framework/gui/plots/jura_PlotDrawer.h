#ifndef jura_PlotDrawer_h
#define jura_PlotDrawer_h


/** This class is used for drawing a plot. */

class JUCE_API rsPlotDrawer
{

public:


  rsPlotDrawer(const rsPlotSettings& plotSettings, const PlotColourScheme& colorScheme);


  virtual void drawPlot(Graphics& g, double x, double y, double w, double h);


  virtual void drawPlot(XmlElement* svg, double x, double y, double w, double h);


  //static void drawPlot(Graphics& g, const rsPlotSettings& settings, 
  //  const rsPlotColourScheme& colors, double x, double y, double w, double h);

protected:

  void setupMapper(double x, double y, double w, double h);

  const rsPlotSettings& settings;
  const PlotColourScheme& colors;

  RAPT::rsCoordinateMapper2D<double> mapper;

};

#endif