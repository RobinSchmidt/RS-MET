#ifndef jura_PlotDrawer_h
#define jura_PlotDrawer_h


/** This class is used for drawing a plot. */

class JUCE_API rsPlotDrawer
{

public:

  static void drawPlot(Graphics& g, const rsPlotSettings& settings, 
    const rsPlotColourScheme& colors, double x, double y, double w, double h);

};

#endif