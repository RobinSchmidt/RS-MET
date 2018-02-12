#ifndef jura_FunctionPlot_h
#define jura_FunctionPlot_h

/** A class to plot one or more functions based on std::function. You can add functions that should
be plotted by calling addFunction.

\todo maybe factor out a class rsCurvePlot or rsGraphPlot that has the stuff that is common to
data and function plots...or maybe integrate that stuff in rsPlot or rsPlotSettings
*/

class JUCE_API rsFunctionPlot : public rsPlot
{

public:

  rsFunctionPlot() {}
  virtual ~rsFunctionPlot() {}

  virtual void paint(Graphics &g) override;
  //virtual void resized() override;


  virtual void addFunction(std::function<double(double)> function);


protected:

  std::vector<std::function<double(double)>> functions;

};

#endif