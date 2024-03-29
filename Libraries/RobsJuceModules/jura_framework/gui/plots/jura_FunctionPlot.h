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



  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Puts the object back into it's initial state which it has after construction. */
  virtual void init();

  /** Adds a function to be plotted. It will receive the x-value (in model-coordinates) and should
  return the corresponding y-value (also in model-coordinates). For example, it could receive a 
  frequency in Hz and output a magnitude level in decibels. */
  virtual void addFunction(std::function<double(double)> function);

  /** When plotting a function, we sample the function output at x-values that correspond to the 
  pixel x-coordinates. For some functions with certain localized features such as poles/resonances 
  or zeros/notches, our sampling may miss these features leading to artifacts in the plot, such as 
  a wrong resonance amplitude. With this function, you can set some additional evaluation points 
  that may not be covered (well enough) by the regular smapling. */
  virtual void setSpecialEvaluationPoint(size_t funcIndex, size_t pointIndex, double xValue);

  /** Sets up the plot sttings for a nice looking frequency response plot with a logarithmic 
  frequency axis (as x-axis) and a linear dB y-axis (of course, dB itself is a logarithmic measure,
  "linear" here means that the dB values are shown lineraly, so it's also logarithmic amplitude 
  axis at the end of the day). */
  virtual void setupForDecibelsAgainstLogFrequency(double minFreq, double maxFreq, 
    double minDb, double maxDb, double yGridSpacing = 6);

  /** Sets the number of functions that should be shown in the plot. The idea is that you may not 
  want to show all of those functions that have been added via addFunction all the time but only 
  the first few. */
  virtual void setNumFunctionsToPlot(size_t newNumber) { numFunctionsToPlot = newNumber; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  virtual Colour getGraphColor(size_t index);

  virtual float getGraphThickness(size_t index);

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  virtual void paint(Graphics &g) override;
  //virtual void resized() override;




protected:

  std::vector<std::function<double(double)>> functions;
  std::vector<std::vector<double>> specialPoints;
  size_t numFunctionsToPlot = 0;

};

#endif