
void rsFunctionPlot::addFunction(std::function<double(double)> function)
{
  append(functions, function);
}

void rsFunctionPlot::setupForDecibelsAgainstLogFrequency(double minFreq, double maxFreq,
  double minDb, double maxDb, double yGridSpacing)
{
  // maybe switch off re-rendering during setup?

  setHorizontalCoarseGrid(yGridSpacing,  true);
  setHorizontalFineGrid(  1.0,           false);
  setVerticalCoarseGrid(  2.0,           true);
  setVerticalFineGrid(pow(2.0, 1.0/3.0), false);
  setMaximumRange(minFreq, maxFreq, minDb, maxDb);
  setCurrentRange(minFreq, maxFreq, minDb, maxDb);
  useLogarithmicScaleX(true);
  useLogarithmicScaleY(false);
  setAxisPositionX(rsPlotSettings::BOTTOM);
  setAxisLabelX(String::empty);
  setAxisPositionY(rsPlotSettings::LEFT);
  setAxisValuesPositionY(rsPlotSettings::LEFT);
  setAxisLabelY(String::empty);
  //setStringConversionForInfoLineX(hertzToStringWithUnitTotal5);  
  //setStringConversionForInfoLineX(rojue::frequencyToNoteString);  
  setStringConversionForInfoLineX(frequencyInHzAndAsNote); 
  setStringConversionForInfoLineY(decibelsToStringWithUnit2);
}

void rsFunctionPlot::setSpecialEvaluationPoint(size_t funcIndex, size_t pointIndex, double xValue)
{
  //jassert(funcIndex  < specialPoints.size());            // funcIndex out of range
  //jassert(pointIndex < specialPoints[funcIndex].size()); // pointIndex out of range

  while(funcIndex >= specialPoints.size())
    specialPoints.push_back(std::vector<double>());
  while(pointIndex >= specialPoints[funcIndex].size())
    specialPoints[funcIndex].push_back(0.0);
  specialPoints[funcIndex][pointIndex] = xValue;
}

Colour rsFunctionPlot::getGraphColour(int index)
{
  return Colours::white;  // preliminary
}

float rsFunctionPlot::getGraphThickness(int index)
{
  return 2.f;
}

void rsFunctionPlot::paint(Graphics &g)
{
  rsPlot::paint(g); // draws background (axes, grids, etc.)

  // todo: create drawer object and draw the functions (maybe, we should have an rsPlotDrawer 
  // member, so we don't have to create a new object in each subclass in order to do the subclass
  // specific additional drawing

  rsPlotDrawer drawer(plotSettings, plotColourScheme, 0, 0, getWidth(), getHeight());
  for(size_t i = 0; i < functions.size(); i++)
  {
    g.setColour(getGraphcColour(i));
    if(specialPoints.size() <= i)
      drawer.drawWithLines(g, functions[i], 1.0, getGraphThickness(i));
    else
      drawer.drawWithLines(g, functions[i], specialPoints[i], 1.0, getGraphThickness(i));
    // maybe change order of thickness and increment - thickness is more likely to be customized
  }
}