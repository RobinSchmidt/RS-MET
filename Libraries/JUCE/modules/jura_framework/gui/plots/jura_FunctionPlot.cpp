
void rsFunctionPlot::addFunction(std::function<double(double)> function)
{
  append(functions, function);
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
    // todo: select color
    drawer.drawWithLines(g, functions[i]);
  }
}