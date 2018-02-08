rsPlotDrawer::rsPlotDrawer(const rsPlotSettings& plotSettings,
  const PlotColourScheme& colorScheme)
  : settings(plotSettings), colors(colorScheme)
{

}

void rsPlotDrawer::drawPlot(Graphics& g, double x, double y, double w, double h)
{
  mapper.mapperX.setLogScaled(settings.logScaledX);
  mapper.mapperY.setLogScaled(settings.logScaledY);
  mapper.setInputRange(
    settings.currentRange.getMinX(), settings.currentRange.getMaxX(),
    settings.currentRange.getMinY(), settings.currentRange.getMaxY());
  mapper.setOutputRange(x+0.5, jmax(x+w-0.5, x+1.0), jmax(y+h-0.5, y+1.0), y+0.5);

  if(settings.horizontalFineGridIsVisible)
  {
    g.setColour(colors.fineGrid);
    drawHorizontalGrid(g, mapper, settings.horizontalFineGridInterval, 1.f);

    //drawHorizontalGrid(g, settings.horizontalFineGridInterval, settings.logScaledY,
    //  colors.fineGrid, 1.0f);
  }




  int dummy = 0;
}