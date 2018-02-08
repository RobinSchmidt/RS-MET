rsPlotDrawer::rsPlotDrawer(const rsPlotSettings& plotSettings,
  const PlotColourScheme& colorScheme)
  : settings(plotSettings), colors(colorScheme)
{

}

void rsPlotDrawer::drawPlot(Graphics& g, double x, double y, double w, double h)
{
  setupMapper(x, y, w, h);

  // fine grids:
  if(settings.horizontalFineGridIsVisible) {
    g.setColour(colors.fineGrid);
    drawHorizontalGrid(g, mapper, settings.horizontalFineGridInterval, 1.f); }
  if(settings.verticalFineGridIsVisible) {
    g.setColour(colors.fineGrid);
    drawVerticalGrid(g, mapper, settings.verticalFineGridInterval, 1.f); }
  if(settings.radialFineGridIsVisible) {
    g.setColour(colors.fineGrid);
    drawRadialGrid(g, mapper, settings.radialFineGridInterval, 1.f); }
  if(settings.angularFineGridIsVisible) {
    g.setColour(colors.fineGrid);
    drawAngularGrid(g, mapper, settings.angularFineGridInterval, 1.f); }

  // coarse grids:
  if(settings.horizontalCoarseGridIsVisible) {
    g.setColour(colors.coarseGrid);
    drawHorizontalGrid(g, mapper, settings.horizontalCoarseGridInterval, 1.f); }
  if(settings.verticalCoarseGridIsVisible) {
    g.setColour(colors.coarseGrid);
    drawVerticalGrid(g, mapper, settings.verticalCoarseGridInterval, 1.f); }
  if(settings.radialCoarseGridIsVisible) {
    g.setColour(colors.coarseGrid);
    drawRadialGrid(g, mapper, settings.radialCoarseGridInterval, 1.f); }
  if(settings.angularCoarseGridIsVisible) {
    g.setColour(colors.coarseGrid);
    drawAngularGrid(g, mapper, settings.angularCoarseGridInterval, 1.f); }

  // coordinate axes:
  if(settings.axisPositionX != rsPlotSettings::INVISIBLE) {
    g.setColour(colors.axes);
    drawAxisX(g, mapper, getHorizontalAxisY(), settings.axisLabelX, colors.text); }
  if(settings.axisPositionY != rsPlotSettings::INVISIBLE) {
    g.setColour(colors.axes);
    jura::drawAxisY(g, mapper, getVerticalAxisX(), settings.axisLabelY, colors.text); }


  int dummy = 0;
}

void rsPlotDrawer::drawPlot(XmlElement* svg, double x, double y, double w, double h)
{
  setupMapper(x, y, w, h);

  // fine grids:
  if(settings.horizontalFineGridIsVisible)
    drawHorizontalGrid(svg, mapper, settings.horizontalFineGridInterval, 1.f, colors.fineGrid);
  if(settings.verticalFineGridIsVisible)
    drawVerticalGrid(svg, mapper, settings.verticalFineGridInterval, 1.f, colors.fineGrid);
  if(settings.radialFineGridIsVisible)
    drawRadialGrid(svg, mapper, settings.radialFineGridInterval, 1.f, colors.fineGrid);
  if(settings.angularFineGridIsVisible)
    drawAngularGrid(svg, mapper, settings.angularFineGridInterval, 1.f, colors.fineGrid);

  // coarse grids:
  if(settings.horizontalCoarseGridIsVisible)
    drawHorizontalGrid(svg, mapper, settings.horizontalCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.verticalCoarseGridIsVisible)
    drawVerticalGrid(svg, mapper, settings.verticalCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.radialCoarseGridIsVisible)
    drawRadialGrid(svg, mapper, settings.radialCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.angularCoarseGridIsVisible)
    drawAngularGrid(svg, mapper, settings.angularCoarseGridInterval, 1.f, colors.coarseGrid);

  // coordinate axes:
  if(settings.axisPositionX != rsPlotSettings::INVISIBLE)
    drawAxisX(svg, mapper, getHorizontalAxisY(), settings.axisLabelX, colors.axes);
  if(settings.axisPositionY != rsPlotSettings::INVISIBLE) 
    drawAxisY(svg, mapper, getVerticalAxisX(),   settings.axisLabelY, colors.axes);


  int dummy = 0;
}

void rsPlotDrawer::setupMapper(double x, double y, double w, double h)
{
  mapper.mapperX.setLogScaled(settings.logScaledX);
  mapper.mapperY.setLogScaled(settings.logScaledY);
  mapper.setInputRange(
    settings.currentRange.getMinX(), settings.currentRange.getMaxX(),
    settings.currentRange.getMinY(), settings.currentRange.getMaxY());
  mapper.setOutputRange(x+0.5, jmax(x+w-0.5, x+1.0), jmax(y+h-0.5, y+1.0), y+0.5);
}

double rsPlotDrawer::getVerticalAxisX()
{
  double x = 0.0;
  if( settings.axisPositionY == rsPlotSettings::LEFT )
    x = mapper.unmapX(jmin(8., mapper.getOutMaxX())); // maybe use margin parameter instead of 8
  else if( settings.axisPositionY == rsPlotSettings::RIGHT )
    x = mapper.unmapX(jmax(mapper.getOutMaxX()-8, 0.));
  return x;
}

double rsPlotDrawer::getHorizontalAxisY()
{
  double y = 0.0;
  if( settings.axisPositionX == rsPlotSettings::BOTTOM )
    y = mapper.unmapY(jmax(mapper.getOutMinY()-8, 0.)); // is outMinX bcs of reversal of y-axis
  else if( settings.axisPositionX == rsPlotSettings::TOP )
    y = mapper.unmapY(jmin(8., mapper.getOutMinY()));
  return y;
}
