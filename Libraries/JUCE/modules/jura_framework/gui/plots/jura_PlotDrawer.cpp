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


  int dummy = 0;
}

void rsPlotDrawer::drawPlot(XmlElement* svg, double x, double y, double w, double h)
{
  setupMapper(x, y, w, h);

  if(settings.horizontalFineGridIsVisible)
    drawHorizontalGrid(svg, mapper, settings.horizontalFineGridInterval, 1.f, colors.fineGrid);
  if(settings.verticalFineGridIsVisible)
    drawVerticalGrid(svg, mapper, settings.verticalFineGridInterval, 1.f, colors.fineGrid);
  if(settings.radialFineGridIsVisible)
    drawRadialGrid(svg, mapper, settings.radialFineGridInterval, 1.f, colors.fineGrid);
  if(settings.angularFineGridIsVisible)
    drawAngularGrid(svg, mapper, settings.angularFineGridInterval, 1.f, colors.fineGrid);

  if(settings.horizontalCoarseGridIsVisible)
    drawHorizontalGrid(svg, mapper, settings.horizontalCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.verticalCoarseGridIsVisible)
    drawVerticalGrid(svg, mapper, settings.verticalCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.radialCoarseGridIsVisible)
    drawRadialGrid(svg, mapper, settings.radialCoarseGridInterval, 1.f, colors.coarseGrid);
  if(settings.angularCoarseGridIsVisible)
    drawAngularGrid(svg, mapper, settings.angularCoarseGridInterval, 1.f, colors.coarseGrid);


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