rsPlotDrawer::rsPlotDrawer(const rsPlotSettings& plotSettings,
  const PlotColourScheme& colorScheme, double xLeft, double yTop, double width, double height)
  : settings(plotSettings), colors(colorScheme)
{
  x = xLeft;
  y = yTop;
  w = width;
  h = height;
  setupMapper();
}

void rsPlotDrawer::drawPlot(Graphics& g)
{
  // split into drawPlotBackground / drawPlotForeground

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
    drawAxisY(g, mapper, getVerticalAxisX(),   settings.axisLabelY, colors.text); }
  if(settings.axisPositionX != rsPlotSettings::INVISIBLE 
    && settings.axisValuesPositionX != rsPlotSettings::NO_ANNOTATION) {
    g.setColour(colors.axes);
    drawAxisValuesX(g, mapper, settings.verticalCoarseGridInterval, getHorizontalAxisY(), 
      settings.stringConversionForAxisX, colors.text); }
  if(settings.axisPositionY != rsPlotSettings::INVISIBLE 
    && settings.axisValuesPositionY != rsPlotSettings::NO_ANNOTATION)  {
    g.setColour(colors.axes);
    jura::drawAxisValuesY(g, mapper, settings.horizontalCoarseGridInterval, 
      getVerticalAxisX(), settings.stringConversionForAxisY, colors.text); }





  // caption/headline (positioning formulas needs test):
  static const BitmapFont *font = &BitmapFontRoundedBoldA10D0::instance;
  float cw = (float) font->getTextPixelWidth(settings.captionString);
  float ch = (float) font->getFontHeight();
  switch( settings.captionPosition )
  {
  case rsPlotSettings::NO_CAPTION: break;
  case rsPlotSettings::TOP_CENTER: {
    drawBitmapText(g, settings.captionString, x+0.5f*w-0.5f*cw, y+16, w, 16, font, 
      Justification::centred, colors.text);
  } break;
  case rsPlotSettings::CENTER: {
    float cx = (float) (x + (w-cw) * .5);
    float cy = (float) (y + (h-ch) * .5);
    drawBitmapText(g, settings.captionString, cx, cy, w, 16, font, 
      Justification::topLeft, colors.text);
  } break;
  }

  // outline:
  g.setColour(colors.outline);
  g.drawRect(float(x), float(y), float(w), float(h), 2.f);
}

void rsPlotDrawer::drawPlotBackground(Graphics& g)
{

}

void rsPlotDrawer::drawPlotForeground(Graphics& g)
{

}

void rsPlotDrawer::drawWithLines(Graphics& g, int numValues, float* valuesX, float* valuesY)
{

}

void rsPlotDrawer::drawAsDots(Graphics& g, int numValues, float* valuesX, float* valuesY)
{

}

//-------------------------------------------------------------------------------------------------
// SVG stuff:

void rsPlotDrawer::drawPlot(XmlElement* svg)
{
  setupMapper();

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
  if(settings.axisPositionX != rsPlotSettings::INVISIBLE 
    && settings.axisValuesPositionX != rsPlotSettings::NO_ANNOTATION) {
    drawAxisValuesX(svg, mapper, settings.verticalCoarseGridInterval, 
      getHorizontalAxisY(), settings.stringConversionForAxisX, colors.axes); }
  if(settings.axisPositionY != rsPlotSettings::INVISIBLE 
    && settings.axisValuesPositionY != rsPlotSettings::NO_ANNOTATION) {
    drawAxisValuesY(svg, mapper, settings.horizontalCoarseGridInterval, 
      getVerticalAxisX(), settings.stringConversionForAxisY, colors.axes); }
}

//-------------------------------------------------------------------------------------------------

void rsPlotDrawer::setupMapper()
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
