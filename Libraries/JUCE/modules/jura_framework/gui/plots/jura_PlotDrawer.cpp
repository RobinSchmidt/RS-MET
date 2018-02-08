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
  // this function doesn't make much sense anymore unless we use a template method approach to draw
  // the data in between these two calls:
  drawPlotBackground(g);
  //drawPlotContent(g);  // could be a template method
  drawPlotForeground(g);
}

void rsPlotDrawer::drawPlotBackground(Graphics& g)
{
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
}

void rsPlotDrawer::drawPlotForeground(Graphics& g)
{
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

void rsPlotDrawer::drawWithLines(Graphics& g, int numValues, float* valuesX, float* valuesY)
{
  float thickness = 2.f; // make parameter
  for(int i = 0; i < numValues-1; i++) 
  {
    float x1 = (float) mapper.mapX(valuesX[i]);
    float y1 = (float) mapper.mapY(valuesY[i]);
    float x2 = (float) mapper.mapX(valuesX[i+1]);
    float y2 = (float) mapper.mapY(valuesY[i+1]);
    g.drawLine(x1, y1, x2, y2, thickness);
  }
}

void rsPlotDrawer::filledFunction(Graphics& g, int N, float* x, float* y)
{
  // not yet tested
  Path path;
  path.startNewSubPath((float) mapper.mapX(x[0]), (float) mapper.mapY(y[0]));
  for(int i = 1; i < N; i++)
    path.lineTo((float) mapper.mapX(x[i]), (float) mapper.mapY(y[i]));
  path.lineTo((float) mapper.mapX(x[N-1]), (float) mapper.mapY(0));
  path.lineTo((float) mapper.mapX(x[0]),   (float) mapper.mapY(0));
  path.lineTo((float) mapper.mapX(x[0]),   (float) mapper.mapY(y[0]));
  path.closeSubPath();
  g.fillPath(path);
}

void rsPlotDrawer::drawAsDots(Graphics& g, int numValues, float* valuesX, float* valuesY)
{
  //g.setColour(graphColor); should be done ouside
  // todo: make dot-size adjustable, templatize to make it work for double (and maybe int), too

  // make parameters - to be used to draw lines to x- and/or y-axis:
  bool lineToAxisX = false;
  bool lineToAxisY = false;
  float x0 = (float) mapper.mapX(0);
  float y0 = (float) mapper.mapY(0);

  float size  = 4.f; // make parameter
  float size2 = 0.5f * size;

  float x, y;	   // current x and y value
  for(int i = 0; i < numValues; i++)
  {
    // read out the tables:
    x = (float) mapper.mapX(valuesX[i]);
    y = (float) mapper.mapY(valuesY[i]);

    // add a dot at postion x, y:
    //g.fillEllipse(x-1, y-1, 3, 3); // is this correct?
    g.fillEllipse(x-size2, y-size2, size, size); // is this correct?

    // draw lines to x- and y-axis if the option is selected:
    if(lineToAxisX) g.drawLine(x,  y, x, y0); // needs testing
    if(lineToAxisY) g.drawLine(x0, y, x, y);
  }
}

//-------------------------------------------------------------------------------------------------
// SVG stuff:

void rsPlotDrawer::drawPlotBackground(XmlElement* svg)
{
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

void rsPlotDrawer::drawWithLines(XmlElement* svg, int numValues, float* valuesX, float* valuesY)
{
  String pathString;
  for(int i = 0; i < numValues; i++)
  {
    float x = (float) mapper.mapX(valuesX[i]);
    float y = (float) mapper.mapY(valuesY[i]);
    pathString += String(x) + " " + String(y) + ", ";
  }
  Colour colour = Colours::black; // preliminary
  XmlElement* curvePath = new XmlElement("polyline");
  curvePath->setAttribute("points", pathString);
  curvePath->setAttribute("style", "stroke-width: " + String(1.0) + "; stroke: #" 
    + colour.toString().substring(2) + "; fill: none;" );
  svg->addChildElement(curvePath);
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
