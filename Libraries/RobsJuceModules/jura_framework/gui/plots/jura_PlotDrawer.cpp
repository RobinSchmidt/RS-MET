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
  // background:
  fillRectWithBilinearGradient(g, 0, 0, int(w), int(h),
    colors.topLeft, colors.topRight, colors.bottomLeft, colors.bottomRight);

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

template<class T>
void rsPlotDrawer::drawWithLines(Graphics& g, int numValues, T* valuesX, T* valuesY)
{
  float thickness = 2.f; // make parameter
  for(int i = 0; i < numValues-1; i++) 
  {
    float x1 = (float) mapper.mapX(valuesX[i]);
    float y1 = (float) mapper.mapY(valuesY[i]);
    float x2 = (float) mapper.mapX(valuesX[i+1]);
    float y2 = (float) mapper.mapY(valuesY[i+1]);
    g.drawLine(x1, y1, x2, y2, thickness); // drawLine is more efficient than using a Path object
  }
}

template<class T>
void rsPlotDrawer::drawWithLines(Graphics& g, std::function<T(T)>& func, double inc, float thickness)
{
  double x1 = x;
  double y1 = mapper.mapY(func(mapper.unmapX(x1)));
  double x2, y2;
  while(x1 <= x+w)
  {
    x2 = x1 + inc;
    y2 = mapper.mapY(func(mapper.unmapX(x2)));
    g.drawLine((float)x1, (float)y1, (float)x2, (float)y2, thickness);
    x1 = x2;
    y1 = y2;
  }
  // maybe try to optimize with a function iterator that creates successive x-coordinates (in model
  // coordinates), saves one unmapX operation per pixel at the cost of a function iterator update 
  // (which is cheaper)
}

template<class T>
void rsPlotDrawer::drawWithLines(Graphics& g, std::function<T(T)>& func, 
  const std::vector<T>& specialValues, double inc, float thickness)
{
  if(specialValues.size() == 0)
    drawWithLines(g, func, inc, thickness);

  // suffixes: p: pixel-coordinate, m: model-coordinate
  size_t is  = 0;                 // index of special value
  double xsm = specialValues[is]; // current special value
  double x1p = x;
  double y1p = mapper.mapY(func( mapper.unmapX(x1p) ));
  double x2p, x2m, y2p;
  while(x1p <= x+w)
  {
    x2p = x1p + inc;
    x2m = mapper.unmapX(x2p);

    if(x2m >= xsm)
    {
      double xsp = mapper.mapX(xsm);
      double ysp = mapper.mapY(func(xsm));
      g.drawLine((float)x1p, (float)y1p, (float)xsp, (float)ysp, thickness);
      is++;
      if(is < specialValues.size())
        xsm = specialValues[is];
      else
        xsm = INF;
      x1p = xsp; 
      y1p = ysp;

      // debug-stuff:
      //g.drawLine(0, ysp, 500, ysp, 1.f);
      //continue; // test
    }

    y2p = mapper.mapY(func(x2m));
    g.drawLine((float)x1p, (float)y1p, (float)x2p, (float)y2p, thickness);
    x1p = x2p;
    y1p = y2p;
  }
}

template<class T>
void rsPlotDrawer::fillFunction(Graphics& g, int N, T* x, T* y)
{
  T zero = 0;  // todo: allow a baseline other than zero
  Path path;
  path.startNewSubPath((float) mapper.mapX(x[0]), (float) mapper.mapY(y[0]));
  for(int i = 1; i < N; i++)
    path.lineTo((float) mapper.mapX(x[i]), (float) mapper.mapY(y[i]));
  path.lineTo((float) mapper.mapX(x[N-1]), (float) mapper.mapY(zero));
  path.lineTo((float) mapper.mapX(x[0]),   (float) mapper.mapY(zero));
  path.lineTo((float) mapper.mapX(x[0]),   (float) mapper.mapY(y[0]));
  path.closeSubPath();
  g.fillPath(path);
}

template<class T>
void rsPlotDrawer::drawAsDots(Graphics& g, int numValues, T* valuesX, T* valuesY, float size, 
  bool stemsX, bool stemsY)
{
  // make parameters - to be used to draw lines to x- and/or y-axis:
  float x0 = (float) mapper.mapX(0);
  float y0 = (float) mapper.mapY(0);

  //float size  = 4.f; // make parameter
  float size2 = 0.5f * size;

  float x, y;	   // current x and y value
  for(int i = 0; i < numValues; i++)
  {
    // read out the tables:
    x = (float) mapper.mapX(valuesX[i]);
    y = (float) mapper.mapY(valuesY[i]);

    // add a dot at postion x, y:
    g.fillEllipse(x-size2, y-size2, size, size);

    // draw lines to x- and y-axis if the option is selected (todo: move into separate function
    // drawAsStemsToX, drawAsStemsToY):
    if(stemsX) g.drawLine(x,  y, x, y0);
    if(stemsY) g.drawLine(x0, y, x, y);
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

template<class T>
void rsPlotDrawer::drawWithLines(XmlElement* svg, int numValues, T* valuesX, T* valuesY)
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
    settings.getCurrentRangeMinX(), settings.getCurrentRangeMaxX(),
    settings.getCurrentRangeMinY(), settings.getCurrentRangeMaxY());
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
