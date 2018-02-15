rsVectorPad::rsVectorPad()
{
  notifyPreSmoothing(true);
  notifyPostSmoothing(true);
  dummyParam = new Parameter("Dummy", -1.0, +1.0, 0.0);
  paramX = paramY = dummyParam;
}

rsVectorPad::~rsVectorPad()
{
  if(paramX != nullptr) paramX->deRegisterParameterObserver(this);
  if(paramY != nullptr) paramY->deRegisterParameterObserver(this);
  delete dummyParam;
}

void rsVectorPad::assignParameterX(Parameter* newParameterX)
{
  paramX->deRegisterParameterObserver(this);
  paramX = newParameterX;
  if(paramX)
    paramX->registerParameterObserver(this);
  else
    paramX = dummyParam;
}

void rsVectorPad::assignParameterY(Parameter* newParameterY)
{
  paramY->deRegisterParameterObserver(this);
  paramY = newParameterY;
  if(paramY)
    paramY->registerParameterObserver(this);
  else
    paramY = dummyParam;
  // get rid of duplication
}

void rsVectorPad::parameterChanged(Parameter* p)
{
  repaintOnMessageThread();
}

void rsVectorPad::adjustMarginsToPlotX(rsPlot* plot)
{
  // retrieve min/max values of paramX and map them to pixel positions according to the mapper
  // in the plot, from the pixel-positions, compute the margins

  jassert(getBounds() == plot->getBounds()); // makes sense only when plot is used as background

  double x = paramX->getMinValue();
  x = plot->getCoordinateMapper()->mapX(x);
  leftMargin = x - 0.5;

  x = paramX->getMaxValue();
  x = plot->getCoordinateMapper()->mapX(x);
  rightMargin = getWidth() - x - 0.5;
}

void rsVectorPad::adjustMarginsToPlotY(rsPlot* plot)
{
  jassert(getBounds() == plot->getBounds()); // makes sense only when plot is used as background

  double y = paramY->getMinValue();
  y = plot->getCoordinateMapper()->mapY(y);
  bottomMargin = getHeight() - y - 0.5;

  y = paramY->getMaxValue();
  y = plot->getCoordinateMapper()->mapY(y);
  topMargin = y - 0.5;
}

void rsVectorPad::setBackgroundPlot(rsPlot* newBackgroundPlot)
{
  jassertfalse; // not yet implemented
  /*
  if(backgroundPlot != nullptr)
    backgroundPlot->getPlotSettings()->deRegisterObserver(this);
  backgroundPlot = newBackgroundPlot;
  adjustMarginsToPlot(backgroundPlot);
  if(backgroundPlot != nullptr)
    backgroundPlot->getPlotSettings()->registerObserver(this);
  */
}

void rsVectorPad::paint(Graphics& g)
{
  if(paintBackground)
    g.fillAll(getBackgroundColour());
  float x = (float) normalizedToPixelX(paramX->getNormalizedValue()); 
  float y = (float) normalizedToPixelY(paramY->getNormalizedValue()); 
  g.setColour(getHandleColour());
  g.fillEllipse(x-0.5f*dotSize, y-0.5f*dotSize, dotSize, dotSize);
}

void rsVectorPad::mouseDown(const MouseEvent& e)
{
  if(e.mods.isCommandDown()) {
    paramX->resetToDefaultValue(true, true);
    paramY->resetToDefaultValue(true, true); }
  else
    setParametersFromMouseEvent(e);
}

void rsVectorPad::mouseDrag(const MouseEvent& e)
{
  setParametersFromMouseEvent(e);
}

double rsVectorPad::pixelToNormalizedX(double x)
{
  return RAPT::rsLinToLin(x, leftMargin + 0.5, getWidth() - rightMargin - 0.5, 0, 1);
}

double rsVectorPad::pixelToNormalizedY(double y)
{
  return RAPT::rsLinToLin(y, getHeight() - bottomMargin - 0.5, topMargin + 0.5, 0, 1);
}

double rsVectorPad::normalizedToPixelX(double x)
{
  return RAPT::rsLinToLin(x, 0, 1, leftMargin + 0.5, getWidth() - rightMargin - 0.5);
}

double rsVectorPad::normalizedToPixelY(double y)
{
  return RAPT::rsLinToLin(y, 0, 1, getHeight() - bottomMargin - 0.5, topMargin + 0.5);
}

void rsVectorPad::setParametersFromMouseEvent(const MouseEvent& e)
{
  setParametersXY(pixelToNormalizedX(e.x), pixelToNormalizedY(e.y));
}

void rsVectorPad::setParametersXY(double x, double y) // rename to setNormalizedParameters
{
  paramX->setNormalizedValue(clip(x, 0, 1), true, true);
  paramY->setNormalizedValue(clip(y, 0, 1), true, true);
}