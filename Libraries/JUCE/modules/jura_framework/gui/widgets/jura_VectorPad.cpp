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

void rsVectorPad::paint(Graphics& g)
{
  if(paintBackground)
    g.fillAll(getBackgroundColour());
  float x, y;

  x = (float) RAPT::rsLinToLin(paramX->getNormalizedValue(),  0, 1, 0.5, double(getWidth()-0.5));
  y = (float) RAPT::rsLinToLin(paramY->getNormalizedValue(), 0, 1, double(getHeight()-0.5), 0.5);

    // maybe use a function normalizedToPixelCoords(&x, &y) . we should also use 
    // 0.5...getWidth()-0.5 to be consistent with the plots

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
  // todo: use margins, clip
  return RAPT::rsLinToLin(x, 0.5, getWidth()-0.5,  0, 1);
}

double rsVectorPad::pixelToNormalizedY(double y)
{
  return RAPT::rsLinToLin(y, getHeight()-0.5, 0.5, 0, 1);
}

double rsVectorPad::normalizedToPixelX(double x)
{
  return RAPT::rsLinToLin(x, 0, 1, 0.5, getWidth()-0.5);
}

double rsVectorPad::normalizedToPixelY(double y)
{
  return RAPT::rsLinToLin(y, 0, 1, getHeight()-0.5, 0.5);
}

void rsVectorPad::setParametersFromMouseEvent(const MouseEvent& e)
{
  double x, y;

  x = RAPT::rsLinToLin(double(e.x), 0.5, double(getWidth()-0.5),  0, 1);
  y = RAPT::rsLinToLin(double(e.y), double(getHeight()-0.5), 0.5, 0, 1);

  setParametersXY(x, y);
}

void rsVectorPad::setParametersXY(double x, double y) // rename to setNormalizedParameters
{
  paramX->setNormalizedValue(clip(x, 0, 1), true, true);
  paramY->setNormalizedValue(clip(y, 0, 1), true, true);
}