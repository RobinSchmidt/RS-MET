rsVectorPad::rsVectorPad()
{
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
  g.fillAll(getBackgroundColour());
  float x, y;
  x = (float) RAPT::rsLinToLin(paramX->getValue(),  xMin, xMax, 0.0, double(getWidth()-1));
  y = (float) RAPT::rsLinToLin(paramY->getValue(), yMin, yMax, double(getHeight()-1), 0.0);
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

void rsVectorPad::setParametersFromMouseEvent(const MouseEvent& e)
{
  double x, y;
  x = RAPT::rsLinToLin(double(e.x), 0.0, double(getWidth()-1),  xMin, xMax);
  y = RAPT::rsLinToLin(double(e.y), double(getHeight()-1), 0.0, yMin, yMax);
  setParametersXY(x, y);
}

void rsVectorPad::setParametersXY(double x, double y)
{
  paramX->setValue(x, true, true);
  paramY->setValue(y, true, true);
}