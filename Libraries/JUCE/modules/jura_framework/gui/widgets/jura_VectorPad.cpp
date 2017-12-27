rsVectorPad::rsVectorPad()
{

}

rsVectorPad::~rsVectorPad()
{
  if(paramX != nullptr) paramX->deRegisterParameterObserver(this);
  if(paramY != nullptr) paramY->deRegisterParameterObserver(this);
}

void rsVectorPad::assignParameterX(Parameter* newParameterX)
{
  paramX = newParameterX;
  paramX->registerParameterObserver(this);
}

void rsVectorPad::assignParameterY(Parameter* newParameterY)
{
  paramY = newParameterY;
  paramY->registerParameterObserver(this);
}

void rsVectorPad::parameterChanged(Parameter* p)
{
  repaintOnMessageThread();
}

void rsVectorPad::paint(Graphics& g)
{
  g.fillAll(Colours::black); // use background color
}

void rsVectorPad::mouseDown(const MouseEvent& e)
{
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