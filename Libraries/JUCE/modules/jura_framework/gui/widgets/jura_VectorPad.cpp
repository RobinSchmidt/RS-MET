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
