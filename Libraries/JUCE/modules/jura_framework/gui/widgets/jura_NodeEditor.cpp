rsNodeEditor::rsNodeEditor()
{
  notifyPreSmoothing(true);
  notifyPostSmoothing(true);
}

rsNodeEditor::~rsNodeEditor()
{

}

void rsVectorPad::parameterChanged(Parameter* p)
{
  repaintOnMessageThread();
}

void rsVectorPad::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  /*
  float x, y;
  x = (float) RAPT::rsLinToLin(paramX->getValue(),  xMin, xMax, 0.0, double(getWidth()-1));
  y = (float) RAPT::rsLinToLin(paramY->getValue(), yMin, yMax, double(getHeight()-1), 0.0);
  g.setColour(getHandleColour());
  g.fillEllipse(x-0.5f*dotSize, y-0.5f*dotSize, dotSize, dotSize);
  */
}

void rsVectorPad::mouseDown(const MouseEvent& e)
{

}

void rsVectorPad::mouseDrag(const MouseEvent& e)
{

}
