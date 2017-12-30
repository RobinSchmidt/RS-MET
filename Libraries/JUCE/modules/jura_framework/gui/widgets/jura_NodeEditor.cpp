rsDraggableNode::rsDraggableNode(rsNodeEditor* editor)
{
  jassert(editor != nullptr);
  nodeEditor = editor;
}

rsDraggableNode::~rsDraggableNode()
{
  if(paramX) paramX->deRegisterParameterObserver(this);
  if(paramY) paramY->deRegisterParameterObserver(this);
}

void rsDraggableNode::assignParameterX(Parameter* newParameterX)
{
  if(paramX)
    paramX->deRegisterParameterObserver(this);
  paramX = newParameterX;
  if(paramX)
    paramX->registerParameterObserver(this);
}

void rsDraggableNode::assignParameterY(Parameter* newParameterY)
{
  if(paramY)
    paramY->deRegisterParameterObserver(this);
  paramY = newParameterY;
  if(paramY)
    paramY->registerParameterObserver(this);
}

void rsDraggableNode::parameterChanged(Parameter* p)
{
  nodeEditor->nodeChanged(this);
}

//=================================================================================================

rsNodeEditor::rsNodeEditor()
{
  notifyPreSmoothing(true);
  notifyPostSmoothing(true);
}

rsNodeEditor::~rsNodeEditor()
{

}

void rsNodeEditor::parameterChanged(Parameter* p)
{
  //repaintOnMessageThread();
}

void rsNodeEditor::paint(Graphics& g)
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

void rsNodeEditor::mouseDown(const MouseEvent& e)
{

}

void rsNodeEditor::mouseDrag(const MouseEvent& e)
{

}

void rsNodeEditor::nodeChanged(const rsDraggableNode* node)
{
  repaintOnMessageThread();
}
