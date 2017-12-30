rsDraggableNode::rsDraggableNode(rsNodeEditor* editor, double x, double y)
{
  jassert(editor != nullptr);
  nodeEditor = editor;
  pixelX = x;
  pixelY = y;
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
  for(size_t i = 0; i < nodes.size(); i++)
    delete nodes[i];
}

void rsNodeEditor::addNode(double pixelX, double pixelY)
{
  nodes.push_back(new rsDraggableNode(this, pixelX, pixelY));
}

void rsNodeEditor::parameterChanged(Parameter* p)
{
  //repaintOnMessageThread();
}

void rsNodeEditor::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getHandleColour());
  drawNodes(g);
}

void rsNodeEditor::mouseDown(const MouseEvent& e)
{
  addNode((float)e.x, (float)e.y); 
  // todo: figure out, if there's currently a node under the mouse - if so, don't add a new node
  // but let the user drag around the existing node

  repaint();
}

void rsNodeEditor::mouseDrag(const MouseEvent& e)
{

}

void rsNodeEditor::nodeChanged(const rsDraggableNode* node)
{
  repaintOnMessageThread();
}

void rsNodeEditor::drawNodes(Graphics& g)
{
  for(size_t i = 0; i < nodes.size(); i++)
  {
    float x, y;
    x = (float)nodes[i]->pixelX;
    y = (float)nodes[i]->pixelY;
    g.fillEllipse(x-0.5f*dotSize, y-0.5f*dotSize, dotSize, dotSize);
  }
}
