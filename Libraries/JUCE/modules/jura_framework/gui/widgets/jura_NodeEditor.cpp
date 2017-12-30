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

void rsDraggableNode::setPixelPosition(double newX, double newY)
{
  pixelX = newX;
  pixelY = newY;
  nodeEditor->nodeChanged(this);
  // todo: do not call nodeChanged directly here - instead, set up the paramX, paramY parameters
  // according to the new pixel position. this will trigger a call to our parameterChanged function
  // which will in turn call nodeEditor->nodeChanged(this); ...hmm...or maybe that's not so good
  // because nodeChanged will then get called twice (once for x and once for y)
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

// setup:

rsDraggableNode* rsNodeEditor::addNode(double pixelX, double pixelY)
{
  rsDraggableNode* newNode = new rsDraggableNode(this, pixelX, pixelY);
  nodes.push_back(newNode);
  return newNode;
}

void rsNodeEditor::removeNodeAt(int pixelX, int pixelY)
{
  int i = getNodeIndexAt(pixelX, pixelY);
  if(i != -1)
  {
    delete nodes[i];
    remove(nodes, i);
  }
}

void rsNodeEditor::setDotSize(float newDotSize)
{
  dotSize = newDotSize;
  repaint();
}

// inquiry:

rsDraggableNode* rsNodeEditor::getNoteAt(int pixelX, int pixelY)
{
  int i = getNodeIndexAt(pixelX, pixelY);
  if(i != -1)
    return nodes[i];
  return nullptr;
}

int rsNodeEditor::getNodeIndexAt(int pixelX, int pixelY)
{
  float x  = (float)pixelX;
  float y  = (float)pixelY;
  float r2 = (float)(dotSize*dotSize);  // radius^2 of circle to check
  for(size_t i = 0; i < nodes.size(); i++)
  {
    float dx = x - (float)nodes[i]->pixelX;
    float dy = y - (float)nodes[i]->pixelY;
    float d2 = dx*dx + dy*dy;
    if(d2 <= r2)
      return (int)i;
  }
  return -1;
}

// callbacks:

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
  draggedNode = getNoteAt(e.x, e.y);
  if(e.mods.isLeftButtonDown())
  {
    if(draggedNode == nullptr)
    {
      draggedNode = addNode((float)e.x, (float)e.y);
      repaint();
    }
  }
  else if(e.mods.isRightButtonDown())
  {
    removeNodeAt(e.x, e.y);
    draggedNode = nullptr;
    repaint();
  }
}

void rsNodeEditor::mouseDrag(const MouseEvent& e)
{
  if(draggedNode != nullptr)
    draggedNode->setPixelPosition(e.x, e.y); // will indirectly trigger a repaint
}

void rsNodeEditor::mouseUp(const MouseEvent &e)
{
  draggedNode = nullptr;
}

void rsNodeEditor::mouseMove(const MouseEvent &e)
{
  rsDraggableNode* nodeUnderMouse = getNoteAt(e.x, e.y);
  if(nodeUnderMouse != nullptr)
    setMouseCursor(MouseCursor(MouseCursor::PointingHandCursor));
  else
    setMouseCursor(MouseCursor(MouseCursor::NormalCursor));
}

void rsNodeEditor::nodeChanged(const rsDraggableNode* node)
{
  repaintOnMessageThread();
}

// misc:

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

//=================================================================================================

rsNodeBasedFunctionEditor::rsNodeBasedFunctionEditor(
  RAPT::rsInterpolatingFunction<double>* functionMapper)
{
  mapper = functionMapper;
}

void rsNodeBasedFunctionEditor::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getHandleColour());
  if(nodes.size() > 1)
  {
    for(size_t i = 0; i < nodes.size()-1; i++)
    {
      float x1, y1, x2, y2;
      x1 = (float)nodes[i]->getPixelX();
      y1 = (float)nodes[i]->getPixelY();
      x2 = (float)nodes[i+1]->getPixelX();
      y2 = (float)nodes[i+1]->getPixelY();
      g.drawLine(x1, y1, x2, y2, 2.f);
    }
  }
  drawNodes(g);
}

rsDraggableNode* rsNodeBasedFunctionEditor::addNode(double pixelX, double pixelY)
{
  double x, y;
   
  // preliminary - todo: apply mapping between pixel and model coordinates:
  x = pixelX;
  y = pixelY;

  if(mapper != nullptr)
  {
    mapper->addDataPoint(x, y);
    // maybe we need to set up callbacks here?
  }

  return rsNodeEditor::addNode(pixelX, pixelY);
}

void rsNodeBasedFunctionEditor::removeNodeAt(int pixelX, int pixelY)
{
  rsNodeEditor::removeNodeAt(pixelX, pixelY);
}

void rsNodeBasedFunctionEditor::nodeChanged(const rsDraggableNode* node)
{
  rsNodeEditor::nodeChanged(node);
}
