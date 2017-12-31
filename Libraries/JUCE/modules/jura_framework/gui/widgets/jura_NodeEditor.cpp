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

void rsDraggableNode::setPixelPosition(double newX, double newY, bool callNodeChanged)
{
  pixelX = newX;
  pixelY = newY;
  if(callNodeChanged)
    nodeEditor->nodeChanged(index);
  //nodeEditor->nodeChanged(this);
  // todo: do not call nodeChanged directly here - instead, set up the paramX, paramY parameters
  // according to the new pixel position. this will trigger a call to our parameterChanged function
  // which will in turn call nodeEditor->nodeChanged(this); ...hmm...or maybe that's not so good
  // because nodeChanged will then get called twice (once for x and once for y)
}

void rsDraggableNode::parameterChanged(Parameter* p)
{
  nodeEditor->nodeChanged(this->index);
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
/*
rsDraggableNode* rsNodeEditor::addNode(double pixelX, double pixelY)
{
  rsDraggableNode* newNode = new rsDraggableNode(this, pixelX, pixelY);
  nodes.push_back(newNode);
  return newNode;
}
*/

int rsNodeEditor::addNode(double pixelX, double pixelY)
{
  rsDraggableNode* newNode = new rsDraggableNode(this, pixelX, pixelY);
  nodes.push_back(newNode);
  int i = size(nodes)-1;
  nodes[i]->index = i;
  return i;
}

void rsNodeEditor::removeNode(int i)
{
  delete nodes[i];  
  remove(nodes, i);
  for(i = i; i < size(nodes); i++)
    nodes[i]->decrementIndex();
  repaint();
}

void rsNodeEditor::removeNodeAt(int pixelX, int pixelY)
{
  int i = getNodeIndexAt(pixelX, pixelY);
  if(i != -1)
    removeNode(i);
}

int rsNodeEditor::moveNodeTo(int index, int pixelX, int pixelY)
{
  nodes[index]->setPixelPosition(pixelX, pixelY); // will indirectly trigger a repaint
  return index;
}

void rsNodeEditor::reIndexNode(int oldIndex, int newIndex)
{
  if(draggedNodeIndex == oldIndex)
    draggedNodeIndex = newIndex;
  while(oldIndex < newIndex)
  {
    //swapNodes(oldIndex, newIndex);
    swapNodes(oldIndex, oldIndex+1);
    oldIndex++;
  }
  while(oldIndex > newIndex)
  {
    //swapNodes(oldIndex, newIndex);
    swapNodes(oldIndex, oldIndex-1);
    oldIndex--;
  }
}

void rsNodeEditor::swapNodes(int i, int j)
{
  nodes[i]->index = j;
  nodes[j]->index = i;
  RAPT::rsSwap(nodes[i], nodes[j]);
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
  draggedNodeIndex = getNodeIndexAt(e.x, e.y);
  if(e.mods.isLeftButtonDown())
  {
    if(draggedNodeIndex == -1)
    {
      draggedNodeIndex = addNode((float)e.x, (float)e.y);
      repaint();
    }
  }
  else if(e.mods.isRightButtonDown())
  {
    removeNodeAt(e.x, e.y);
    draggedNodeIndex = -1;
    repaint();
  }
}

void rsNodeEditor::mouseDrag(const MouseEvent& e)
{
  if(draggedNodeIndex != -1)
    draggedNodeIndex = moveNodeTo(draggedNodeIndex, e.x, e.y);
}

void rsNodeEditor::mouseUp(const MouseEvent &e)
{
  draggedNodeIndex = -1;
}

void rsNodeEditor::mouseMove(const MouseEvent &e)
{
  rsDraggableNode* nodeUnderMouse = getNoteAt(e.x, e.y);
  if(nodeUnderMouse != nullptr)
    setMouseCursor(MouseCursor(MouseCursor::PointingHandCursor));
  else
    setMouseCursor(MouseCursor(MouseCursor::NormalCursor));
}

int rsNodeEditor::nodeChanged(int i)
{
  repaintOnMessageThread();
  return i;
}

// misc:

void rsNodeEditor::drawNodes(Graphics& g)
{
  bool drawNodeInfo = true;
  for(size_t i = 0; i < nodes.size(); i++)
  {
    float x, y;
    x = (float)nodes[i]->pixelX;
    y = (float)nodes[i]->pixelY;
    g.fillEllipse(x-0.5f*dotSize, y-0.5f*dotSize, dotSize, dotSize);
    if(drawNodeInfo)
    {
      String str = String(i) + ": x=" + String(x) + ", y=" + String(y);
      //String str = String(i) + "," + String(x) + "," + String(y);
      drawBitmapFontText(g, nodes[i]->pixelX, nodes[i]->pixelY-10, str, 
        &normalFont7px, getTextColour(), -1, Justification::centred);
    }
  }
}

//=================================================================================================

rsNodeBasedFunctionEditor::rsNodeBasedFunctionEditor(
  RAPT::rsInterpolatingFunction<double>* functionMapper)
{
  mapper = functionMapper;
}

double rsNodeBasedFunctionEditor::toPixelX(double modelX)
{
  return RAPT::rsLinToLin(modelX, xMin, xMax, 0.0, double(getWidth()-1));
}

double rsNodeBasedFunctionEditor::toPixelY(double modelY)
{
  return RAPT::rsLinToLin(modelY, yMin, yMax, double(getHeight()-1), 0.0);
}

double rsNodeBasedFunctionEditor::toModelX(double pixelX)
{
  return RAPT::rsLinToLin(pixelX, 0.0, double(getWidth()-1),  xMin, xMax);
}

double rsNodeBasedFunctionEditor::toModelY(double pixelY)
{
  return RAPT::rsLinToLin(pixelY, double(getHeight()-1), 0.0, yMin, yMax);
}

void rsNodeBasedFunctionEditor::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getHandleColour());
  for(int i = 0; i < getWidth()-1; i++)
  {
    double x1 = double(i);
    double x2 = double(i+1);
    double y1 = mapper->getValue(toModelX(x1));
    double y2 = mapper->getValue(toModelX(x2));
    g.drawLine((float)x1, (float)toPixelY(y1), (float)x2, (float)toPixelY(y2), 2.f);
  }
  /*
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
  */
  drawNodes(g);
}

int rsNodeBasedFunctionEditor::addNode(double x, double y)
{
  int i = (int) mapper->addDataPoint(toModelX(x), toModelY(y));
  rsDraggableNode* newNode = new rsDraggableNode(this, x, y);
  insert(nodes, newNode, i);
  nodes[i]->setIndex(i);
  for(int j = i+1; j < size(nodes); j++)
    nodes[j]->incrementIndex();
  return i;
}

void rsNodeBasedFunctionEditor::removeNode(int i)
{
  mapper->removeDataPoint(i);
  rsNodeEditor::removeNode(i);
}

int rsNodeBasedFunctionEditor::moveNodeTo(int index, int pixelX, int pixelY)
{
  int newIndex = (int)mapper->moveDataPoint(index, toModelX(pixelX), toModelY(pixelY));
  reIndexNode(index, newIndex);
  nodes[newIndex]->setPixelPosition(pixelX, pixelY, true); 
  repaint();
  return newIndex;
}

int rsNodeBasedFunctionEditor::nodeChanged(int nodeIndex)
{
  double x = nodes[nodeIndex]->getPixelX();
  double y = nodes[nodeIndex]->getPixelY();
  int newIndex = (int)mapper->moveDataPoint(nodeIndex, toModelX(x), toModelY(y));
  reIndexNode(nodeIndex, newIndex);
  rsNodeEditor::nodeChanged(nodeIndex);
  return newIndex;
}
