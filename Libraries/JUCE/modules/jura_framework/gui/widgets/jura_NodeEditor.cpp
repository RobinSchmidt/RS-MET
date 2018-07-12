rsDraggableNode::rsDraggableNode(rsNodeEditor* editor, double _x, double _y)
{
  jassert(editor != nullptr);
  nodeEditor = editor;
  x = _x;
  y = _y;
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

void rsDraggableNode::setPosition(double newX, double newY, bool callNodeChanged)
{
  x = newX;
  y = newY;
  if(callNodeChanged)
    nodeEditor->nodeChanged(index);
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
  //xyMapper.setInputRange(0, 1, 0, 1);
  xyMapper.setInputRange(-1, 1, -1, 1);
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

int rsNodeEditor::addNode(double x, double y)
{
  xyMapper.unmap(&x, &y);  // map from pixel- to model coords
  rsDraggableNode* newNode = new rsDraggableNode(this, x, y);
  nodes.push_back(newNode);
  int i = size(nodes)-1;
  nodes[i]->setIndex(i);
  return constrainNode(i);
  //return i;
}

bool rsNodeEditor::removeNode(int i)
{
  if(!isNodeRemovable(i))
    return false; // maybe return false?
  delete nodes[i];  
  remove(nodes, i);
  for(i = i; i < size(nodes); i++)
    nodes[i]->decrementIndex();
  repaint();
  return true;
}

bool rsNodeEditor::removeNodeAt(int pixelX, int pixelY)
{
  int i = getNodeIndexAt(pixelX, pixelY);
  if(i != -1)
    return removeNode(i);
  return false;
}

int rsNodeEditor::moveNodeTo(int index, int pixelX, int pixelY)
{
  double x = xyMapper.unmapX(pixelX);
  double y = xyMapper.unmapY(pixelY);
  nodes[index]->setPosition(x, y); // will indirectly trigger a repaint
  return constrainNode(index);
  //return index;
}

void rsNodeEditor::reIndexNode(int oldIndex, int newIndex)
{
  if(draggedNodeIndex == oldIndex)
    draggedNodeIndex = newIndex;
  while(oldIndex < newIndex) {
    swapNodes(oldIndex, oldIndex+1); oldIndex++; }
  while(oldIndex > newIndex) {
    swapNodes(oldIndex, oldIndex-1); oldIndex--; }
}

void rsNodeEditor::swapNodes(int i, int j)
{
  nodes[i]->setIndex(j);
  nodes[j]->setIndex(i);
  RAPT::rsSwap(nodes[i], nodes[j]);
}

void rsNodeEditor::setDotSize(float newDotSize)
{
  dotSize = newDotSize;
  repaint();
}

// inquiry:

rsDraggableNode* rsNodeEditor::getNodeAt(int pixelX, int pixelY)
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
    float dx = x - getPixelX(nodes[i]);
    float dy = y - getPixelY(nodes[i]);
    float d2 = dx*dx + dy*dy;
    if(d2 <= r2)
      return (int)i;
  }
  return -1;
}

float rsNodeEditor::getPixelX(const rsDraggableNode* node)
{
  return (float) xyMapper.mapX(node->getX());
}

float rsNodeEditor::getPixelY(const rsDraggableNode* node)
{
  return (float) xyMapper.mapY(node->getY());
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

void rsNodeEditor::resized()
{
  xyMapper.setOutputRange(0, getWidth()-1, getHeight()-1, 0);
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
  rsDraggableNode* nodeUnderMouse = getNodeAt(e.x, e.y);
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

// observation:

void rsNodeEditor::sendNodeAddNotification(int nodeIndex)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->nodeWasAdded(this, nodeIndex);
}

void rsNodeEditor::sendNodeRemoveNotification(int nodeIndex)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->nodeWillBeRemoved(this, nodeIndex);
}

void rsNodeEditor::sendNodeMoveNotification(int nodeIndex)
{
  for(size_t i = 0; i < observers.size(); i++)
    observers[i]->nodeWasMoved(this, nodeIndex);
}

// misc:

void rsNodeEditor::drawNodes(Graphics& g)
{
  bool drawNodeInfo = true;
  for(size_t i = 0; i < nodes.size(); i++)
  {
    float x, y;
    x = getPixelX(nodes[i]);
    y = getPixelY(nodes[i]);

    juce::String xStr = String(x);
    juce::String yStr = String(y);

    g.fillEllipse(x-0.5f*dotSize, y-0.5f*dotSize, dotSize, dotSize);
    if(drawNodeInfo)
    {
      String str = String(i) + ": x=" + xStr + ", y=" + yStr;  // verbose
      //String str = String(i) + "," + String(x) + "," + String(y); // compact



      drawBitmapFontText(g, (int)getPixelX(nodes[i]), (int)getPixelY(nodes[i])-10, str, 
        &normalFont7px, getTextColour(), -1, Justification::centred);
    }
  }
}

//=================================================================================================

rsNodeBasedFunctionEditor::rsNodeBasedFunctionEditor(
  RAPT::rsNodeBasedFunction<double>* functionMapper, CriticalSection* lockToUse)
{
  lock = lockToUse;
  valueMapper = functionMapper;
}

void rsNodeBasedFunctionEditor::setFunctionToEdit(RAPT::rsNodeBasedFunction<double>* func) 
{ 
  valueMapper = func;
  updateDraggableNodesArray();
}

void rsNodeBasedFunctionEditor::paint(Graphics& g)
{
  ScopedPointerLock spl(lock);
  g.fillAll(getBackgroundColour());
  g.setColour(getHandleColour());
  for(int i = 0; i < getWidth()-1; i++)
  {
    double x1 = double(i);
    double x2 = double(i+1);
    double y1 = valueMapper->getValue(xyMapper.unmapX(x1));
    double y2 = valueMapper->getValue(xyMapper.unmapX(x2));
    g.drawLine((float)x1, (float) xyMapper.mapY(y1), (float)x2, (float) xyMapper.mapY(y2), 2.f);
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
  ScopedPointerLock spl(lock);
  xyMapper.unmap(&x, &y);
  clipIfDesired(&x, &y);
  int i = (int) valueMapper->addNode(x, y);
  rsDraggableNode* newNode = new rsDraggableNode(this, x, y);
  insert(nodes, newNode, i);
  nodes[i]->setIndex(i);
  for(int j = i+1; j < size(nodes); j++)
    nodes[j]->incrementIndex();
  return i;
}

bool rsNodeBasedFunctionEditor::removeNode(int i)
{
  ScopedPointerLock spl(lock);
  if(rsNodeEditor::removeNode(i))
  {
    return valueMapper->removeNode(i);
  }
  return false;
}

int rsNodeBasedFunctionEditor::moveNodeTo(int index, int pixelX, int pixelY)
{
  ScopedPointerLock spl(lock);
  double x = xyMapper.unmapX(pixelX);
  double y = xyMapper.unmapY(pixelY);
  clipIfDesired(&x, &y);
  int newIndex = (int)valueMapper->moveNode(index, x, y);
  reIndexNode(index, newIndex);
  //nodes[newIndex]->setPosition(x, y, true); 
  nodes[newIndex]->setPosition(
    valueMapper->getNodeX(newIndex), valueMapper->getNodeY(newIndex), true);
  repaint();
  return newIndex;
}

int rsNodeBasedFunctionEditor::nodeChanged(int nodeIndex)
{
  ScopedPointerLock spl(lock);
  double x = getPixelX(nodes[nodeIndex]);
  double y = getPixelY(nodes[nodeIndex]);
  xyMapper.unmap(&x, &y);
  //clipIfDesired(&x, &y); // needed?
  int newIndex = (int)valueMapper->moveNode(nodeIndex, x, y);
  reIndexNode(nodeIndex, newIndex);
  rsNodeEditor::nodeChanged(nodeIndex);
  return newIndex;
}

void rsNodeBasedFunctionEditor::updateDraggableNodesArray()
{
  nodes.clear();
  if(valueMapper == nullptr)
    return;
  const std::vector<RAPT::rsFunctionNode<double>> funcNodes = valueMapper->getNodes();
  for(size_t i = 0; i < funcNodes.size(); i++) 
  {
    rsDraggableNode* node = new rsDraggableNode(this, funcNodes[i].getX(), funcNodes[i].getY());
    node->setIndex((int)i);
    nodes.push_back(node); 
  }
}

void rsNodeBasedFunctionEditor::clipIfDesired(double* x, double* y)
{
  if(clipRanges) {
    *x = clip(*x, xyMapper.getInMinX(), xyMapper.getInMaxX());
    *y = clip(*y, xyMapper.getInMinY(), xyMapper.getInMaxY()); }
}