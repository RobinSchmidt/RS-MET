//#include "rojue_RTreeView.h"
//using namespace rojue;

RTreeViewNode::RTreeViewNode(const juce::String& _nodeText, int _identifier,
  const juce::String& _description, bool _isEnabled, bool _isTicked, bool _isOpen,
  RTreeViewNode *_parentNode)
{
  this->identifier  = _identifier;
  this->nodeText    = _nodeText;
  this->description = _description;
  this->isEnabled   = _isEnabled;
  this->isTicked    = _isTicked;
  this->isOpen      = _isOpen;
  this->parentNode  = _parentNode;
  this->userData    = NULL;
  deleteChildNodesOnDestruction = false;
}

RTreeViewNode::~RTreeViewNode()
{
  if( deleteChildNodesOnDestruction == true )
    deleteChildNodesRecursively();
}

void RTreeViewNode::setAllNodesUnticked()
{
  isTicked = false;
  for(int i = 0; i < size(childNodes); i++)
    childNodes[i]->setAllNodesUnticked();
}

/*
void RTreeViewNode::setAllNodesWithMatchingIdentifierTicked(int identifierToMatch)
{
  for(int i=0; i<childNodes.size(); i++)
  {
    if( childNodes[i]->getNodeIdentifier() == identifierToMatch )
      childNodes[i]->setTicked(true);
  }
  for(int i=0; i<childNodes.size(); i++)
    childNodes[i]->setAllNodesWithMatchingIdentifierTicked(identifierToMatch);
}
*/

void RTreeViewNode::addChildNode(RTreeViewNode *nodeToAdd)
{
  //childNodes.add(nodeToAdd);
  childNodes.push_back(nodeToAdd);
  nodeToAdd->setParentNode(this);
}

void RTreeViewNode::deleteChildNodesRecursively()
{
  for(int i = 0; i < size(childNodes); i++)
  {
    childNodes[i]->deleteChildNodesRecursively();
    delete childNodes[i];
  }
  childNodes.clear();
}

int RTreeViewNode::getNumLeafNodes() const
{
  if(this->isLeafNode())
    return 1;
  int leafs = 0;
  for(int i = 0; i < size(childNodes); i++)
    leafs += childNodes[i]->getNumLeafNodes();
  return leafs;
}

bool RTreeViewNode::hasAnyChildNodeChildNodes() const
{
  for(int i = 0; i < size(childNodes); i++)
  {
    if( childNodes[i]->hasChildNodes() )
      return true;
  }
  return false;
}

int RTreeViewNode::getLevel() const
{
  if( parentNode == NULL )
    return 0;
  else return parentNode->getLevel() + 1;
}

int RTreeViewNode::findIndexForIdentifier(int identifierToFind) const
{
  for(int i=0; i < size(childNodes); i++)
  {
    if( identifierToFind == childNodes[i]->identifier )
      return i;
  }
  return -1;
}

RTreeViewNode* RTreeViewNode::findNodeByData(void *dataPointerToFind)
{
  if( userData == dataPointerToFind )
    return this;
  else if( isLeafNode() )
    return NULL;
  else
  {
    RTreeViewNode *result;
    for(int i = 0; i < size(childNodes); i++)
    {
      result = childNodes[i]->findNodeByData(dataPointerToFind);
      if( result != NULL )
        return result;
    }
  }
  return NULL;
}

RTreeViewNode* RTreeViewNode::findNodeByIndentifier(int identifierToFind)
{
  if( identifier == identifierToFind )
    return this;
  else if( isLeafNode() )
    return NULL;
  else
  {
    RTreeViewNode *result;
    for(int i = 0; i < size(childNodes); i++)
    {
      result = childNodes[i]->findNodeByIndentifier(identifierToFind);
      if( result != NULL )
        return result;
    }
  }
  return NULL;
}

RTreeViewNode* RTreeViewNode::findNodeByText(const juce::String& textToFind)
{
  // the code duplication here really sucks - we have 3 times exactly the same logic - how can
  // this be avoided?
  if( nodeText == textToFind )
    return this;
  else if( isLeafNode() )
    return NULL;
  else
  {
    RTreeViewNode *result;
    for(int i = 0; i < size(childNodes); i++)
    {
      result = childNodes[i]->findNodeByText(textToFind);
      if( result != NULL )
        return result;
    }
  }
  return NULL;
}

RTreeViewNode* RTreeViewNode::getCopy()
{
  RTreeViewNode *copiedNode = new RTreeViewNode(nodeText, identifier, description, isEnabled,
    isTicked, isOpen, parentNode);
  for(int i = 0; i < size(childNodes); i++)
  {
    RTreeViewNode *copiedChild = childNodes[i]->getCopy();
    copiedNode->addChildNode(copiedChild);
  }
  return copiedNode;
}

//=================================================================================================
// class RTreeView:

RTreeView::RTreeView()
{
  rootNode                = NULL;
  drawRootNode            = true;
  openOrCloseNodesOnClick = true;
  xOffset                 = 0;
  yOffset                 = 0;
  setOpaque(true);

  addChildWidget( leftRightScrollBar = new RScrollBar(false) );
  leftRightScrollBar->setSingleStepSize(getNodeHeight());
  leftRightScrollBar->addListener(this);

  addChildWidget( upDownScrollBar = new RScrollBar(true) );
  upDownScrollBar->setSingleStepSize(getNodeHeight());
  upDownScrollBar->addListener(this);
}

RTreeView::~RTreeView()
{
  deleteAllChildren();
  dismissIfModal();
}

void RTreeView::setRootNode(RTreeViewNode *newRootNode)
{
  rootNode = newRootNode;
  repaint();  // shouldn't the be repaintOnMainThread or soemthing?
}

void RTreeView::setDrawRootNode(bool shouldBeDrawn)
{
  drawRootNode = shouldBeDrawn;
  repaint();
}

void RTreeView::registerTreeViewObserver(RTreeViewObserver* const observerToRegister)
{
  treeViewObservers.add(observerToRegister);
}

void RTreeView::deRegisterTreeViewObserver(RTreeViewObserver* const observerToDeRegister)
{
  treeViewObservers.removeFirstMatchingValue(observerToDeRegister);
}

int RTreeView::getRequiredWidth(bool ignoreOpenness) const
{
  /*
  if( drawRootNode == true )
    return getNodeWidth(rootNode) + 2 * (textMargin + outlineThickness);
  else
    return getNodeWidth(rootNode) + 2 * (textMargin + outlineThickness) - textMargin - plusMinusSize;
  */

  if( rootNode == NULL )
    return 0;

  int result = getNodeWidth(rootNode, ignoreOpenness) + 2 * (textMargin + outlineThickness);
  if( drawRootNode == false )
  {
    result -= (textMargin + plusMinusSize);
    if( !rootNode->hasAnyChildNodeChildNodes() )
      result -= (textMargin + plusMinusSize + 2);
      // the +2 is still somewhat arbitrary (it comes from the fact that our text-margin here is 2
      // pixels less than inside a combobox)
      // \todo: clean this up
  }
  return result;
}

int RTreeView::getRequiredHeight(bool ignoreOpenness) const
{
  if( rootNode == NULL )
    return 0;

  if( drawRootNode == true )
    return getNodeHeight(rootNode, ignoreOpenness) + 2 * (textMargin + outlineThickness);
  else
    return getNodeHeight(rootNode, ignoreOpenness) + 2 * (textMargin + outlineThickness) - getNodeHeight();
}

int RTreeView::getRootNodeX() const
{
  if( rootNode == NULL )
    return 0;

  int result = xOffset + outlineThickness + textMargin + plusMinusSize;

  if( drawRootNode == false && !rootNode->hasAnyChildNodeChildNodes() )
    result -= (plusMinusSize + textMargin + 2);
    // the +2 is still somewhat arbitrary (it comes from the fact that our text-margin here is 2
    // pixels less than inside a combobox)
    // \todo: clean this up

  return result;
}

int RTreeView::getRootNodeY() const
{
  return yOffset + outlineThickness + textMargin;
}

int RTreeView::getNodeWidth(const RTreeViewNode *node, bool ignoreOpenness) const
{
  int textWidth = font->getTextPixelWidth(node->nodeText, font->getDefaultKerning());
  int childWidth = 0;
  int maxChildWidth = 0;
  if( node->hasChildNodes() && (node->isOpen || ignoreOpenness) )
  {
    for(int i = 0; i < size(node->childNodes); i++)
    {
      childWidth = getNodeWidth(node->childNodes[i], ignoreOpenness);
      if( childWidth > maxChildWidth )
        maxChildWidth = childWidth;
    }
  }
  return jmax(textWidth+plusMinusSize+textMargin, maxChildWidth+plusMinusSize+textMargin);
}

int RTreeView::getNodeHeight() const
{
  return font->getFontHeight() + lineSpacing;
}

int RTreeView::getNodeHeight(const RTreeViewNode *node, bool ignoreOpenness) const
{
  int result = getNodeHeight();
  if( node->isOpen || ignoreOpenness )
  {
    for(int i = 0; i < size(node->childNodes); i++)
      result += getNodeHeight(node->childNodes[i], ignoreOpenness);
  }
  return result;
}

bool RTreeView::isPointClickable(int x, int y) const
{
  if( x < outlineThickness || y < outlineThickness+textMargin ||
    x > getWidth()-outlineThickness || y > getHeight()-outlineThickness )
    return false;

  if( upDownScrollBar->isVisible() && x > getWidth() - scrollBarThickness )
    return false;
  else if( leftRightScrollBar->isVisible() && y > getHeight() - scrollBarThickness )
    return false;
  else
    return true;
}

bool RTreeView::hasAnyNodeChildNodes() const
{
  if( rootNode == NULL )
    return false;

  for(int i=0; i < rootNode->getNumChildNodes(); i++)
  {
    if( rootNode->childNodes[i]->hasChildNodes() )
      return true;
  }
  return false;
}

bool RTreeView::isAnyNodeTicked() const
{
  if( rootNode == NULL )
    return false;

  for(int i=0; i < rootNode->getNumChildNodes(); i++)
  {
    if( rootNode->childNodes[i]->isTicked )
      return true;
  }
  return false;
}

void RTreeView::inputAttemptWhenModal()
{
  dismissIfModal();
}

void RTreeView::scrollBarMoved(RScrollBar* scrollBarThatHasMoved, const double newRangeStart)
{
  if( scrollBarThatHasMoved == upDownScrollBar )
    yOffset = -roundToInt(newRangeStart);
  else if( scrollBarThatHasMoved == leftRightScrollBar )
    xOffset = -roundToInt(newRangeStart);
  repaint();
}

void RTreeView::mouseEnter(const MouseEvent &e)
{
  RWidget::mouseEnter(e);
  repaint();
}

void RTreeView::mouseExit(const MouseEvent &e)
{
  RWidget::mouseExit(e);
  repaint();
}

void RTreeView::mouseDown(const MouseEvent &e)
{
  if( !isPointClickable(e.x, e.y) )
    return;

  RTreeViewNode *nodeUnderMouse = getNodeAtY(e.y);
  if( nodeUnderMouse != NULL )
    nodeClicked(nodeUnderMouse, e, getNodeClickPosition(nodeUnderMouse, e.x));

  dismissIfModal(); // move into ModalComponent
}

void RTreeView::mouseMove(const MouseEvent &e)
{
  repaint();
}

void RTreeView::mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel)
{
  if( upDownScrollBar->isVisible() )
    upDownScrollBar->moveScrollbarInSteps((int)-sign(wheel.deltaY));
}
//void RTreeView::mouseWheelMove(const MouseEvent &e, float wheelIncrementX, float wheelIncrementY)
//{
//  if( upDownScrollBar->isVisible() )
//    upDownScrollBar->moveScrollbarInSteps((int)-sign(wheelIncrementY));
//}

void RTreeView::paint(Graphics &g)
{
  g.fillAll(getBackgroundColour());
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), RWidget::outlineThickness);

  if( rootNode == nullptr )
    return;

  int x = getRootNodeX();
  int y = getRootNodeY();
  if( drawRootNode == true )
    y = drawNode(g, x - plusMinusSize, y, rootNode);
  else
  {
    for(int i=0; i<rootNode->getNumChildNodes(); i++)
      y = drawNode(g, x - plusMinusSize, y, rootNode->childNodes[i]);
  }

  // hide area that leaks through the scrollbar bottom-right corner:
  if( upDownScrollBar->isVisible() && leftRightScrollBar->isVisible() )
  {
    Rectangle<int> r(leftRightScrollBar->getRight()-outlineThickness,
      upDownScrollBar->getBottom()-outlineThickness, scrollBarThickness, scrollBarThickness);
    g.setColour(getBackgroundColour());
    g.fillRect(r);
    g.setColour(getOutlineColour());
    g.drawRect(r, outlineThickness);
  }
}

void RTreeView::resized()
{
  updateScrollBarBoundsAndVisibility();
}

void RTreeView::showModallyAt(int screenX, int screenY, int width, int height)
{
  if( width < 1 )
    width = getRequiredWidth(false);
  if( height < 1 )
    height = getRequiredHeight(false);

  setBounds(screenX, screenY, width, height);
  Desktop::getInstance().addGlobalMouseListener(this);
  addToDesktop(ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary );
  runModalLoop();
}

void RTreeView::showModallyAtMousePosition(int width, int height)
{
  Point<int> mousePosition = Desktop::getMousePosition();
  showModallyAt(mousePosition.getX(), mousePosition.getY(), width, height);
}

void RTreeView::dismissIfModal()
{
  if( isCurrentlyModal() )
  {
    Desktop::getInstance().removeGlobalMouseListener(this);
    removeFromDesktop();
    exitModalState(0);
  }
}

void RTreeView::updateScrollBarBoundsAndVisibility()
{
  // the logic to determine the available width/height and necessity for scrollbars are
  // interrelated in a somewhat messy way - maybe someday we should factor that stuff out into a
  // class Scrollable or a function getAvailableWidthAndHeight(int&, int&) or something:
  int availableHeight = getHeight();
  int requiredHeight  = getRequiredHeight(false);
  int availableWidth  = getWidth();
  int requiredWidth   = getRequiredWidth(false);

  bool upDownBarNeeded     = requiredHeight > availableHeight;
  bool leftRightBarNeeded  = requiredWidth  > availableWidth;
  bool leftRightBarNeeded2 = false;
  bool upDownBarNeeded2    = false;

  if( upDownBarNeeded ) {
    availableWidth      -= scrollBarThickness;
    leftRightBarNeeded2  = requiredWidth  > availableWidth;
  }
  if( leftRightBarNeeded ) {
    availableHeight  -= scrollBarThickness;
    upDownBarNeeded2  = requiredHeight > availableHeight;
  }

  if( (upDownBarNeeded == false) &&  (upDownBarNeeded2 == true) ) {
    upDownBarNeeded  = true;
    availableWidth  -= scrollBarThickness;
  }
  if( (leftRightBarNeeded == false) && (leftRightBarNeeded2 == true) ) {
    leftRightBarNeeded  = true;
    availableHeight    -= scrollBarThickness;
  }

  // OK, done with the mess - now set up the scrollbars:
  if( upDownBarNeeded ) {
    upDownScrollBar->setVisible(true);
    upDownScrollBar->setBounds(getWidth()-scrollBarThickness, 0, scrollBarThickness,
      availableHeight);
    upDownScrollBar->setRangeLimits(0.0, requiredHeight);
    upDownScrollBar->setCurrentRange(-yOffset, availableHeight);
  } else {
    yOffset = 0;
    upDownScrollBar->setVisible(false);
  }

  if( leftRightBarNeeded ) {
    leftRightScrollBar->setVisible(true);
    leftRightScrollBar->setBounds(0, getHeight()-scrollBarThickness, availableWidth,
      scrollBarThickness);
    leftRightScrollBar->setRangeLimits(0.0, requiredWidth);
    leftRightScrollBar->setCurrentRange(-xOffset, availableWidth);
  } else {
    xOffset = 0;
    leftRightScrollBar->setVisible(false);
  }

  if( upDownBarNeeded && leftRightBarNeeded ) {
    upDownScrollBar->setBounds(getWidth()-scrollBarThickness, 0, scrollBarThickness,
      availableHeight+outlineThickness);
    leftRightScrollBar->setBounds(0, getHeight()-scrollBarThickness,
      availableWidth+outlineThickness, scrollBarThickness);
  }
}

int RTreeView::drawNode(Graphics &g, int x, int y, const RTreeViewNode *nodeToDraw)
{
  // Avoid drawing nodes that are invisible because they are too high up or too low down:
  //if(y < -getNodeHeight() || y > getHeight() ) 
  //  return y + getNodeHeight(); 
  // This doesn't work correctly yet because it doesn't take into account the position of the 
  // scrollbar. When scrolling down, we see no nodes anymore. If determinin the visibility becoms 
  // more complex, maybe write a function isNodeVisible or isNodeWithinVisibleRange or something
  // like that (the name isNodeVisible may be confused with the visibility setting)


  // Highlight background for ticked nodes:
  if( nodeToDraw->isTicked )
  {
    g.setColour(getHandleColour());
    g.fillRect(outlineThickness, y, getWidth()-2*outlineThickness, getNodeHeight()-lineSpacing);
  }

  // Semi-highlight background for nodes where the mouse is over:
  Point<int> mousePosition = getMouseXYRelative();
  if( contains(mousePosition)
    && mousePosition.getY() >= y-lineSpacing/2
    && mousePosition.getY() <  y+getNodeHeight()-lineSpacing/2
    && nodeToDraw->isEnabled )
  {
    //g.setColour(getHandleColour().withMultipliedAlpha(0.325f));
    g.setColour(getHandleColour().withMultipliedAlpha(0.625f));
    g.fillRect(outlineThickness, y, getWidth()-2*outlineThickness, getNodeHeight()-lineSpacing);
  }

  // Draw the plus/minus button, if appropriate:
  if( nodeToDraw->hasChildNodes() )
  {
    float yOffset = 0.5f * (getNodeHeight() - plusMinusSize - lineSpacing);
    float yTmp    = y + yOffset;

    g.setColour(getTextColour());
    g.drawRect((float)x, (float)yTmp, (float)plusMinusSize, (float)plusMinusSize, 1.f);

    if( nodeToDraw->isOpen )
      drawBitmapFontText(g, x+2, (int) yTmp, "-", font, getTextColour());
    else
      drawBitmapFontText(g, x+2, (int) yTmp, "+", font, getTextColour());

    x += plusMinusSize + textMargin;
  }
  else
    x += plusMinusSize + textMargin;

  // Draw the node text:
  Colour textColour = getTextColour();
  if( !nodeToDraw->isEnabled )
    textColour = textColour.withMultipliedAlpha(0.625f);
  drawBitmapFontText(g, x, y, nodeToDraw->nodeText, font, textColour,
    font->getDefaultKerning(), Justification::topLeft);

  // Draw child nodes recursively:
  if( nodeToDraw->hasChildNodes() && nodeToDraw->isOpen )
  {
    y += getNodeHeight();
    for(int i = 0; i < nodeToDraw->getNumChildNodes(); ++i)
      y = drawNode(g, x, y, nodeToDraw->childNodes[i]);
    return y;
  }
  else
  {
    y += getNodeHeight();
    return y;
  }

  // ToDo:
  // -For larger TreeViews, the drawing performance is really bad. It becomes very unresponsive. 
  //  Maybe we should check x and y on entry and return early if they are such that we would draw
  //  outside the visible area. Maybe check, if JUCE has a TreeView and use that
  // -Maybe compare with juce::TreeView, see https://docs.juce.com/master/classTreeView.html
}

RTreeViewNode* RTreeView::getNodeAtY(int y)
{
  if( rootNode == NULL )
    return NULL;

  int yStart = getRootNodeY();

  if( drawRootNode == false )
    yStart -= getNodeHeight();

  return getNodeAtY(y, yStart, rootNode);
}

RTreeViewNode* RTreeView::getNodeAtY(int y, int &yStart, RTreeViewNode* nodeToStartWith)
{
  if( y >= yStart-lineSpacing/2 && y < yStart+getNodeHeight()-lineSpacing/2 )
    return nodeToStartWith;
  else if( nodeToStartWith->hasChildNodes() && nodeToStartWith->isOpen )
  {
    for(int i = 0; i < size(nodeToStartWith->childNodes); i++)
    {
      yStart += getNodeHeight(); // yStart is a reference, so it gets incremented in the recursion, too
      RTreeViewNode* currentNode = getNodeAtY(y, yStart, nodeToStartWith->childNodes[i]);
      if( currentNode != NULL )
        return currentNode;
    }
    return NULL;
  }
  else
    return NULL;
}

int RTreeView::getNodeClickPosition(RTreeViewNode *node, int pixelPositionX)
{
  int level      = node->getLevel();
  int xPlusMinus = outlineThickness + textMargin + level * (plusMinusSize + textMargin);
  if( !drawRootNode )
   xPlusMinus -= (plusMinusSize + textMargin);
  int xText      = xPlusMinus + (plusMinusSize + textMargin);
  if( pixelPositionX < xPlusMinus-textMargin )
    return LEFT_TO_PLUSMINUS;
  else if( pixelPositionX < xText )
    return ON_PLUSMINUS;
  else
    return ON_TEXT;
}

void RTreeView::nodeClicked(RTreeViewNode *nodeThatWasClicked, const MouseEvent &mouseEvent,
  int clickPosition)
{
  sendNodeClickNotification(nodeThatWasClicked, mouseEvent, clickPosition);
  if( openOrCloseNodesOnClick == true )
  {
    if( nodeThatWasClicked->hasChildNodes() )
      nodeThatWasClicked->isOpen = !nodeThatWasClicked->isOpen;
    updateScrollBarBoundsAndVisibility();
    repaint();
  }
}

void RTreeView::sendNodeClickNotification(RTreeViewNode *nodeThatWasClicked,
  const MouseEvent &mouseEvent, int clickPosition)
{
  for(int i=0; i<treeViewObservers.size(); i++)
    treeViewObservers[i]->treeNodeClicked(this, nodeThatWasClicked, mouseEvent, clickPosition);
}

void RTreeView::sendNodeChangeNotification(RTreeViewNode *nodeThatHasChanged)
{
  for(int i=0; i<treeViewObservers.size(); i++)
    treeViewObservers[i]->treeNodeChanged(this, nodeThatHasChanged);
}

//=================================================================================================
// class RTreeLeafNodeSelector:

RTreeLeafNodeSelector::RTreeLeafNodeSelector()
{
  selectedNode = NULL;
}

void RTreeLeafNodeSelector::deSelectNode()
{
  selectedNode = NULL;
}

void RTreeLeafNodeSelector::selectNodeByIndex(int indexToSelect, bool sendNotification)
{
  if( rootNode == NULL )
    return;

  rootNode->setAllNodesUnticked();
  if( indexToSelect >= 0  && indexToSelect < rootNode->getNumChildNodes() )
  {
    rootNode->childNodes[indexToSelect]->setTicked(true);
    selectedNode = rootNode->childNodes[indexToSelect];
  }
  else
    selectedNode = NULL;

  if( sendNotification == true )
    sendNodeChangeNotification(rootNode->childNodes[indexToSelect]);
}

void RTreeLeafNodeSelector::selectNodeByIdentifier(int nodeIdentifierToSelect,
  bool sendNotification)
{
  if( rootNode == NULL )
    return;

  rootNode->setAllNodesUnticked();
  //RTreeViewNode *selectedNode = rootNode->findNodeByIndentifier(nodeIdentifierToSelect); // old
  selectedNode = rootNode->findNodeByIndentifier(nodeIdentifierToSelect);                  // new
  if( selectedNode != NULL )
    selectedNode->setTicked(true);

  if( sendNotification == true )
    sendNodeChangeNotification(selectedNode);
}

void RTreeLeafNodeSelector::selectNodeByText(const juce::String& textToSelect,
  bool sendNotification)
{
  if( rootNode == NULL )
    return;

  rootNode->setAllNodesUnticked();
  //RTreeViewNode *selectedNode = rootNode->findNodeByText(textToSelect); // old
  selectedNode = rootNode->findNodeByText(textToSelect);                  // new
  if( selectedNode != NULL )
    selectedNode->setTicked(true);

  if( sendNotification == true )
    sendNodeChangeNotification(selectedNode);
}

void RTreeLeafNodeSelector::nodeClicked(RTreeViewNode *nodeThatWasClicked,
  const MouseEvent &mouseEvent, int clickPosition)
{
  if( rootNode == NULL )
    return;

  if( !nodeThatWasClicked->hasChildNodes() )
  {
    rootNode->setAllNodesUnticked();
    nodeThatWasClicked->setTicked(true);
    selectedNode = nodeThatWasClicked;
  }
  RTreeView::nodeClicked(nodeThatWasClicked, mouseEvent, clickPosition);
}
