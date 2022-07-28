#ifndef jura_RTreeView_h
#define jura_RTreeView_h

/** This class is used by RTreeView to represent the nodes that appear in the menu. 

\todo: allow not only for text-nodes but also for Component nodes

\todo: maybe factor out a general TreeNode class that can be used for anything that is organized in
a tree. Maybe a NodeData class must be factored out ...

*/

class JUCE_API RTreeViewNode
{

  friend class RTreeView;
  friend class RTreeLeafNodeSelector;
  friend class RPopUpMenu;  // it is not very clean to have this as friend - maybe change this later

public:

  /** Constructor. */
  RTreeViewNode(const juce::String& nodeText = juce::String(), int identifier = -1,
    const juce::String& description = juce::String(), bool isEnabled = true, 
    bool isTicked = false, bool isOpen = true, RTreeViewNode *parentNode = NULL);

  /** Destructor. */
  virtual ~RTreeViewNode();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets a new text for this node. */
  virtual void setNodeText(const juce::String& newText) { nodeText = newText; }

  /** Enables or disables this node. */
  virtual void setEnabled(bool shouldBeEnabled) { isEnabled = shouldBeEnabled; }

  /** Sets this node as ticked (or not). */
  virtual void setTicked(bool shouldBeTicked) { isTicked = shouldBeTicked; }

  /** Sets this node open (or not). Relevant only for non-leaf nodes. */
  virtual void setOpen(bool shouldBeOpen) { isOpen = shouldBeOpen; }

  /** Associates arbitrary data with this node - its up to the user to interpret this data. When 
  the node is deleted, the data pointed to will not be freed here (\todo make this optionally 
  possible). */
  virtual void setUserData(void *newDataPointer) { userData = newDataPointer; }

  /** Sets this and all the child-nodes into unticked state. */
  virtual void setAllNodesUnticked();

  /** Sets all nodes into ticked state that have an identifier that matches the passed 
  "identifierToMatch" - works recursively on the child-nodes. */
  //virtual void setAllNodesWithMatchingIdentifierTicked(int identifierToMatch);

  /** Selects, whether or not the child nodes will be deleted (recursively) on destruction. By 
  default, they are not. If you set this to true for the root node of your tree, the destruction of
  the root node will clean up the memory for the whole tree. */
  virtual void setDeleteChildNodesOnDestruction(bool shouldDelete) 
  { 
    deleteChildNodesOnDestruction = shouldDelete; 
  }

  /** Sets up the parent node for this node. */
  virtual void setParentNode(RTreeViewNode *newParentNode) { parentNode = newParentNode; }

  /** Adds an node to our array of sub-nodes. */
  virtual void addChildNode(RTreeViewNode *nodeToAdd);

  /** Deletes the child-nodes recursively. Used internally in the destructor to clean up memory 
  when the deleteChildNodesOnDestruction flag is true. */
  virtual void deleteChildNodesRecursively();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of (direct) child nodes. */
  virtual int getNumChildNodes() const { return (int) childNodes.size(); }

  /** Returns the number of leaf nodes in this node. That number is 1, if the node is itself a leaf 
  node (i.e. has no child nodes), otherwise it will recursively call getNumLeafNodes on all child 
  nodes and add them up. */
  virtual int getNumLeafNodes() const;

  /** Returns true when this node has one or more child-nodes, false otherwise. */
  virtual bool hasChildNodes() const { return childNodes.size() > 0; }

  /** Returns true when this node is a leaf node (i.e. does not have any child-nodes), false 
  otherwise. */
  virtual bool isLeafNode() const { return !hasChildNodes(); }

  /** Returns true when this node is a root node (i.e. does not have a parent node), false 
  otherwise. */
  virtual bool isRootNode() const { return parentNode == NULL; }

  /** Returns true if any child-node of this node has itself one or more child-nodes, false 
  otherwise. False is also returned when this node does not have any child-nodes. */
  virtual bool hasAnyChildNodeChildNodes() const;

  /** Returns the level of the node in the tree - root nodes have level 0, their direct child-nodes
  level 1 and so on. */
  virtual int getLevel() const;

  /** Returns the index in the child-nodes array that has the given identifier. If no child-node 
  has this identifier, it will return -1, if more than one child nodes have the given identifier, 
  the first index will be returned. This function only checks the direct child-nodes, it does not 
  search recursively. */
  virtual int findIndexForIdentifier(int identifierToFind) const;

  /** Returns the text for this node. */
  virtual juce::String getNodeText() const { return nodeText; }

  /** Returns the identifier for this node. */
  virtual int getNodeIdentifier() const { return identifier; }

  /** Returns the data that was previously associated with this node (via setUserData), if any, 
  NULL otherwise. Its up to the user to interpret this data. */
  virtual void* getUserData() const { return userData; }

  /** Searches inside this node and its descendants for a particular pointer to user-data. If a 
  node containing this data-pointer is found, this function will return a pointer to the 
  respective node. If none is found, it returns a NULL pointer. If several nodes with the 
  data-pointer exist, it returns the first one that is found. */
  RTreeViewNode* findNodeByData(void *dataPointerToFind);

  /** Similar to findNodeByData, but searches for a particular identifier. */
  RTreeViewNode* findNodeByIndentifier(int identifierToFind);

  /** Similar to findNodeByData, but searches for a particular text. */
  RTreeViewNode* findNodeByText(const juce::String& textToFind);

  /** Returns the text of this node. */
  const juce::String& getText() const { return nodeText; }

  /** From some given node, it returns the next leaf node in the tree. The second argument decides 
  whether we consider only sibling-nodes which are those nodes that are on the same level and 
  branch of the tree. If false is passed, the function will possibly crawl up and/or down the 
  branches in order to find the next leaf node. This is useful for skipping through nodes one after 
  another. */
  //virtual RTreeViewNode* getNextLeafNode(RTreeViewNode *node, bool considerOnlySiblings);
   // \todo: implement RTreeViewNode::getNextLeafNode

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Creates a deep copy of this node and returns a pointer to it. The caller is responsible for 
  deleting it eventually. */
  virtual RTreeViewNode* getCopy();
  // maybe rename to createDeepCopy


protected:

  // node-data:
  juce::String nodeText, description;
  int  identifier;  //, level;
  bool isEnabled, isTicked, isOpen;  // , isHidden  - maybe remove the "is" prefix, reserve it for functions like isTicked()
  void *userData; // arbitrary data that is associated with this node
  // bool userDataOwned = false;


  // data for maniging the tree-structure:
  RTreeViewNode               *parentNode;
  std::vector<RTreeViewNode*> childNodes;
  bool deleteChildNodesOnDestruction;

private:

  RTreeViewNode(const RTreeViewNode& other);
  RTreeViewNode& operator=(const RTreeViewNode& other);

  juce_UseDebuggingNewOperator;
};


//=================================================================================================
// class RTreeView and asscoiated observer class:

class RTreeView;

/** Observer class for RTreeViews */

class JUCE_API RTreeViewObserver
{

public:

  virtual ~RTreeViewObserver() {}

  /** Callback that gets called when the user has clicked a node in an RTreeView. The clickPosition 
  argument is one of the values defined in the RTreeView::nodeClickPositions enumeration and may be 
  used to trigger different actions depending on where exactly the node was clicked. */
  virtual void treeNodeClicked(RTreeView *treeView, RTreeViewNode *nodeThatWasClicked, 
    const MouseEvent &mouseEvent, int clickPosition) = 0;

  /** Callback that gets called when a node in the tree was changed somehow. */
  virtual void treeNodeChanged(RTreeView *treeView, RTreeViewNode *nodeThatHasChanged) = 0;

};

//=================================================================================================

/**

This class implements a tree view that can be used to display hierarchical structures of choices. 
The class also provides means to show up modally on the desktop (for example, it is used by the 
RComboBox class that way).

maybe let the user choose between clicking behaviors:

\todo factor out a class Scrollable (or ScrollableComponentContainer): getFullWidth, 
getSectionWidth, same for height...maybe use two Rectangles to implement it: 
(0, 0, fullWidth, fullHeight), (x, y, sectionWidth, sectionHeight)

\todo maybe make a template Observer and Observee class to get rid of the boilerplate 
de/registerObserver code

*/

class JUCE_API RTreeView : public RWidget, public RScrollBarListener
{

  friend class RPopUpMenu;

public:

  /** An enumeration of the distinguished x-positions where the user can click. These values used 
  as additional info in the node-click callback functions. */
  enum nodeClickPositions
  {
    LEFT_TO_PLUSMINUS,
    ON_PLUSMINUS,
    ON_TEXT
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates an empty tree view. */
  RTreeView();

  /** Destructor. */
  virtual ~RTreeView();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the root node for the tree to be shown. The RTreeView will NOT delete the root-node in 
  its destructor - that is: it does not take over ownership over the node. If the node is not 
  needed anymore, you may call setRootNode with a NULL pointer as argument and delete the node 
  yourself. */
  virtual void setRootNode(RTreeViewNode *newRootNode);

  /** Selects whether or not the root node should be drawn. */
  virtual void setDrawRootNode(bool shouldBeDrawn);

  /** Selects, whether or not a click on a node automatically triggers opening/closing of the node 
  (and subsequent repainting). */
  virtual void setOpenOrCloseNodesOnClick(bool shouldOpenOrClose) 
  { 
    openOrCloseNodesOnClick = shouldOpenOrClose; 
  }

  /** Registers an observer that will be called when a node is clicked. */
  virtual void registerTreeViewObserver(RTreeViewObserver* const observerToRegister);

  /** Deregisters a previously-registered observer. */
  virtual void deRegisterTreeViewObserver(RTreeViewObserver* const observerToDeRegister);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns a pointer to the selected node, if any, otherwise NULL will be returned. The returned 
  pointer is a direct reference to our data here, so the caller is not supposed to delete it. */
  //virtual RTreeViewNode* getSelectedNode() const;

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void inputAttemptWhenModal();
  //virtual bool canModalEventBeSentToComponent(const Component* targetComponent);

  virtual void scrollBarMoved(RScrollBar* scrollBarThatHasMoved, const double newRangeStart);

  virtual void mouseEnter(const MouseEvent &e);  // overriden to make sure that highlighting appears, if appropriate
  virtual void mouseExit(const MouseEvent &e);   // overriden to make sure that highlighting disappears
  virtual void mouseDown(const MouseEvent &e);   // overriden to update the selected node
  virtual void mouseMove(const MouseEvent &e);   // overriden to update the highlighting
  virtual void mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel);
  //virtual void mouseWheelMove(const MouseEvent &e, float wheelIncrementX, float wheelIncrementY);


  virtual void paint(Graphics&);
  virtual void resized();

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Shows the menu modally at a specified location on the screen. The passed coordinates are used 
  to position the top-left pixel of the menu. You may also pass a width for the menu - if you leave 
  that at the default value (zero), the menu will work out the required width itself (by using 
  getRequiredWidth()). The same goes likewise for the height. When the function has returned, the 
  user can retrieve the selected node via getSelectedNode(). */
  virtual void showModallyAt(int screenX, int screenY, int width = 0, int height = 0);

  /** Similar to showAt, but uses the current position of the mouse for the top-left corner. */
  virtual void showModallyAtMousePosition(int width = 0, int height = 0);

  // \todo: get rid of these showModally functions - instead, write a class ModalComponent with a function setContentComponent that 
  // managaes all the aspects of modality



protected:

  /** Returns the required width such that all nodes fit horizontally into the menu. The 
  ignoreOpenness parameter determines whether or not only open nodes shall count. */
  virtual int getRequiredWidth(bool ignoreOpenness) const;

  /** Returns the required width such that all nodes fit vertically into the menu. The 
  ignoreOpenness parameter determines whether or not only open nodes shall count. */
  virtual int getRequiredHeight(bool ignoreOpenness) const;

  /** Returns the x-coordinate for the root node taking into account the outline, and 
  scroll-position. */
  virtual int getRootNodeX() const;

  /** Returns the y-coordinate for the root node taking into account the outline, and 
  scroll-position. */
  virtual int getRootNodeY() const;

  /** Returns the width for a specific node, taking into account child-nodes, and their children, 
  etc. (recursively). When the ignoreOpenness parameter is false, only open nodes count. */
  virtual int getNodeWidth(const RTreeViewNode *node, bool ignoreOpenness) const;

  /** Returns the height for one generic node without taking into account child-itmes. */
  virtual int getNodeHeight() const;

  /** Returns the height for a specific node, taking into account child-nodes, and their children, 
  etc. (recursively). When the ignoreOpenness parameter is false, only open nodes count. */
  virtual int getNodeHeight(const RTreeViewNode *node, bool ignoreOpenness) const;

  /** Returns true when the point given by coordinates x and y is considered clickable. It might 
  not be due to being behind a scrollbar or inside the bottom right area where the scrollbars meet 
  when they are both visible. */
  virtual bool isPointClickable(int x, int y) const;

  /** Returns true when at least one of the nodes has sub-nodes (in which case we need to display 
  an arrow and possible a corresponding sub-menu. */
  virtual bool hasAnyNodeChildNodes() const;

  /** Returns true when at least one of the nodes is ticked (in which case we need to display a 
  tick-mark). */
  virtual bool isAnyNodeTicked() const;

  /** Removes the menu from the desktop and exits the modal loop. */
  virtual void dismissIfModal();

  /** Updates the visibility and bounds of the horizontal and vertical scrollbars depending on 
  whether either of the two or both are actually needed. */
  virtual void updateScrollBarBoundsAndVisibility();

  /** Draws an node at the given position. If the node has open child nodes, it will draw them too 
  (recursively). The return value is the appropriate y-coordinate for the next node to draw (we 
  use it that way, because the appropriate y will depend on the number of child-nodes that will be 
  drawn in this call). */
  virtual int drawNode(Graphics &g, int x, int y, const RTreeViewNode *nodeToDraw);
    // todo: remove the return value and pass x and y as references

  /** Returns a pointer to the node that is at the given y-coordinate (or NULL, if none) .*/
  virtual RTreeViewNode* getNodeAtY(int y);

  /** Used internally by getNodeAtY(int) for recursion. */
  virtual RTreeViewNode* getNodeAtY(int y, int &yStart, RTreeViewNode* nodeToStartWith);

  /** Takes a node and an x-coodinate in pixels and returns one of the values defined in the 
  nodeClickPositions enumeration according to whether the x-coordinate is to the left, on or to the 
  right of the plus/minus button. The plus/minus button has a bit of tolerance here: clicks that 
  are inside the margin between the plusminus and the text are considered to be on the plusminus, 
  too. Likewise, clicks that are within this this margin-distance to the left of the button are 
  also considered as on the button. */
  virtual int getNodeClickPosition(RTreeViewNode *node, int pixelPositionX);

  /** Internal callback that gets called when a node was clicked - intended to be overriden in 
  subclasses in order to do stuff when a node was clicked (like opening/closing, ticking, editing, 
  etc.). The clickPosition argument is one of the values defined in the nodeClickPositions 
  enumeration and may be used to trigger different actions depending on where exactly the node was 
  clicked. The baseclass implementation notifies our listeners and opens or closes the node (unless 
  it's a leaf node). */
  virtual void nodeClicked(RTreeViewNode *nodeThatWasClicked, const MouseEvent &mouseEvent, 
    int clickPosition);

  /** Calls treeNodeClicked on all of our observers. */
  virtual void sendNodeClickNotification(RTreeViewNode *nodeThatWasClicked, 
    const MouseEvent &mouseEvent, int clickPosition);

  /** Calls treeNodeChanged on all of our observers. */
  virtual void sendNodeChangeNotification(RTreeViewNode *nodeThatHasChanged);


  // data members:

  static const int textMargin         = 4;   // margin in pixels between the text and the outline
  static const int lineSpacing        = 2;   // vertical spacing in pixels between two nodes
  static const int plusMinusSize      = 10;  // size for the plus or minus that opens/closes the node
  static const int scrollBarThickness = 16;  // thickness of the scrollbars

  RTreeViewNode *rootNode; //, *selectedNode;    // factor the selected node stuff into a subclass RTreeNodeSelector

  RScrollBar *leftRightScrollBar, *upDownScrollBar;

  bool drawRootNode;  // this is not really used yet
  bool openOrCloseNodesOnClick;

  int xOffset, yOffset;  // offset for the root node in pixels for both coordinates

  juce::Array<RTreeViewObserver*> treeViewObservers;


  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// class RTreeLeafNodeSelector:

/** This class is a subclass of RTreeView that allows for selecting one of the leaf-nodes. This is
useful for popup-menus that are used by comboboxes. */

class JUCE_API RTreeLeafNodeSelector : public RTreeView
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Creates an empty RTreeLeafNodeSelector - the pointer to the selected node is initialized to 
  NULL. */
  RTreeLeafNodeSelector();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Selects the (direct) child node of this node that has the given index. The second argument 
  decides whether a node-change notification should be sent out to our observers. */
  virtual void selectNodeByIndex(int indexToSelect, bool sendNodeChangeNotification);
    // maybe delete - doesn't seem to be useful

  /** Selects the 1st leaf node that has an identifier that matches the "nodeIdentifierToSelect" 
  parameter. */
  virtual void selectNodeByIdentifier(int nodeIdentifierToSelect, bool sendNodeChangeNotification);

  /** Selects the 1st leaf node that has an identifier that matches the "textToSelect" 
  parameter. */
  virtual void selectNodeByText(const juce::String& textToSelect, bool sendNodeChangeNotification);

  /** Deselects any node that might be currently selected. If no node is selected, it does 
  nothing. */
  virtual void deSelectNode();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns a pointer to the selected node, if any, otherwise NULL will be returned. The returned 
  pointer is a direct reference to our data here, so the caller is not supposed to delete it. */
  virtual RTreeViewNode* getSelectedNode() const { return selectedNode; };

protected:

  /** Overriden to implement the (exclusive) selection of a leaf-node. */
  virtual void nodeClicked(RTreeViewNode *nodeThatWasClicked, const MouseEvent &mouseEvent, 
    int clickPosition);

  RTreeViewNode *selectedNode;    // pointer to our currently selected node

};

#endif  
