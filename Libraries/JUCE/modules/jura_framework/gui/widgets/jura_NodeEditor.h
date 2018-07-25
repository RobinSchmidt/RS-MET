#ifndef jura_NodeEditor_h
#define jura_NodeEditor_h  

class rsNodeEditor;

/** A class for nodes that can be dragged around in a node-editor and my have Parameter objects
associated with the x- and y coordinates. */

class rsDraggableNode : public ParameterObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  /** Constructor */
  rsDraggableNode(rsNodeEditor* editorContainingThisNode, double x, double y);

  /** Destructor */
  ~rsDraggableNode();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Assigns the parameter object associated with the x-coordinate of this node. */
  virtual void assignParameterX(Parameter* newParameterX);

  /** Assigns the parameter object associated with the y-coordinate of this node. */
  virtual void assignParameterY(Parameter* newParameterY);

  /** With this function, you can associate additional parameters with each node that controls 
  additional features of a node such as a bandwidth for an eq-band, a slope, a shape, whatever. */
  virtual void addNodeParameter(Parameter* newParameter);

  /** Clears the array of additional node feature parameters and optionally deletes the parameter 
  objects. */
  virtual void clearNodeParameters(bool deleteObjects);

  /** Sets the position of this node (in model coordinates) calls the nodeChanged function of 
  our editor. */
  virtual void setPosition(double newX, double newY, bool callNodeChanged = true);

  /** This function is used to set up the "index" member variable which should always reflect the
  index in the "nodes" array in the editor. You need to take care of this in subclasses when 
  adding, dragging and removing nodes. */
  inline void setIndex(int newIndex) { index = newIndex; }

  /** Increments the index of this node by 1. */
  inline void incrementIndex() { index++; }

  /** Decrements the index of this node by 1. */
  inline void decrementIndex() { index--; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the x-coordinate of this node. */
  inline double getX() const { return x; }

  /** Returns the y-coordinate of this node. */
  inline double getY() const { return y; }

  /** Returns a pointer to the parameter object that is associated with the node's x-variable 
  (maybe a nullptr, if none is assigned). */
  Parameter* getParameterX() { return paramX; }

  /** @see getParameterX */
  Parameter* getParameterY() { return paramY; }

  /** Returns a pointer to the node feature parameter with given index, i.e. one of those added via
  addNodeParameter. If the index is out of range, this will lead to an access violation. */
  Parameter* getNodeParameter(int index) { return nodeParams[index]; }

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  /** Overriden in order to cause an update/repaint in our editor. */
  virtual void parameterChanged(Parameter* p) override;

protected:
   
  int index = -1;
  double x, y;                                    // (model) coordinates of the node 
  Parameter *paramX = nullptr, *paramY = nullptr; // parameter objects associated with x and y
  std::vector<Parameter*> nodeParams;             // additional node feature parameters
  rsNodeEditor* nodeEditor;                       // editor which edits this node 
                                                  // todo: maybe allow more than one editor

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsDraggableNode)
};

//=================================================================================================

/** Baseclass for classes that must keep track of the state of an rsNodeEditor object. */

class JUCE_API rsNodeEditorObserver
{

public:

  virtual ~rsNodeEditorObserver() {}

  /** Callback that will be called after a new node was added. */
  virtual void nodeWasAdded(rsNodeEditor* editor, int nodeIndex) = 0;

  /** Callback that will be called before an existing node will be removed. */
  virtual void nodeWillBeRemoved(rsNodeEditor* editor, int nodeIndex) = 0;

  /** Called when a node was selected. The passed value may be -1, indicating that no node was 
  selected, i.e. a previously selcetd node was de-selected. */
  virtual void nodeWasSelected(rsNodeEditor* editor, int nodeIndex) = 0;

  /** Callback that will be called after an existing node was moved. Your subclass may choose to
  not override this, if it is already an observer of the parameters associated with the node's x/y
  coordinates. */
  virtual void nodeWasMoved(rsNodeEditor* editor, int nodeIndex) {};

  // virtual void nodeParameterChanged(rsNodeEditor* editor, int nodeIndex, int paramIndex) {}
  // ...maybe

};

//=================================================================================================

/** This is a baseclass for editors that need to create/remove/drag-around a number of nodes. */

class JUCE_API rsNodeEditor : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction
  
  rsNodeEditor();
  virtual ~rsNodeEditor();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adds a new node at the given pixel position and returns the index at which it was inserted 
  into our node array. The new node will use nullptrs for the Parameter objects associated with the 
  x- and y-coordinate. Subclasses may want to override this in order to create and assign actual 
  Parameter objects with associated callbacks. */
  virtual int addNode(double pixelX, double pixelY);

  /** Removes the node with given index and returns if this was successful. */
  virtual bool removeNode(int index);

  /** If there is a node at the given pixel position, this function will remove it (otherwise it 
  will have no effect). */
  virtual bool removeNodeAt(int pixelX, int pixelY);

  /** Moves the node with given index to the given new pixel coordinates. The return value is 
  the new index of the node. In the baseclass implementation, that index is the same as the input 
  parameter "index" but subclasses may want to re-order nodes depending their positions. */
  virtual int moveNodeTo(int index, int pixelX, int pixelY);

  /** Moves the node that is currently at oldIndex to the newIndex in our array of nodes. */
  virtual void reIndexNode(int oldIndex, int newIndex);

  /** Swaps the two nodes at the given indices in out node-array. */
  virtual void swapNodes(int index1, int index2);

  /** Sets the size of the dots in pixels. */
  void setDotSize(float newDotSize);

  /** Sets the range of values between which nodes can be moved (in model coordinates). */
  void setValueRange(double minX, double maxX, double minY, double maxY)
  { xyMapper.setInputRange(minX, maxX, minY, maxY); }

  /** Chooses one of the nodes to be the selected node. Useful in order to attach other widgets
  with the node to edit the node settings more precisely and get access to advanced node 
  paraeters, if any (for example: the curve setting in a node-based function). To deselect the
  node, you may pass -1. */
  virtual void selectNode(int index);

  /** Deselects the selected node (if any). */
  inline void deSelectNode() { selectNode(-1); }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns a pointer to a node at the given position, if any. If there's no node at the given
  position, it returns a nullptr. If there are several nodes at the given position, it will return
  the first one that is found. 
  todo: if there are several nodes, return the first one among those that are closest to the given 
  pixel position */
  rsDraggableNode* getNodeAt(int pixelX, int pixelY);
    
  /** Returns the index in our array of nodes a node at the given position. If there is none, it 
  will return -1. */
  int getNodeIndexAt(int pixelX, int pixelY);

  /** Returns a pointer to the node with given index (or a nullptr if the index is out of 
  range). */
  rsDraggableNode* getNode(int nodeIndex);

  /** Returns the index of the given node in our array nodes. Returns -1, if the node isn't 
  found. */
  int getNodeIndex(rsDraggableNode* node);

  /** Returns the x-coordinate in pixels of the given node. */
  float getPixelX(const rsDraggableNode* node);

  /** Returns the y-coordinate in pixels of the given node. */
  float getPixelY(const rsDraggableNode* node);



  //-----------------------------------------------------------------------------------------------
  // \name Hooks

  /** This function will be called before an attempt to remove a node and will not remove it, if 
  that function returns false. The baseclass implementation just returns true but you can override 
  it in a subclass if you subclass - for example - requires a certain minimum number of nodes. */
  virtual bool isNodeRemovable(int index) { return true; }

  /** This function is called after a node has been inserted or moved to a new position. You can
  override it, if you need to apply costraints on the positions of nodes, like for example that the
  coordinates of the first and/or last node must hae certain values. Because changing the position 
  of a node may change its index, too, this function should return the new index of the after the
  constraints have been applied. The baseclass implementation does nothing and just returns the 
  same index that you give to it. */
  virtual int constrainNode(int index) { return index; }

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void paint(Graphics& g) override;
  virtual void resized() override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;
  virtual void mouseUp(const MouseEvent &e) override;
  virtual void mouseMove(const MouseEvent &e) override;
  //virtual void nodeChanged(const rsDraggableNode* node);

  /** A callback that will be called when the node with given index has changed. If you override 
  this you should return the new index of the node (you may want to re-order nodes depending on 
  their new positions and the baseclass needs to keep track of that). */
  virtual int nodeChanged(int nodeIndex);

  //-----------------------------------------------------------------------------------------------
  // \name Oberservation

  void registerObserver(rsNodeEditorObserver* obs) { appendIfNotAlreadyThere(observers, obs); }

  void deRegisterObserver(rsNodeEditorObserver* obs) { removeFirstOccurrence(observers, obs); }

  void sendNodeAddNotification(   int nodeIndex);
  void sendNodeRemoveNotification(int nodeIndex);
  void sendNodeMoveNotification(  int nodeIndex);
  void sendNodeSelectNotification(int nodeIndex);

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Darws all the nodes. You may want to override it if you need special drawing. */
  virtual void drawNodes(Graphics& g);

  
protected:

  std::vector<rsDraggableNode*> nodes;
  //int draggedNodeIndex = -1;  // -1 is code for "none"
  int selectedNodeIndex = -1; // -1 is code for "none"
  float dotSize = 8;

  std::vector<rsNodeEditorObserver*> observers;

  RAPT::rsCoordinateMapper2D<double> xyMapper; // converts to/from pixel-coodinates

  // functions to convert node x,y-coordinates to text for display:
  juce::String (*xToString) (double x) = valueToString2;
  juce::String (*yToString) (double y) = valueToString2;
  // todo: make selectable from client code (setStringConversionX/Y)
  // int valuePlotMode = 0; // todo: none/compact/verbose

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeEditor)
};

//=================================================================================================

// i think, we need a mutex lock...yes, definitely

class rsNodeBasedFunctionEditor : public rsNodeEditor
{

public:

  rsNodeBasedFunctionEditor(RAPT::rsNodeBasedFunction<double>* functionMapper = nullptr, 
    CriticalSection* lockToUse = nullptr);

  ~rsNodeBasedFunctionEditor();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the rsNodeBasedFunction object to be edited. */
  virtual void setFunctionToEdit(RAPT::rsNodeBasedFunction<double>* func);

  /** Sets the mutex lock that we use to access the rsNodeBasedFunction when inserting/removing 
  and moving around nodes. */
  virtual void setMutexToUse(CriticalSection* newMutex) { lock = newMutex; }

  /** Decides whether or not the x,y coordinates should be clipped at their respective min/max
  values. */
  void setClipCoordinatesToRange(bool shouldClip) { clipRanges = shouldClip; }

  /** Updates the inherited array of draggable nodes in order to in sync with the 
  rsNodeBasedFunction object that is edited. */
  void updateDraggableNodesArray();

  /** Sets the shape type for the given node. */
  void setNodeShapeType(rsDraggableNode* node, int newType);

  /** Sets the shape parameter for the given node */
  void setNodeShapeParam(rsDraggableNode* node, double newParam);

  //-----------------------------------------------------------------------------------------------
  // \name Overrides

  virtual void paint(Graphics& g) override;
  //virtual rsDraggableNode* addNode(double pixelX, double pixelY) override;
  virtual int addNode(double pixelX, double pixelY) override;
  virtual bool removeNode(int index) override;
  virtual int moveNodeTo(int index, int pixelX, int pixelY) override;
  virtual int nodeChanged(int nodeIndex) override;
  //virtual void nodeChanged(const rsDraggableNode* node) override; // obsolete?
  virtual bool isNodeRemovable(int i) override { return (int)valueMapper->isNodeRemovable((size_t)i); }
  virtual int constrainNode(int i) override { return (int)valueMapper->constrainNode((size_t)i); }

  //-----------------------------------------------------------------------------------------------
  // \name 




protected:

  /** Structure, used to represent the parameters (x- and y-coordinates, shape-settings) of a node.
  These settings are controlled via Parameter objects. */
  struct NodeParameterSet
  { 
    NodeParameterSet(rsDraggableNode* _node)
    {
      node = _node;
      x = new Parameter("X");
      y = new Parameter("Y");
      //shapeAmount = new Parameter("ShapeValue", -1, +1, 0, Parameter::LINEAR);
      shapeAmount = new Parameter("ShapeValue", -0.99, +0.99, 0, Parameter::LINEAR);
      shapeType = new Parameter("ShapeType", 0, 1, 0, Parameter::STRING, 1);
      shapeType->addStringValue("Linear");
      shapeType->addStringValue("Exponential");
      shapeType->addStringValue("Rational");
    }
    ~NodeParameterSet()
    {
      delete x;
      delete y;
      delete shapeType;
      delete shapeAmount;
    }
    Parameter *x, *y, *shapeType, *shapeAmount; 
    rsDraggableNode* node; // pointer to the node to which the parameters belong
  };




  /** Clips the x,y coordinates (given as model coordinates) to their respective min/max values as
  set up in our xyMapper, if this clipping option is selected (via setClipCoordinatesToRange). */
  void clipIfDesired(double* x, double* y);

  /** Creates the Parameter objects associated with the given node (assumed to have been just 
  added) and adds them to our array. It returns a pointer to the new parameter set to allow 
  subclasses do additional stuff with it in their overrides. */
  //virtual void addNodeParameters(rsDraggableNode* node);
  virtual NodeParameterSet* addNodeParameters(rsDraggableNode* node); 

  /** Called from addNodeParameters to set up the callback functions to be called when the node 
  shape settings are changed. */
  void setupNodeParameterCallbacks(rsDraggableNode* node, NodeParameterSet* params);

  /** Removes the parameters for the given node. */
  void removeNodeParameters(rsDraggableNode* node);

  /** Removes the node parameters at a given index in our nodeParams array. Note that this index
  may not match the node's index - this is because when nodes are moved around, they may have to be
  re-indexe. hmmm...this is ugly ..maybe override reIndexNode in order to keep the arrays in 
  sync. */
  void removeNodeParameters(int indexInParameterArray);

  /** Adds parameters for all nodes that do not yet have parameters associated with them. Called in 
  constructor. */
  void addParametersForAllNodes();

  /** Calls removeNodeParameters for all nodes. */
  void removeParametersForAllNodes();

  /** Clears the nodes array and does all associated clean up work. */
  void clearNodes();

  /** Maps from a value from RAPT::rsFunctionNode::shapes to the corresponding index in our 
  shapeOptions array. Its a bit clunky but we need a different set of shape options here and there,
  so this mapping is needed. */
  int getShapeOptionIndex(int shapeOption);

  // data:

  std::vector<NodeParameterSet*> nodeParams;  // array of parameters associated with each node
  // maybe move the parameter handling into baseclass or into intermediate rsParameterNodeEditor 
  // class

  RAPT::rsNodeBasedFunction<double>* valueMapper = nullptr;

  std::vector<int> shapeOptions;

  bool clipRanges = false; // clip x,y to their min/max values (todo: maybe have 4 separate flags
                           // for low/high x/y clipping)
  // maybe move to baseclass, have also fixFirstNodeX, fixFirstNodeY, fixLastNodeX, fixLastNodeY
  // firstNodeRemovable, lastNodeRemovable

  CriticalSection* lock = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeBasedFunctionEditor)
};




#endif
