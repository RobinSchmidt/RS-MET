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

  /** Sets the position of this node (in model coordinates) calls the nodeChanged function of 
  our editor. */
  virtual void setPosition(double newX, double newY, bool callNodeChanged = true);

  /** This function is used to set up the "index" member variable which should always reflect the
  index in the "nodes" array in the editor. You need to take care of this in subclasses when 
  adding, dragging and removing nodes. */
  inline void setIndex(int newIndex) { index = newIndex; }

  inline void incrementIndex() { index++; }

  inline void decrementIndex() { index--; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  inline double getX() const { return x; }

  inline double getY() const { return y; }

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  /** Overriden in order to cause an update/repaint in our editor. */
  virtual void parameterChanged(Parameter* p) override;

protected:
   
  int index = -1;
  double x, y;                                    // (model) coordinates of the node 
  Parameter *paramX = nullptr, *paramY = nullptr; // parameter objects associated with x and y
  rsNodeEditor* nodeEditor;                       // editor which edits this node 
                                                  // todo: maybe allow more than one editor

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsDraggableNode)
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

  /** Removes the node with given index. */
  virtual void removeNode(int index);

  /** If there is a node at the given pixel position, this function will remove it (otherwise it 
  will have no effect). */
  virtual void removeNodeAt(int pixelX, int pixelY);

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

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns a pointer to a node at the given position, if any. If there's no node at the given
  position, it returns a nullptr. If there are several nodes at the given position, it will return
  the first one that is found. 
  todo: if there are several nodes, return the first one among those that are closest to the given 
  pixel position */
  rsDraggableNode* getNoteAt(int pixelX, int pixelY);

  /** Returns the index in our array of nodes a node at the given position. If there is none, it 
  will return -1. */
  int getNodeIndexAt(int pixelX, int pixelY);

  /** Returns the x-coordinate in pixels of the given node. */
  float getPixelX(const rsDraggableNode* node);

  /** Returns the y-coordinate in pixels of the given node. */
  float getPixelY(const rsDraggableNode* node);

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
  // \name Misc

  /** Darws all the nodes. You may want to override it if you need special drawing. */
  virtual void drawNodes(Graphics& g);

  
protected:

  std::vector<rsDraggableNode*> nodes;
  int draggedNodeIndex = -1;  // -1 is code for "none"
  float dotSize = 8;

  RAPT::rsCoordinateMapper2D<double> xyMapper; // converts to/from pixel-coodinates

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeEditor)
};

//=================================================================================================

// i think, we need a mutex lock...yes, definitely

class rsNodeBasedFunctionEditor : public rsNodeEditor
{

public:

  rsNodeBasedFunctionEditor(RAPT::rsNodeBasedFunction<double>* functionMapper = nullptr, 
    CriticalSection* lockToUse = nullptr);

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

  //-----------------------------------------------------------------------------------------------
  // \name Overrides

  virtual void paint(Graphics& g) override;
  //virtual rsDraggableNode* addNode(double pixelX, double pixelY) override;
  virtual int addNode(double pixelX, double pixelY) override;
  virtual void removeNode(int index) override;
  int moveNodeTo(int index, int pixelX, int pixelY) override;
  //virtual void nodeChanged(const rsDraggableNode* node) override;
  virtual int nodeChanged(int nodeIndex) override;


protected:

  /** Updates the inherited array of draggable nodes in order to in sync with the 
  rsNodeBasedFunction object that is edited. */
  void updateDraggableNodesArray();

  /** Clips the x,y coordinates (given as model coordinates) to their respective min/max values as
  set up in our xyMapper, if this clipping option is selected (via setClipCoordinatesToRange). */
  void clipIfDesired(double* x, double* y);

  RAPT::rsNodeBasedFunction<double>* valueMapper = nullptr;

  bool clipRanges = false; // clip x,y to their min/max values (todo: maybe have 4 separate flags
                           // for low/high x/y clipping)
  // maybe move to baseclass, have also fixFirstNodeX, fixFirstNodeY, fixLastNodeX, fixLastNodeY
  // firstNodeRemovable, lastNodeRemovable

  CriticalSection* lock = nullptr;



  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeBasedFunctionEditor)
};




#endif
