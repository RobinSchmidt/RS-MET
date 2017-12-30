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
  rsDraggableNode(rsNodeEditor* editorContainingThisNode, double pixelX, double pixelY);

  /** Destructor */
  ~rsDraggableNode();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Assigns the parameter object associated with the x-coordinate of this node. */
  virtual void assignParameterX(Parameter* newParameterX);

  /** Assigns the parameter object associated with the y-coordinate of this node. */
  virtual void assignParameterY(Parameter* newParameterY);

  /** Sets up a new pixel position for this node. */
  virtual void setPixelPosition(double newX, double newY);

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  /** Overriden in order to cause an update/repaint in our editor. */
  virtual void parameterChanged(Parameter* p) override;

protected:

  double pixelX = 0, pixelY = 0;                  // pixel coordinates of the node
  Parameter *paramX = nullptr, *paramY = nullptr; // parameter objects associated with x and y
  rsNodeEditor* nodeEditor;                       // editor which edits this node 
                                                  // todo: maybe allow more than one editor

  friend class rsNodeEditor;
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

  /** Adds a new node at the given pixel position a returns a pointer to it. The new node will use 
  nullptrs for the Parameter objects associated with the x- and y-coordinate. Subclasses may want 
  to override this in order to create and assign actual Parameter objects with associated 
  callbacks. */
  virtual rsDraggableNode* addNode(double pixelX, double pixelY);
   // maybe use int for x,y

  /** If there is a node at the given pixel position, this function will remove it (otherwise it 
  will have no effect). */
  virtual void removeNodeAt(int pixelX, int pixelY);

  /** Sets the size of the dots in pixels. */
  void setDotSize(float newDotSize);

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

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void paint(Graphics& g) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;
  virtual void mouseUp(const MouseEvent &e) override;
  virtual void mouseMove(const MouseEvent &e) override;
  virtual void nodeChanged(const rsDraggableNode* node);

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Darws all the nodes. You may want to override it if you need special drawing. */
  virtual void drawNodes(Graphics& g);

  
protected:


  std::vector<rsDraggableNode*> nodes;
  rsDraggableNode* draggedNode = nullptr; // node that is currently dragged

  // these variables also appear in rsVectorPad - maybe we can factor out a baseclass:
  //double xMin = -1, xMax = +1, yMin = -1, yMax = +1; // not yet used - maybe not needed
  float dotSize = 8;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeEditor)
};

//=================================================================================================

class rsNodeBasedFunctionEditor : public rsNodeEditor
{

public:

  //rsNodeBasedFunctionEditor() = default;

  rsNodeBasedFunctionEditor(RAPT::rsInterpolatingFunction<double>* functionMapper);


  //-----------------------------------------------------------------------------------------------
  // \name Overrides:

  virtual rsDraggableNode* addNode(double pixelX, double pixelY) override;
  virtual void removeNodeAt(int pixelX, int pixelY) override;
  virtual void nodeChanged(const rsDraggableNode* node);


protected:

  RAPT::rsInterpolatingFunction<double>* mapper = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeBasedFunctionEditor)
};

#endif
