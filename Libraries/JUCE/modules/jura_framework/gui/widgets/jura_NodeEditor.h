#ifndef jura_NodeEditor_h
#define jura_NodeEditor_h  


class rsNodeEditor;


class rsDraggableNode : public ParameterObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction

  /** Constructor */
  rsDraggableNode(rsNodeEditor* editorContainingThisNode);

  /** Destructor */
  ~rsDraggableNode();



  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Assigns the parameter object associated with the x-coordinate of this node. */
  virtual void assignParameterX(Parameter* newParameterX);

  /** Assigns the parameter object associated with the y-coordinate of this node. */
  virtual void assignParameterY(Parameter* newParameterY);

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


/** This is a baseclass for editors that need to create/remove/drag-around a number of nodes. */

class JUCE_API rsNodeEditor : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction
  
  rsNodeEditor();
  virtual ~rsNodeEditor();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:



  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void paint(Graphics& g) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;
  virtual void nodeChanged(const rsDraggableNode* node);



  
protected:




  // these variables also appear in rsVectorPad - maybe we can factor out a baseclass:
  double xMin = -1, xMax = +1, yMin = -1, yMax = +1;
  float dotSize = 16;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeEditor)
};


#endif
