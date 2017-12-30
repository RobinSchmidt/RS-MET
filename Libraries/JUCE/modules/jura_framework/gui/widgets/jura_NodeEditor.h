#ifndef jura_NodeEditor_h
#define jura_NodeEditor_h  


/** This is a baseclass for editors that need to create/remove/drag-around a number of nodes. */

class JUCE_API rsNodeEditor : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:
  
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



  
protected:




  // these variables also appear in rsVectorPad - maybe we can factor out a baseclass:
  double xMin = -1, xMax = +1, yMin = -1, yMax = +1;
  float dotSize = 16;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNodeEditor)
};


#endif
