#ifndef jura_RDraggableNumber_h
#define jura_RDraggableNumber_h

class RDraggableNumber : public RSlider
{

public:

  RDraggableNumber(const juce::String& componentName = juce::String());
  virtual ~RDraggableNumber();

  // others:
  virtual void paint(Graphics& g);
  virtual void mouseDown (const MouseEvent& e);
  //virtual void mouseUp (const MouseEvent& e);
  virtual void mouseDrag (const MouseEvent& e);
  //virtual void mouseDoubleClick (const MouseEvent& e);
  //virtual void mouseWheelMove (const MouseEvent& e, float wheelIncrementX, float wheelIncrementY);

protected:

  double valueOnMouseDown;

  RDraggableNumber (const RDraggableNumber&);
  const RDraggableNumber& operator= (const RDraggableNumber&);

  juce_UseDebuggingNewOperator;
};

#endif   
