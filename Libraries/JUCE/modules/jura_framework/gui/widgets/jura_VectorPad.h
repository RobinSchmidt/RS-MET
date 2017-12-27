#ifndef jura_VectorPad_h
#define jura_VectorPad_h  


/** This is a class for ... */

class JUCE_API rsVectorPad : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:
  
  rsVectorPad();
  virtual ~rsVectorPad();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  virtual void assignParameterX(Parameter* newParameterX);
  
  virtual void assignParameterY(Parameter* newParameterY);

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void paint(Graphics& g) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;



  
protected:

  void setParametersFromMouseEvent(const MouseEvent& e);

  void setParametersXY(double x, double y);


  Parameter *paramX = nullptr, *paramY = nullptr;

  double xMin = -1, xMax = +1, yMin = -1, yMax = +1;

  float dotSize = 16;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVectorPad)
};


#endif
