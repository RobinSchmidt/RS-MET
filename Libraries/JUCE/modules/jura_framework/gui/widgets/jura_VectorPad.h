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
  
  
protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVectorPad)
};


#endif
