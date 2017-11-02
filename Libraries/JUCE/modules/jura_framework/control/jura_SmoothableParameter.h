#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h



class JUCE_API SmoothingManager
{


};

//=================================================================================================

class JUCE_API SmoothableParameter : public Parameter
{

public:

  virtual ~SmoothableParameter(){};


  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // maybe we need to override setProportionalValue too? ..and may some others?

protected:


  SmoothingManager* smoothingManager = nullptr;

};

#endif