#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

//=================================================================================================

class JUCE_API rsSmoother
{

public:

  void update()
  {
    //currentValue = smoothingFilter.getSample(targetValue);

    // \todo:
    // -call a value-change callback, 
    // -maybe return a bool that indicates if target has been reached which can be 
    //  used by rsSmoothingManager to deactivate the smoother
  }

protected:


  double currentValue, targetValue;

  //RAPT::rsSmoothingFilter<double, double> smoothingFilter;
    // ahh..damn...we don't have access to rapt here...smoothing stuff needs to be moved to
    // jura_processors

};


//=================================================================================================

class JUCE_API rsSmoothingManager
{


protected:

  std::vector<rsSmoother*> smootherPool;
  std::vector<rsSmoother*> usedSmoothers;

};

//=================================================================================================

class JUCE_API rsSmoothableParameter : public Parameter
{

public:

  virtual ~rsSmoothableParameter(){};


  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // maybe we need to override setProportionalValue too? ..and may some others?

protected:


  rsSmoothingManager* smoothingManager = nullptr;

};

#endif