#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

//=================================================================================================

class JUCE_API rsSmoother
{

public:

  // setSmoothingTime, setSmoothingOrder, setSmoothingShape

  void setTargetValue(double newTargetValue)
  {
    targetValue = newTargetValue;
    tolerance   = jmax(targetValue * 1.e-6, 1.e-12); // ad hoc - maybe find some better formula
  }

  inline void callCallback()
  {
    if(valueChangeCallback != nullptr)
      valueChangeCallback->call(currentValue);
  }

  bool updateValue()
  {
    currentValue = smoothingFilter.getSample(targetValue);
    if(fabs(currentValue-targetValue) < tolerance) 
    {
      currentValue = targetValue;
      callCallback();
      return true;
    }
    else
    {
      callCallback();
      return false;
    }
  }

protected:

    
  double currentValue, targetValue;
  double tolerance = 1.e-8;    // maybe make user-settable, determines when we consider the
                               // current value close enough to the target value to be considered
                               // equal, i.e. the target has been reached

  RAPT::rsSmoothingFilter<double, double> smoothingFilter;

  GenericMemberFunctionCallback1<void, double> *valueChangeCallback = nullptr;

};


//=================================================================================================

class JUCE_API rsSmoothingManager
{



  void updateSmoothedValues()
  {
    for(int i = 0; i < size(usedSmoothers); i++)
    {
      bool targetReached = usedSmoothers[i]->updateValue();
      if(targetReached)
      {
        // remove smoother from usedSmoothers and put it back to the smootherPool
      }
    }
  }

protected:

  std::vector<rsSmoother*> smootherPool;
  std::vector<rsSmoother*> usedSmoothers;

};

//=================================================================================================

class JUCE_API rsSmoothableParameter : public MetaControlledParameter
{

public:

  virtual ~rsSmoothableParameter(){};


  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // maybe we need to override setProportionalValue too? ..and may some others?

protected:


  rsSmoothingManager* smoothingManager = nullptr;

};

#endif