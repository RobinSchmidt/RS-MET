#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

/*

concepts/ideas:
-the smoother must call ModulatedParameter::setUnmodulatedValue
 ->SmoothableParameter should derive from ModulatableParameter?

*/


//=================================================================================================

class JUCE_API rsSmoother
{

public:

  // setSmoothingOrder, setSmoothingShape, 
  

  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    smoothingFilter.setTimeConstantAndSampleRate(timeConstant, sampleRate);
  }

  void setValueChangeCallback(GenericMemberFunctionCallback1<void, double>* newCallback)
  {
    valueChangeCallback = newCallback;
  }

  void setTargetValue(double newTargetValue)
  {
    targetValue = newTargetValue;
    tolerance   = jmax(fabs(targetValue) * relativeTolerance, absoluteTolerance);
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

  RAPT::rsSmoothingFilter<double, double> smoothingFilter;

  GenericMemberFunctionCallback1<void, double>* valueChangeCallback = nullptr;
    // maybe use std::function?


  // these values determine when we consider the current value close enough to the target value to
  // be considered equal, i.e. the target has been reached:
  static double relativeTolerance;
  static double absoluteTolerance;
  double tolerance;
  
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
        // remove smoother from usedSmoothers and put it back into the smootherPool
      }
    }
  }

protected:

  std::vector<rsSmoother*> smootherPool;
  std::vector<rsSmoother*> usedSmoothers;

  double smoothingTime = 100;  // milliseconds
  double sampleRate = 44100;

};

//=================================================================================================

// maybe we need to change the inheritance hierarchy to
// Parameter < ModulatableParameter < SmoothableParameter < MetaControlledParameter
// but where will then a PolyphonicParameter go? ...we'll see

class JUCE_API rsSmoothableParameter : public ModulatableParameter //: public MetaControlledParameter
{

public:


  /** Constructor */
  rsSmoothableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);

  virtual ~rsSmoothableParameter();


  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // maybe we need to override setProportionalValue too? ..and maybe some others?

  /** This is used as target function for the callback in rsSmoother. */
  void setSmoothedValue(double newValue);

protected:




  rsSmoothingManager* smoothingManager = nullptr;

  GenericMemberFunctionCallback1<void, double>* smootherCallbackTarget = nullptr; 
    // maybe use std::function?

};

#endif