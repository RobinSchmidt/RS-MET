#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

/*

concepts/ideas:
-the smoother must call ModulatedParameter::setUnmodulatedValue
 ->SmoothableParameter should derive from ModulatableParameter?

*/


class rsSmoothingManager;

/** Baseclass for all smoothing targets. An rsSmoother object gets passed a pointer to an object of 
(some subclass of) this class and repeatedly calls setSmoothedValue, which your subclass must 
override and take appropriate action inside the overriden function. */

class JUCE_API rsSmoothingTarget
{

public:


  virtual ~rsSmoothingTarget() = default;

  /** Subclasses must ovveride this to accept a new, smoothed value (supposedly coming out of some
  lowpass filter). */
  virtual void setSmoothedValue(double newValue) = 0;

  //virtual void smoothingWillStart();
  //virtual void smoothingHasEnded();

  /** Sets up the smoothing manager to be used. This functions should be called soon after 
  construction. */
  void setSmoothingManager(rsSmoothingManager* newManager)
  {
    smoothingManager = newManager;
  }

  inline double getSmoothingTime() { return smoothingTime; }

protected:

  rsSmoothingManager* smoothingManager = nullptr;

  double smoothingTime = 100.0; // in milliseconds

};

//=================================================================================================

/** */

class JUCE_API rsSmoother
{

public:

  // setSmoothingOrder, setSmoothingShape, 
  

  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    smoothingFilter.setTimeConstantAndSampleRate(0.001*timeConstant, sampleRate);
  }

  /** Assigns this Smoother to a new SmoothingTarget object on which setSmoothedValue will be 
  called. */
  void setSmoothingTarget(rsSmoothingTarget* newTarget)
  {
    target = newTarget;
  }

  void setTargetValue(double newTargetValue)
  {
    targetValue = newTargetValue;
    tolerance   = jmax(fabs(targetValue) * relativeTolerance, absoluteTolerance);
  }

  /** Updates the current value via the smoothing-filter and returns true, if the target-value has 
  been reached. If so, the SmoothingManager may remove the smoother object form the array of active 
  smoothers. */
  bool updateValue()
  {
    currentValue = smoothingFilter.getSample(targetValue);
    if(fabs(currentValue-targetValue) < tolerance) 
    {
      currentValue = targetValue;
      target->setSmoothedValue(currentValue);
      return true;
    }
    else
    {
      target->setSmoothedValue(currentValue);
      return false;
    }
  }

protected:

    
  double currentValue, targetValue;

  RAPT::rsSmoothingFilter<double, double> smoothingFilter;

  rsSmoothingTarget* target = nullptr;

  // these values determine when we consider the current value close enough to the target value to
  // be considered equal, i.e. the target has been reached:
  static double relativeTolerance;
  static double absoluteTolerance;
  double tolerance;
  
};

//=================================================================================================

class JUCE_API rsSmoothingManager
{

public:

  ~rsSmoothingManager();

  void addSmootherFor(rsSmoothingTarget* target, double targetValue);
    // maybe rename to startSmoothing

  /** Removes a smoother from the usedSmoothers and puts it back into the smootherPool. */ 
  void removeSmoother(int index);

  void updateSmoothedValues()
  {
    // i think, we should lock a mutex here

    for(int i = 0; i < size(usedSmoothers); i++)
    {
      bool targetReached = usedSmoothers[i]->updateValue();
      if(targetReached) 
      {
        removeSmoother(i);
        i--; 
      }
    }
  }

protected:

  std::vector<rsSmoother*> smootherPool;
  std::vector<rsSmoother*> usedSmoothers;

  double sampleRate = 44100;

};

//=================================================================================================

// maybe we need to change the inheritance hierarchy to
// Parameter < ModulatableParameter < SmoothableParameter < MetaControlledParameter
// but where will then a PolyphonicParameter go? ...we'll see

class JUCE_API rsSmoothableParameter : public ModulatableParameter, public rsSmoothingTarget
{

public:


  /** Constructor */
  rsSmoothableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);

  virtual ~rsSmoothableParameter();

  /** Oevvrides setValue in order to use the passed newValue as target-value for smoothing instead 
  of immediatly setting it and calling the callback. */
  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // maybe we need to override setProportionalValue too? ..and maybe some others?


  virtual void setSmoothedValue(double newValue) override;

protected:

};

#endif