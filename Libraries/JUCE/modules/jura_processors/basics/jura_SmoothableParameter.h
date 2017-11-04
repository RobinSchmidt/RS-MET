#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

/*

concepts/ideas:
-the smoother must call ModulatedParameter::setUnmodulatedValue
 ->SmoothableParameter should derive from ModulatableParameter?

*/


class rsSmoothingManager;
class rsSmoother;

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

  /** Function that will be called when smoothing will start on this target. The baseclass 
  implementation just sets the isSmoothing flag to true. if you need additional actions to be 
  performed, you may override it. */
  virtual void smoothingWillStart() { isSmoothing = true; }

  /** Function that will be called when smoothing has ended on this target.
  @see smoothingWillStart */
  virtual void smoothingHasEnded()  { isSmoothing = false; }

  /** Sets up the smoothing manager to be used. This functions should be called soon after 
  construction. */
  void setSmoothingManager(rsSmoothingManager* newManager) { smoothingManager = newManager; }

  /** Returns the desired smoothing time for this object (in milliseconds). */
  inline double getSmoothingTime() { return smoothingTime; }

protected:

  double smoothingTime = 100.0; // in milliseconds
  bool isSmoothing = false;
  rsSmoothingManager* smoothingManager = nullptr;

private:

  rsSmoother* smoother = nullptr;

  friend class rsSmoothingManager;
  friend class rsSmoother;

};

//=================================================================================================

/** An object that takes a SmoothingTarget and repeatedly calls setSmoothedValue on it with a value 
that is obtained from applying a lowpass filter to the target value. */

class JUCE_API rsSmoother
{

public:

  /** Constructor. */
  rsSmoother();

  /** Assigns this Smoother to a new SmoothingTarget object on which setSmoothedValue will be 
  called. */
  void setSmoothingTarget(rsSmoothingTarget* newTarget)
  {
    target = newTarget;
    target->smoother = this;
  }

  /** Sets the target value that should be reached in our SmoothingTarget. */
  void setTargetValue(double newTargetValue)
  {
    targetValue = newTargetValue;
    tolerance   = jmax(fabs(targetValue) * relativeTolerance, absoluteTolerance);
  }

  /** Sets time-constant (in milliseconds) and sampleRate for the underlying smoothin filter. */
  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    smoothingFilter.setTimeConstantAndSampleRate(0.001*timeConstant, sampleRate);
  }

  // setSmoothingOrder, setSmoothingShape, 


  /** Returns a pointer to our SmoothingTarget. */
  rsSmoothingTarget* getSmoothingTarget() { return target; }

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

  /** Iterates through our array of active smoothers and lets each of them perform its smoothing
  update operation. Should be called by outside code once per sample before the actual dsp-code for 
  that same sample is computed. If modulation is also desired, it should be called before 
  applyModulations is called on the ModualtionManager object. */
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

/** Subclass of ModulatableParameter that is also a subclass of rsSmoothingTarget to allows 
smoothing of user input to be performed. This smoothed user input is used to set up the 
unmodulatedValue in the ModulationTarget baseclass. */

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