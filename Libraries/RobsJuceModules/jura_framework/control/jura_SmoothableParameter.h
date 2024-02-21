#ifndef jura_SmoothableParameter_h
#define jura_SmoothableParameter_h

class rsSmoothingManager;
class rsSmoother;

/** Baseclass for all smoothing targets. An rsSmoother object gets passed a pointer to an object of 
(some subclass of) this class and repeatedly calls setSmoothedValue, which your subclass must 
override and take appropriate action inside the overriden function. */

class JUCE_API rsSmoothingTarget
{

public:

  rsSmoothingTarget() = default;

  /** Destructor. It will remove the smoother for this target from the smoothinManager (if any), 
  such that it will not attempt to dereference a dangling pointer. */
  virtual ~rsSmoothingTarget();

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

  /** Sets the time-constant (in milliseconds) that will be used in the smoothing filter. 
  Passing 0 switches the smoothing off completely, passing -1 will make the target use a global
  smoothing time defined in the rsSmoothingManager. */
  void setSmoothingTime(double newTime) { smoothingTime = newTime; }

  /** Returns the desired smoothing time for this object (in milliseconds). */
  inline double getSmoothingTime() const { return smoothingTime; }

  /** Returns true when this smoothing target should use the global smoothing time defined in the
  rsSmoothingManager. */
  inline bool shouldUseGlobalSmoothingTime() const { return smoothingTime == -1.0; }


protected:

  double smoothingTime = 0.0; // in milliseconds, -1.0 is code for using global smoothing time
  bool isSmoothing = false;
  rsSmoothingManager* smoothingManager = nullptr; 

private:

  rsSmoother* smoother = nullptr; // maybe this can also serve as isSmoothing flag?

  friend class rsSmoothingManager;
  friend class rsSmoother;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsSmoothingTarget)
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

  /** Sets the currentValue variable. This must be initialized to the old value before the smoother
  is invoked. */
  void setCurrentValue(double newCurrentValue)
  {
    currentValue = newCurrentValue;
    smoothingFilter.setStates(currentValue);
  }

  /** Sets time-constant (in milliseconds) and sampleRate for the underlying smoothing filter. */
  void setTimeConstantAndSampleRate(double timeConstant, double sampleRate)
  {
    smoothingFilter.setTimeConstantAndSampleRate(0.001*timeConstant, sampleRate);
  }

  /** Sets the order of the smoothing lowpass filter. Higher orders make the transition shape more
  sigmoid. An order of 1 gives the typical RC loading curve as transition shape. */
  void setSmoothingOrder(int newOrder) { smoothingFilter.setOrder(newOrder); }

  // setSmoothingShape, 


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
  
  friend class rsSmoothingManager;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsSmoother)
};

//=================================================================================================

class JUCE_API rsSmoothingManager
{

public:

  rsSmoothingManager() = default;

  ~rsSmoothingManager();

  /** Sets the sample rate at which all the smoothing filters should operate. */
  void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; }

  /** Sets a global smoothing time in milliseconds. Whether that global time is actually used, or
  the per-target time value will be used is determined by... */
  void setSmoothingTime(double newSmoothingTime) { smoothingTime = newSmoothingTime; }
   // todo: loop over used smoothers to set the up to the new time

  /** Sets the mutex lock objevt to be used for accessing our smoother arrays in a thread-safe 
  way. */
  void setMutexLock(CriticalSection* newLock) { lock = newLock; }

  /** Turns all smoothing globally off. This may be useful when initializing parameters and/or 
  recalling a state. If any smoothers are currently running at the moment, the bypass is activated,
  it will immediately set all targets to their target values. */
  void setBypassSmoothing(bool shouldBeBypassed);

  /** Adds a smoother (which is basically a little wrapper around a lowpass filter) for the given
  smoothing target to our array of active smoothers. */
  void addSmootherFor(rsSmoothingTarget* target, double targetValue, double oldValue);

  /** Removes the smoother for the given target (if any) and puts it back into the smootherPool. */
  void removeSmootherFor(rsSmoothingTarget* target);

  /** Removes a smoother from the usedSmoothers and puts it back into the smootherPool. */ 
  void removeSmoother(int index);

  /** Returns true, if there are currently any targets being smoothed. You can check this in your
  processBlock callback - if it returns false, you may do away with the per-sample smoother 
  processing. Calling code should already hold the lock for the mutex. */
  bool needsSmoothing() { return usedSmoothers.size() > 0; }

  /** Returns whether or not smoothing is gloabally bypassed. */
  bool isSmoothingBypassed() { return smoothingBypassed; }

  /** Iterates through our array of active smoothers and lets each of them perform its smoothing
  update operation. Should be called by outside code once per sample before the actual dsp-code for 
  that same sample is computed. If modulation is also desired, it should be called before 
  applyModulations is called on the ModualtionManager object. */
  void updateSmoothedValues()
  {
    ScopedLock sl(*lock);
    updateSmoothedValuesNoLock();
  }

  /** Same as updateSmoothedValues but wihtout locking the mutex. Should be used when calling code
  has the mutex already locked. */
  void updateSmoothedValuesNoLock()
  {
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

  /** If currently any smoothers are running, their targets will immediately be set to the target 
  values by calling this function (an the smoothers are removed from the active smoother array and 
  put back into the pool). */
  void flushTargetValues();


  std::vector<rsSmoother*> usedSmoothers;  
    // temporarily made public - Elan needs to loop over them when setting a global smoothing speed
    // maybe we should have a function: useGlobalSmoothingTime(bool, double)

protected:

  std::vector<rsSmoother*> smootherPool;

  double sampleRate = 44100;
  double smoothingTime = 0.0;   // in milliseconds
  CriticalSection* lock = nullptr;
  bool smoothingBypassed = false;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsSmoothingManager)
};

//=================================================================================================

/** Subclass of ModulatableParameter that is also a subclass of rsSmoothingTarget to allow 
smoothing of user input to be performed. This smoothed user input is used to set up the 
unmodulatedValue in the ModulationTarget baseclass. */

class JUCE_API rsSmoothableParameter : public Parameter, public rsSmoothingTarget
{

public:

  /** Constructor */
  rsSmoothableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, Parameter::Scaling scaling = Parameter::Scaling::LINEAR, 
    double interval = 0.0);

  /** Destructor */
  //virtual ~rsSmoothableParameter() = default;

  // override setValue too

  /** Overrides setNormalizedValue in order to use the passed newValue as target-value for smoothing 
  instead of immediatly setting it and calling the callback. */
  virtual void setNormalizedValue(double newValue, bool sendNotification, bool callCallbacks) override;

  /** Overriden from rsSmoothingTarget. This is the per-sample callback. */
  virtual void setSmoothedValue(double newValue) override;

  /** Overrides smoothingHasEnded in order to notify those ParameterObservers that want to receive
  post smoothing parameterChanged notifications. */
  virtual void smoothingHasEnded() override;

  /** Overriden to possibly store the smoothing time, if necessarry. */
  virtual void saveToXml(XmlElement* xml) const override;

  /** Overriden to possibly recall the smoothing time, if necessarry. */
  virtual void recallFromXml(const XmlElement& xml) override;

  /** Returns true, if this parameter needs smoothing. (maybe factor out into rsSmoothingTarget) */
  inline bool needsSmoothing()
  {
    return !(smoothingTime == 0.0 || smoothingManager == nullptr
      || smoothingManager->isSmoothingBypassed());
  }


  void notifyObserversPreSmoothing();
  void notifyObserversPostSmoothing();
    // maybe get rid of these - notify observers always pre-smoothing, maybe use temporary filter 
    // objects for drawing frequency responses that use the target settings

protected:

  /** Sets the smoothing target value as normalized value in the range 0..1. */
  virtual void setNormalizedTargetValue(double newTargetValue, bool sendNotification, 
    bool callCallbacks);

  bool shouldSendNotification = true; // flag to indicate if we should send a post-smoothing notification

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsSmoothableParameter)
};

#endif