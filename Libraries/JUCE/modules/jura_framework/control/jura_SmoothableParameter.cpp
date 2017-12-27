rsSmoothingTarget::~rsSmoothingTarget()
{
  if(smoothingManager)
    smoothingManager->removeSmootherFor(this);
}

//=================================================================================================

double rsSmoother::relativeTolerance = 1.e-6;
double rsSmoother::absoluteTolerance = 1.e-12;

rsSmoother::rsSmoother()
{
  smoothingFilter.setOrder(1);
}

//=================================================================================================

rsSmoothingManager::~rsSmoothingManager()
{
  int i;
  for(i = 0; i < size(usedSmoothers); i++)
    delete usedSmoothers[i];
  for(i = 0; i < size(smootherPool); i++)
    delete smootherPool[i];
}

void rsSmoothingManager::setBypassSmoothing(bool shouldBeBypassed)
{
  ScopedLock sl(*lock);
  if(shouldBeBypassed == true)
    flushTargetValues();
  smoothingBypassed = shouldBeBypassed;
}

void rsSmoothingManager::addSmootherFor(rsSmoothingTarget* target, double targetValue, 
  double oldValue)
{
  ScopedLock sl(*lock);

  if(smoothingBypassed)
  {
    target->setSmoothedValue(targetValue);
    return;
  }

  if(target->isSmoothing)
    target->smoother->setTargetValue(targetValue);
  else
  {
    rsSmoother* smoother;
    if(size(smootherPool) > 0)
      smoother = getAndRemoveLast(smootherPool);
    else
      smoother = new rsSmoother;

    smoother->setSmoothingOrder(3);
    smoother->setTimeConstantAndSampleRate(target->getSmoothingTime(), sampleRate);
      // maybe, we should call these conditionally to avoid computations when it doesn't acually 
      // change anything

    smoother->setSmoothingTarget(target);
    smoother->setTargetValue(targetValue);
    smoother->setCurrentValue(oldValue);
    append(usedSmoothers, smoother);
    target->smoothingWillStart();
  }
}

void rsSmoothingManager::removeSmootherFor(rsSmoothingTarget* target)
{
  ScopedLock sl(*lock);
  for(size_t i = 0; i < usedSmoothers.size(); i++)
    if(usedSmoothers[i]->target == target) {
      removeSmoother(i);
      return; // there should be at most one smoother per target, so we are done
    }
}

void rsSmoothingManager::removeSmoother(int index)
{
  ScopedLock sl(*lock);
  rsSmoother* smoother = usedSmoothers[index];
  smoother->getSmoothingTarget()->smoothingHasEnded();
  smoother->getSmoothingTarget()->smoother = nullptr;
  remove(usedSmoothers, index);
  append(smootherPool, smoother);
}

void rsSmoothingManager::flushTargetValues()
{
  ScopedLock sl(*lock);
  for(int i = 0; i < size(usedSmoothers); i++)
  {
    usedSmoothers[i]->target->setSmoothedValue(usedSmoothers[i]->targetValue);
    removeSmoother(i);
    i--;
  }
}

//=================================================================================================

rsSmoothableParameter::rsSmoothableParameter(const juce::String& name, double min, double max, 
  double defaultValue, int scaling, double interval)
  : Parameter(name, min, max, defaultValue, scaling, interval)
{

}

void rsSmoothableParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  if(value == newValue)
    return;
  if(smoothingTime == 0.0 || smoothingManager == nullptr)
    Parameter::setValue(newValue, sendNotification, callCallbacks);
  else
  {
    double oldValue = getValue();

    // if the observers should be notified, we want to immediately notify MetaParameters, 
    // for example, but GUI elements such as plots should be notified only after smoothing
    // has finished (so, a frequency response plot doesnt get stuck with a graph that shows
    // the curve before smoothing has finished)
    shouldSendNotification = sendNotification;
    Parameter::setValue(newValue, false, false);
    if(sendNotification)
      notifyNonGuiObservers();
      // maybe we need an additional flag wantsNotificationAfterSmoothing ...or 
      // wantsImmediateNotification and two functions notifyImmmediately, notifyDelayed
      // or let ParameterObserver have flags preSmoothNotify, postSmoothNotify which can
      // be set individually - sliders coudl use preSmoothNotify, plotsuse postSmoothNotify
      // or both

    //Parameter::setValue(newValue, sendNotification, false); // old

    smoothingManager->addSmootherFor(this, newValue, oldValue);
  }
}

void rsSmoothableParameter::setSmoothedValue(double newValue)
{
  //modulatedValue = unmodulatedValue = value = newValue;
  value = newValue;
  callValueChangeCallbacks(); // maybe we should call a "NoLock" version of that?
}

void rsSmoothableParameter::smoothingHasEnded()
{
  rsSmoothingTarget::smoothingHasEnded();
  if(shouldSendNotification)
    notifyGuiObservers();
}
