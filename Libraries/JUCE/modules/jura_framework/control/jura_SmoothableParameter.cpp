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
    if(target->shouldUseGlobalSmoothingTime())
      smoother->setTimeConstantAndSampleRate(smoothingTime, sampleRate);
    else
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
  for(int i = 0; i < size(usedSmoothers); i++)
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

void rsSmoothableParameter::setNormalizedValue(double newNormalizedValue, bool sendNotification, bool callCallbacks)
{
  if(normalizedValue == newNormalizedValue)
    return;
  if(!needsSmoothing())
    Parameter::setNormalizedValue(normalizedValue, sendNotification, callCallbacks);
  else
  {
    double oldNormalizedValue = normalizedValue;

    double tol = 1.e-7;
    //double delta = value-newValue;
    if(fabs(oldNormalizedValue-newNormalizedValue) < tol)
    {
      normalizedValue = newNormalizedValue;
      //Parameter::setNormalizedValue(normalizedValue, sendNotification, callCallbacks);
      return;
    }
    // When this parameter has an attached meta, this function gets called twice. First, when 
    // setting the parameter from a slider and a second time from the meta. In the second call,
    // getValue will return a value that is already the new value, but only up to roundoff, so
    // the if(value == newValue) check doesn't trigger. This would effectively disable smoothing,
    // so we need this additional check here.
    // -may need some more thorough checking, especially with regard to the tolerance value and if 
    //  we should also use a relative tolerance..
    // -we may also set the value to newValue and invoke a callback and notification
    // maybe we should do this in MetaControlledParameter



    // if the observers should be notified, we want to immediately notify MetaParameters, 
    // for example, but GUI elements such as plots should be notified only after smoothing
    // has finished (so, a frequency response plot doesnt get stuck with a graph that shows
    // the curve before smoothing has finished)
    shouldSendNotification = sendNotification;
    Parameter::setNormalizedValue(newNormalizedValue, false, false);
    if(sendNotification)
      notifyObserversPreSmoothing();
      // maybe we need an additional flag wantsNotificationAfterSmoothing ...or 
      // wantsImmediateNotification and two functions notifyImmmediately, notifyDelayed
      // or let ParameterObserver have flags preSmoothNotify, postSmoothNotify which can
      // be set individually - sliders could use preSmoothNotify, plotsuse postSmoothNotify
      // or both
      // ...hmm...at the moment we just notify all observers pre-and post smoothing
      // ..check, if that works well, if so, delete the notify(Non)GuiObservers functions
      // hmm - it triggers a jassert when MetaParameters are notified pre and post

    smoothingManager->addSmootherFor(this, normalizedValue, oldNormalizedValue); 
      // normalizedValue may be != newNormalizedValue now due to the fact that 
      // Parameter::setNormalizedValue may have quantized it and we should aim the smoother at the 
      // quantized value
  }
}

void rsSmoothableParameter::setSmoothedValue(double newValue)
{
  //modulatedValue = unmodulatedValue = value = newValue;
  //value = newValue;                  // old
  //value = proportionToValue(newValue); // new
  //callValueChangeCallbacks(value); // maybe we should call a "NoLock" version of that?

  callValueChangeCallbacks(proportionToValue(newValue)); // maybe we should call a "NoLock" version of that?
}

void rsSmoothableParameter::smoothingHasEnded()
{
  rsSmoothingTarget::smoothingHasEnded();
  if(shouldSendNotification)
    notifyObserversPostSmoothing();
}

void rsSmoothableParameter::saveToXml(XmlElement* xml) const
{
  Parameter::saveToXml(xml);
  // todo: if(smoothingTime != 0) ...store smoothing time
}

void rsSmoothableParameter::recallFromXml(const XmlElement& xml)
{
  Parameter::recallFromXml(xml);
  // todo: recall smoothing time, set to 0 if none is stored
}

void rsSmoothableParameter::notifyObserversPreSmoothing()
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < (int)parameterObservers.size(); i++)
  {
    if(parameterObservers[i]->wantsPreSmoothingNotification() &&
      parameterObservers[i]->wantsAutomationNotification())  // do we need this 2nd check?
      parameterObservers[i]->parameterChanged(this);
  }
}

void rsSmoothableParameter::notifyObserversPostSmoothing()
{
  ScopedPointerLock spl(mutex);
  for(int i = 0; i < (int)parameterObservers.size(); i++)
  {
    if(parameterObservers[i]->wantsPostSmoothingNotification() &&
      parameterObservers[i]->wantsAutomationNotification())  // do we need this 2nd check?
      parameterObservers[i]->parameterChanged(this);
  }
}
