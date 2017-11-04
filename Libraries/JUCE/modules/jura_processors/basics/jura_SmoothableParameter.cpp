
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

void rsSmoothingManager::addSmootherFor(rsSmoothingTarget* target, double targetValue)
{
  if(target->isSmoothing)
    target->smoother->setTargetValue(targetValue);
  else
  {
    rsSmoother* smoother;
    if(size(smootherPool) > 0)
      smoother = getAndRemoveLast(smootherPool);
    else
      smoother = new rsSmoother;
    smoother->setSmoothingTarget(target);
    smoother->setTargetValue(targetValue);
    smoother->setTimeConstantAndSampleRate(target->getSmoothingTime(), sampleRate);
    append(usedSmoothers, smoother);
    target->smoothingWillStart();
  }
}

void rsSmoothingManager::removeSmoother(int index)
{
  rsSmoother* smoother = usedSmoothers[index];
  smoother->getSmoothingTarget()->smoothingHasEnded();
  remove(usedSmoothers, index);
  append(smootherPool, smoother);
}

//=================================================================================================

rsSmoothableParameter::rsSmoothableParameter(const juce::String& name, double min, double max, 
  double defaultValue, int scaling, double interval)
  : ModulatableParameter(name, min, max, defaultValue, scaling, interval)
{

}

rsSmoothableParameter::~rsSmoothableParameter()
{

}

void rsSmoothableParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  if(smoothingManager == nullptr)
    ModulatableParameter::setValue(newValue, sendNotification, callCallbacks);
  else
    smoothingManager->addSmootherFor(this, newValue);
}

void rsSmoothableParameter::setSmoothedValue(double newValue)
{
  ModulatableParameter::setUnmodulatedValue(newValue); // is this all?
}