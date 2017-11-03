
double rsSmoother::relativeTolerance = 1.e-6;
double rsSmoother::absoluteTolerance = 1.e-12;



//=================================================================================================

rsSmoothingManager::~rsSmoothingManager()
{
  int i;
  for(i = 0; i < size(usedSmoothers); i++)
    delete usedSmoothers[i];
  for(i = 0; i < size(smootherPool); i++)
    delete smootherPool[i];
}

void rsSmoothingManager::addSmootherFor(GenericMemberFunctionCallback1<void, double>* newCallback,
  double targetValue)
{
  rsSmoother* smoother;
  if(size(smootherPool) > 0)
    smoother = getAndRemoveLast(smootherPool);
  else
    smoother = new rsSmoother;
  smoother->setValueChangeCallback(newCallback);
  smoother->setTargetValue(targetValue);
  append(usedSmoothers, smoother);
}

void rsSmoothingManager::removeSmoother(int index)
{
  rsSmoother* smoother = usedSmoothers[index];
  remove(usedSmoothers, index);
  append(smootherPool, smoother);
}

//=================================================================================================

rsSmoothableParameter::rsSmoothableParameter(const juce::String& name, double min, double max, 
  double defaultValue, int scaling, double interval)
  : ModulatableParameter(name, min, max, defaultValue, scaling, interval)
{
  // create target callback that can be passed to the smoother:
  smootherCallbackTarget = new SpecificMemberFunctionCallback1<rsSmoothableParameter, void, double>
    (this, &rsSmoothableParameter::setSmoothedValue);
}

rsSmoothableParameter::~rsSmoothableParameter()
{
  delete smootherCallbackTarget;
}

void rsSmoothableParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  if(smoothingManager == nullptr)
    ModulatableParameter::setValue(newValue, sendNotification, callCallbacks);
  else
    smoothingManager->addSmootherFor(smootherCallbackTarget, newValue);
}

void rsSmoothableParameter::setSmoothedValue(double newValue)
{
  ModulatableParameter::setUnmodulatedValue(newValue); // is this all?
}