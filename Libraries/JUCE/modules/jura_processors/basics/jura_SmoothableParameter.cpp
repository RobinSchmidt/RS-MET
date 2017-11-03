
double rsSmoother::relativeTolerance = 1.e-6;
double rsSmoother::absoluteTolerance = 1.e-12;



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
  {

  }
}

void rsSmoothableParameter::setSmoothedValue(double newValue)
{
  ModulatableParameter::setUnmodulatedValue(newValue); // is this all?
}