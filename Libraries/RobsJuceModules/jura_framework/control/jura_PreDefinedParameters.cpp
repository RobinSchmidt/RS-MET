
const std::vector<String> ParameterGridInterval::gridIntervalStringArray{
  "1",
  "1/2",
  "1/4",
  "1/8",
  "1/10",
  "1/12",
  "1/16",
  "1/24",
  "1/32",
  "1/64",
  "1/100",
  "1/128",
};

const std::vector<double> ParameterGridInterval::gridIntervalValueArray{
  1.0,
  1.0/2.0,
  1.0/4.0,
  1.0/8.0,
  1.0/10.0,
  1.0/12.0,
  1.0/16.0,
  1.0/24.0,
  1.0/32.0,
  1.0/64.0,
  1.0/100.0,
  1.0/128.0,
};

ParameterGridInterval::ParameterGridInterval(const juce::String& name) : Parameter(name)
{
  setScaling(STRING);
  for (size_t i = 0; i < gridIntervalStringArray.size(); ++i)
    addStringValue(gridIntervalStringArray[i]);
}

void ParameterGridInterval::setStringValue(const juce::String& newString, bool sendNotification,
  bool callCallbacks)
{
  // for debug - later, we may revert to baseclass version:
  Parameter::setStringValue(newString, sendNotification, callCallbacks); 
  String test = getStringValue();
  jassert(test == newString);
  //int dummy = 0;

  /*
  for (int i = 0; i < gridIntervalStringArray.size(); ++i)
    if(newString == gridIntervalStringArray[i])
    {
      setValue(gridIntervalValueArray[i], sendNotification, callCallbacks);
    }
  */
}

double ParameterGridInterval::getValue() const
{
  ScopedPointerLock spl(mutex); 
  int index = (int)value;
  return gridIntervalValueArray[index];
}

//-------------------------------------------------------------------------------------------------

ParameterWithKeyVelScaling::ParameterWithKeyVelScaling(const juce::String& name, double min, 
  double max, double defaultValue, Parameter::Scaling scaling, double interval)
  : Parameter(name, min, max, defaultValue, scaling, interval)
  , keyParam(name + "ByKey", -200.0, 200.0, 0.0, LINEAR, 0.01)
  , velParam(name + "ByVel", -200.0, 200.0, 0.0, LINEAR, 0.01)
{


}

// functions not yet implemented:
void ParameterWithKeyVelScaling::setKeyScaleRange(double minValue, double maxValue, 
  double defaultValue)
{

}

void ParameterWithKeyVelScaling::setVelScaleRange(double minValue, double maxValue, 
  double defaultValue)
{

}

void ParameterWithKeyVelScaling::setKeyScaleCallback()
{

}

void ParameterWithKeyVelScaling::setVelScaleCallback()
{

}