
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
  for (int i = 0; i < gridIntervalStringArray.size(); ++i)
    addStringValue(gridIntervalStringArray[i]);
}