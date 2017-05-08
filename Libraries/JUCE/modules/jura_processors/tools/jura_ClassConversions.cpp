
juce::String rosicToJuce(const rosic::String &stringToConvert)
{
  return juce::String(stringToConvert.getRawString());
}

rosic::String juceToRosic(const juce::String &stringToConvert)
{
  char *rawString = jura::toZeroTerminatedString(stringToConvert);
  rosic::String result(rawString);
  delete[] rawString;
  return result;
}