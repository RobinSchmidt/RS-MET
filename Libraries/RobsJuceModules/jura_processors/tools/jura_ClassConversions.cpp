
juce::String rosicToJuce(const rosic::rsString &stringToConvert)
{
  return juce::String(stringToConvert.getRawString());
}

rosic::rsString juceToRosic(const juce::String &stringToConvert)
{
  char *rawString = jura::toZeroTerminatedString(stringToConvert);
  rosic::rsString result(rawString);
  delete[] rawString;
  return result;
}