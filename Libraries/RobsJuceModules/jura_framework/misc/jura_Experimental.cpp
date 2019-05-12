
void rsChoiceParameter::addStringValue(const juce::String& valueToAdd, int enumValue)
{
  ScopedPointerLock spl(mutex);
  stringValues.addIfNotAlreadyThere(valueToAdd);
  if( enumValue > getMaxValue() )          // if necessarry, extend the value range...
    setMaxValue(enumValue);                // ...otherwise, setValue will receive wrong values
  jassert(!RAPT::rsContains(choices, enumValue)); // the enum values must be unique
  choices.push_back(enumValue);
}

/*
template<class EnumClass>
void rsChoiceParameter::addStringValue(const juce::String& valueToAdd, EnumClass enumValue)
{
  int intValue = static_cast<int>(enumValue);
  addStringValue(valueToAdd, intValue);
}
// we get a linker error when trying to use this
*/

void rsChoiceParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);

  // the valueChangeCallbackFunction is not nullptr - why?

  //jassert(false);
  // not yet implemented - we have something to do here - make sure that the combo-box update
  // also works correctly

  // perliminary - delegate to baseclass:
  Parameter::setValue(newValue, sendNotification, callCallbacks);

  int dummy = 0;
}
// this gets called twice when changing the setting in the dropdown - figure out why and get rid
// of one of the calls. ...and then get rid of the override - it seems, we can keep the 
// baseclass implementation - it seems it works right when just delegating

juce::String rsChoiceParameter::getStringValue() const
{
  return getStringForEnumValue((int) value);
}

juce::String rsChoiceParameter::getDefaultStringValue() const
{
  return getStringForEnumValue((int) defaultValue);
}

juce::String rsChoiceParameter::getStringForEnumValue(int enumValue) const
{
  int index = (int) RAPT::rsFind(choices, enumValue);
  jassert(index < choices.size()); // enum value not found in our array - something is inconsistent
  return stringValues[index];
}

void rsChoiceParameter::setStringValue(const juce::String& newString, bool sendNotification,
  bool callCallbacks)
{
  int index = stringValues.indexOf(newString);
  setValue(choices[index], sendNotification, callCallbacks);
}

/*
-i think, we should deprecate the way, Strings are handled by implicitly mapping them to integers 
by their index in the string array (and in turn mapped to some enums in the dsp objects
-it relies on having the enum arranged in a very specific order which cannot be changed later
-maybe we should have baseclass Parameter with subclasses NumericParameter and ChoiceParameter 
where the choice is represented by a string that somehow maps to a value of some (yet unknown) 
enum-class (NumericParameter could still branch into ContinuousParameter and DiscreteParameter
...but this would probably not work well with the mod-system) but strings are not supposed to be
modulated anyway

*/

