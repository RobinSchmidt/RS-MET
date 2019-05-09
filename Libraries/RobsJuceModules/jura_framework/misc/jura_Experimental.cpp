
void ChoiceParameter::addStringValue(const juce::String& valueToAdd, int enumValue)
{
  Parameter::addStringValue(valueToAdd);
  if( enumValue > getMaxValue() )          // if necessarry, extend the value range...
    setMaxValue(enumValue);                // ...otherwise, setValue will receive wrong values
  jassert(!RAPT::rsContains(choices,enumValue)); // the enum values must be unique
  choices.push_back(enumValue);
}

/*
template<class EnumClass>
void ChoiceParameter::addStringValue(const juce::String& valueToAdd, EnumClass enumValue)
{
  Parameter::addStringValue(valueToAdd);
  int intVal = static_cast<int>(enumValue);
  choices.push_back(intVal);
}
*/

void ChoiceParameter::setValue(double newValue, bool sendNotification, bool callCallbacks)
{
  ScopedPointerLock spl(mutex);

  int dummy = 0;
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

