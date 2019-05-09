#pragma once

/* Experimental code that is under construction and not yet ready for general use */

//=================================================================================================

/** A subclass of parameter that is supposed to be used for multiple choice parameters that are 
presented as strings to the user. The Parameter baseclass handles such parameters by relying on
having the dsp objects take integer values for the choices that are defined in some enums where
the values of the enums must correspond to the index of the string in the array - like if you want
to use an option that is called "BUTTERWORTH" in some enum, the string-array in the parameter must 
have the corresponding value "Butterworth" at the index that maps to BUTTERWORTH when converted to
int. That's actually quite bad, because it implies that can't change the order the enums without 
breaking the connection (we would have to edit all places where the strings are added and bring 
them into the same order as the new enum - a maintenance nightmare)

Instead, the new way should make use of the "enum class" construct of C++11 - we don't want to rely
on implicit conversion to int anymore but instead force client code to define the connection
explicitly.....tbc... */

class JUCE_API ChoiceParameter : public Parameter
{

public:

  template<class EnumClass>
  void addStringValue(const juce::String& valueToAdd, EnumClass enumValue);
  // does it lead to code-bloat to templatize this function or will the compiler generate just
  // one function because under the hood, the enum class maps to int in every instantiation? it 
  // will have to generate another function, when enum-classes with different underlying integer
  // types are used, like enum class Color: char { red, green, blue };  ...but without the "char",
  // EnumClass will always map to "int" by default under the hood, ..hmm...

  // we also need to override setValue/getValue - it now needs to be handled in a totally
  // different way - the enumerated values don't even need to be contiguous anymore - the whole
  // range- and mapping-stuff becomes invalid ..or not?

protected:

  /** Is now protected because it should not be called from client code anymore. Instead use
  addStringValue(const juce::String& valueToAdd, EnumClass enumValue); */
  virtual void addStringValue(const juce::String& valueToAdd);


  std::vector<int> choices;

};