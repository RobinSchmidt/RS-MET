#ifndef jura_ClassConversions_h
#define jura_ClassConversions_h

/** This file contains functions to convert between corresponding datatypes in different class 
libraries. 
todo: maybe get rid of rosic::String - use std::string instead
*/

/** Converts a rosic::String into a juce::String. */
juce::String rosicToJuce(const rosic::String &stringToConvert);

/** Converts a juce::String into a rosic::String. */
rosic::String juceToRosic(const juce::String &stringToConvert);


#endif