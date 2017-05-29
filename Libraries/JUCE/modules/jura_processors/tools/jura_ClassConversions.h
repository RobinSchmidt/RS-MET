#ifndef jura_ClassConversions_h
#define jura_ClassConversions_h

/** This file contains functions to convert between corresponding datatypes in different class 
libraries. 
todo: maybe get rid of rosic::rsString - use std::string instead
*/

/** Converts a rosic::rsString into a juce::String. */
juce::String rosicToJuce(const rosic::rsString &stringToConvert);

/** Converts a juce::String into a rosic::rsString. */
rosic::rsString juceToRosic(const juce::String &stringToConvert);


#endif