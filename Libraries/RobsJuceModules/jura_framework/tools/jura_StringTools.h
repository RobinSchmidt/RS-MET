#ifndef jura_StringTools_h
#define jura_StringTools_h

/** Converts a juce::String into an old school zero terminated c-string - the caller must take
care to free the memory associated with the pointer via delete[], when the c-string is not needed  
anymore.  */
JUCE_API char* toZeroTerminatedString(juce::String stringToConvert);

/** Cuts a string which represents a full path to the filename only. */
JUCE_API juce::String extractFileName(juce::String fullPath);

/** Creates a string with informations about an audiofile that was passed. */
JUCE_API juce::String createAudioFileInfoString(juce::File fileToCreateStringFrom);

/** Returns the sign of the value as String (+ or - or empty) */
JUCE_API juce::String getSignAsString(double value);

/** Returns a string with prescribed minimum number of digits, prepending zeros if necesarry. */
JUCE_API juce::String intToStringWithLeadingZeros(int value, int minNumDigits);

/** Converts a MidiMessage int a String.  */
JUCE_API juce::String midiMessageToString(MidiMessage message, bool addNewLine = false);

/** Converts the string to a double value. We need this because String::getDoubleValue() doesn't 
parse "-inf" correctly (it returns +inf in this case). */
JUCE_API double toDouble(const juce::String& s);

// functions for converting numeric values representing physical units into strings - an attached
// number signifies the number of decimal digits after the point to be displayed, 'Total' 
// signifies a total number of digits to be displayed - maybe they should be moved into a class:

JUCE_API juce::String amplitudeRawAndInDecibels(double amplitudeRaw);
JUCE_API juce::String beatsToStringWithUnit4(double value);
JUCE_API juce::String centsToStringWithUnit2(double value);
JUCE_API juce::String midiNoteToString(double midiNoteNumber);
//JUCE_API juce::String midiControllerToString(int midiControllerNumber);
JUCE_API juce::String decibelsToStringWithUnit(double value);
JUCE_API juce::String decibelsToStringWithUnit1(double value);
JUCE_API juce::String decibelsToStringWithUnit2(double value);
JUCE_API juce::String decibelsPerOctaveToString(double value);
JUCE_API juce::String decibelsPerOctaveToString2(double value);
JUCE_API juce::String degreesToStringWithUnit0(double value);
JUCE_API juce::String frequencyInHzAndAsNote(double frequencyInHz);
JUCE_API juce::String frequencyToNoteString(double frequencyInHz);
JUCE_API juce::String hertzToStringWithUnit1(double value);
JUCE_API juce::String hertzToStringWithUnit2(double value);
JUCE_API juce::String hertzToStringWithUnitTotal5(double value);
JUCE_API juce::String kiloToString0(double value);
JUCE_API juce::String millisecondsToStringWithUnit0(double value);
JUCE_API juce::String millisecondsToStringWithUnit2(double value);
JUCE_API juce::String octavesToStringWithUnit2(double value);
JUCE_API juce::String percentToStringWithUnit0(double value);
JUCE_API juce::String percentToStringWithUnit1(double value);
JUCE_API juce::String percentToStringWithUnit2(double value);
JUCE_API juce::String ratioToString0(double value);
JUCE_API juce::String ratioToString1(double value);
JUCE_API juce::String ratioBothFullAtCenterToString0(double value);
JUCE_API juce::String secondsToStringWithUnit2(double value);
JUCE_API juce::String secondsToStringWithUnit3(double value);
JUCE_API juce::String secondsToStringWithUnit4(double value);
JUCE_API juce::String secondsToStringWithUnitTotal4(double value);
JUCE_API juce::String semitonesToStringWithUnit2(double value);
JUCE_API juce::String semitonesToStringWithUnit1(double value);
JUCE_API juce::String valueToStringWithTotalNumDigits(double value, int totalNumDigits = 3,
  const juce::String& suffix = juce::String());
JUCE_API juce::String valueToString(double value);
JUCE_API juce::String valueToString0(double value);
JUCE_API juce::String valueToString1(double value);
JUCE_API juce::String valueToString2(double value);
JUCE_API juce::String valueToString3(double value);
JUCE_API juce::String valueToString4(double value);
JUCE_API juce::String valueToString5(double value);
JUCE_API juce::String valueToStringTotal5(double value);
JUCE_API juce::String valueToStringWithSign0(double value);
JUCE_API juce::String valueToStringWithSign1(double value);

#endif   