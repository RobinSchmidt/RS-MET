#ifndef jura_RSyncIntervalComboBox_h
#define jura_RSyncIntervalComboBox_h

/** This class is a RComboBox with pre-defined entries for choosing typical sync-intervals in 
terms of whole notes. */

class RSyncIntervalComboBox : public RComboBox
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RSyncIntervalComboBox(
    const juce::String& componentName = juce::String("RSyncIntervalComboBox"));

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Selects the option which corresponds to some numeric value - it is assumed that the passed
  value is very close to one of the available options (within a margin of 1.0/64.0), otherwise it
  will select the default option (1 beat). */
  void setValue(double newValue, bool sendMessage = true);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Converts a string (which represents one of the options) to a value. */
  static double getValueFromString(juce::String stringToConvert);

  /** Converts a value to a string which represents one of the options. */
  static juce::String getStringFromValue(double valueToConvert);

  /** Returns the numeric value which corresponds to the currently selected option. */
  double getValue();

  juce_UseDebuggingNewOperator;

};

#endif  
