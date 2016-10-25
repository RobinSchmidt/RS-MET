#ifndef rojue_RTimeGridComboBox_h
#define rojue_RTimeGridComboBox_h

#include "rojue_RComboBox.h"

namespace rojue
{

  /**

  This class is a RComboBox with pre-defined entries for choosing typical sync-intervals in terms of 
  whole notes.

  */

  class RTimeGridComboBox : public RComboBox
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    RTimeGridComboBox(const juce::String& componentName = juce::String(T("RTimeGridComboBox")));

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Selects the option which corresponds to some numeric value - it is assumed that the passed 
    value is very close to one of the available options (within a margin of 1.0/64.0), otherwise it 
    will select the default option (1 beat). */
    //void setValue(double newValue, bool sendMessage = true);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Converts a string (which represents one of the options) to a value. */
    //static double getValueFromString(String stringToConvert);

    /** Converts a value to a string which represents one of the options. */
    //static String getStringFromValue(double valueToConvert);

    /** Returns the numeric value which corresponds to the currently selected option. */
    double getValue();

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  };

}

#endif  
