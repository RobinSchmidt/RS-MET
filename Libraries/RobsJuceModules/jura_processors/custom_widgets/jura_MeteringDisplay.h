#ifndef jura_MeteringDisplay_h
#define jura_MeteringDisplay_h


//=================================================================================================

/** A component for showing measurent values auch as the current level. */

class MeteringDisplay : public RWidget
{

public:

  enum meterStyles
  {
    levelMeterStyle = 0,    
    // rename to verticalLevelMeter, used for volume metring in dB

    triangularPointerStyle, 
    // rename to bipolar/verticalTriangularPointer - used for correlation meter

    horizontalRatio,
    // Horizontal bar with display of the ratio of current/maximum on the right like 37/128.
    // Used for displaying resource usage.

    // horizontalPercent
    // Like horizontalRatio but instead of prints something like 37/128, it prints 29%.
    // ...maybe they should both use the same drawing code and the used provides a function to 
    // convert the val/max pair into a string. That flexibility is need to allow units to be 
    // displyed like 270MB/16GB for a RAM-usage meter, for example. The call both horizontalBar.
    // The stringify function by default returns an empty string.


    numMeterStyles
  };

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  MeteringDisplay() {}

  /** Constructor. */
  //MeteringDisplay(const juce::String &description);

  /** Destructor. */
  //virtual ~MeteringDisplay();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Chooses one of the styles for displaying the value. */
  virtual void setMeterStyle(int newMeterStyle);

  /** Sets the up the lower and upper value for the meter. */
  virtual void setRange(float newMinimum, float newMaximum);

  /** Sets the up the reference value which maps to the neutral position for the meter. For a 
  volume meter, this could be the 0dB mark (up to that value, we draw a gradient from green to red,
  above it, it goes to magenta). For a correlation meter, this could be zero. */
  virtual void setReferenceValue(float newReferenceValue);

  /** Sets the up the current value for the meter. This is supposed to be called repeatedly (maybe 
  with help of a juce::Timer) to update the metering display in realtime. */
  virtual void setCurrentValue(float newValue);

  //-----------------------------------------------------------------------------------------------
  // \name Callback overrides:

  void paint(Graphics &g) override;

protected:

  float minValue = 0.f, maxValue = 1.f, currentValue = 0.f, referenceValue = 0.5f; 
  int style = levelMeterStyle;

  //juce_UseDebuggingNewOperator;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MeteringDisplay);
};

//=================================================================================================

/** Subclass of MeteringDisplay that additionally displays textual information for the name of the 
quantity being measured and a numeric value that is formatted by a user supplied function. */

class MeteringDisplayWithText : public MeteringDisplay
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime:

  MeteringDisplayWithText();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Sets the name of the measured quantity that this meter displays. This name is shown in the
  display similar to how slider names are shown. */
  virtual void setMeasurementName(const juce::String& newName);

  /** Sets the function that is used to convert from the current value to a formatted string 
  representation of the numeric value for display. The funtion should take the current value itself
  and the maximum as parameters and return the string to be displayed. */
  virtual void setStringConversion(juce::String (*newFunc) (double value, double maxValue));

  //-----------------------------------------------------------------------------------------------
  // \name Callback overrides:

  void paint(Graphics &g) override;

protected:

  /** Name of the measured quantity as shown in the display. */
  juce::String measurementName;

  /** A pointer to the function which converts a value into a juce::String. */
  juce::String (*valueToString) (double value, double maxValue) = portionToString;


};


#endif   