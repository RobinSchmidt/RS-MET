#ifndef jura_MeteringDisplay_h
#define jura_MeteringDisplay_h


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

  /** Constructor. */
  MeteringDisplay(const juce::String &componentName);

  /** Destructor. */
  virtual ~MeteringDisplay();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Chooses one of the styles for displaying the value. */
  virtual void setMeterStyle(int newMeterStyle);

  /** Sets the up the lower and upper value for the meter. */
  virtual void setRange(float newMinimum, float newMaximum);

  /** Sets the up the reference value which maps to the neutral position fo the meter. */
  virtual void setReferenceValue(float newReferenceValue);

  /** Sets the up the current value for the meter. */
  virtual void setCurrentValue(float newValue);

  //-----------------------------------------------------------------------------------------------
  // \name Overrides:

  /** Overrides paint(). */
  void paint(Graphics &g) override;


protected:


  float minValue = -48., maxValue = +6., currentValue = 0., referenceValue = 0.; 
  int style = levelMeterStyle;

  juce_UseDebuggingNewOperator;
};


#endif   