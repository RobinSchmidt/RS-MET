#ifndef jura_MeteringDisplay_h
#define jura_MeteringDisplay_h


/** A component for showing measurent values auch as the current level.

\todo: ovrride resized instead of setBounds  */

class MeteringDisplay : public RWidget
{

public:

  enum meterStyles
  {
    levelMeterStyle = 0,
    triangularPointerStyle,

    numMeterStyles
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MeteringDisplay(const juce::String &componentName);

  /** Destructor. */
  virtual ~MeteringDisplay();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Chooses one of the styles for displaying the value. */
  virtual void setMeterStyle(int newMeterStyle);

  /** Set the up the lower and upper value for the meter. */
  virtual void setRange(double newMinimum, double newMaximum);

  /** Set the up the reference value which maps to the neutral position fo the meter. */
  virtual void setReferenceValue(double newReferenceValue);

  /** Set the up the current value for the meter. */
  virtual void setCurrentValue(double newValue);

  //-----------------------------------------------------------------------------------------------
  // overrides:

  /** Overrides paint(). */
  virtual void paint(Graphics &g);

  /** Overrides setBounds(). */
  virtual void setBounds(int x, int y, int width, int height);

protected:

  int    style;
  bool   verticalMode;
  double minValue, maxValue, currentValue, referenceValue;

  juce_UseDebuggingNewOperator;
};


#endif   