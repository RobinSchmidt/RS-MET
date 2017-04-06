#ifndef jura_RSlider_h
#define jura_RSlider_h

class RSlider;

class JUCE_API RSliderListener // rename to RSliderObserver
{
public:
  virtual void rSliderValueChanged(RSlider* rSlider) = 0;
};

/** This is a class for horizontal sliders....

\todo: if needed, make a subclass VerticalSlider that overrides paint and the mouse-callbacks

*/

//class JUCE_API RSlider : public RWidget, public RPopUpMenuObserver
//class JUCE_API RSlider : public AutomatableWidget
class JUCE_API RSlider : public RWidget
{

public:

  /** An enumeration of the different layouts - this determines where the name of the slider 
  appears. */
  enum layouts
  {
    NAME_ABOVE = 0,
    NAME_LEFT,
    NAME_INSIDE
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  RSlider(const juce::String& componentName = "");

  /** Destructor. */
  virtual ~RSlider();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the name of this slider. This name will appear on the GUI. */
  virtual void setSliderName(const juce::String &newName);

  /** Chooses a layout. @see layouts */
  virtual void setLayout(int newLayout);

  /** This function is used to pass a function-pointer to the slider. This pointer has to be the 
  address of a function which has a double-parameter and a juce::String as return-value. The 
  function will be used to convert the slider-value into a string for display. */
  virtual void setStringConversionFunction(
    juce::String (*newConversionFunction) (double valueToBeConverted));

  /** Switches the slider into vertical mode of operation - the default mode is horizontal. */
  //virtual void setVertical(bool shouldBeVertical);

  /** Sets up the range of the slider, a quantization interval for the values, a default value and 
  optionally initializes the current value to the default value. */
  virtual void setRange(double newMinimum, double newMaximum, double newInterval, 
    double newDefaultValue, bool initToDefault = true);

  /** Sets the scaling behavior of the slider (linear, exponential, etc.) 
  @see Parameter::scalings. */
  virtual void setScaling(int newScaling);

  /** Sets the current value of the slider and optionally sends out a callback message. */
  virtual void setValue(double newValue, const bool sendUpdateMessage = true, 
    const bool sendMessageSynchronously = false);

  /** Overriden from RWidget - sets the current value of the slider from a string and optionally 
  sends out a callback message. */
  virtual void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true);

  /** Sets the value of the slider expressed as proportion of the slider's range, taking into 
  account the scaling behaviour. The value is thus between 0...1. */
  virtual void setProportionalValue(double newProportionalValue, 
    const bool sendUpdateMessage = true, const bool sendMessageSynchronously = false);
    // todo: rename to setNormalizedValue

  /** Sets the single default value that will be used on ctrl-click. */
  virtual void setDefaultValue(double newDefaultValue);

  /** Sets the default values that are accessible via right-click menu */
  virtual void setDefaultValues(juce::Array<double> newDefaultValue);

  /** Sets the current value of the slider to the default value and optionally sends out a callback 
  message.. */
  virtual void setToDefaultValue(const bool sendUpdateMessage = true, 
    const bool sendMessageSynchronously = false);

  /** Overrides a RWidget::assignParameter in order to retrieve some infos from the Parameter (such 
  as range, default-value, etc.) and sets up the slider accordingly. */
  virtual void assignParameter(Parameter* parameterToAssign);

  /** The callback method that will get called when one of our observed parameters has changed its 
  range. */
  virtual void parameterRangeChanged(Parameter* parameterThatHasChangedRange);

  /** Copies the settings (such as the range, scaling etc. but not the assignedParameter) from 
  another RSlider into this one.*/
  virtual void copySettingsFrom(const RSlider* otherSlider);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the name of the slider. */
  virtual const juce::String& getSliderName() const;

  /** Returns the current value of the slider. */
  virtual double getValue() const { return currentValue; };

  /** Returns the value of the slider expressed as proportion of the slider's range, taking into 
  account the scaling behaviour. The value is thus between 0...1. */
  virtual double getProportionalValue() const { return valueToProportionOfLength(currentValue); }

  /** Returns the maximum value of the slider. */
  virtual double getMaximum()  const { return maxValue; }

  /** Returns the minimum value of the slider. */
  virtual double getMinimum()  const { return minValue; }

  /** Returns the quantization interval for the values of the slider. */
  virtual double getInterval() const { return interval; }

  /** Converts a value representing a proportion of the sliders length (assumed to be in 0...1) to 
  the corresponding value of the slider. */
  virtual double proportionOfLengthToValue(double proportion) const;

  /** Converts a value of the slider into a value that represents a proportion of the sliders 
  length (in the range 0...1). */
  virtual double valueToProportionOfLength(double value) const;

  /** Overriden from RWidget - returns the state (defined as the current value) as string. */
  virtual juce::String getStateAsString() const;

  virtual const juce::String getTextFromValue(double value) const;

  virtual double getValueFromText(const juce::String& text) const;

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;
  virtual void mouseWheelMove(const MouseEvent &event, const MouseWheelDetails &wheel) override;
  virtual void mouseDoubleClick(const MouseEvent& e) override;
  virtual void resized() override;
  virtual void paint(Graphics& g) override;

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Adds a listener to this slider whcih will be called back when the slider's value has 
  changed. */
  virtual void addListener(RSliderListener* listener) throw();

  /** Removes a listener form this slider. */
  virtual void removeListener(RSliderListener* listener) throw();

  /** Overrides the method inherited from RWidget */
  virtual void updateWidgetFromAssignedParameter(bool sendChangeMessage = false);


protected:


  /** Returns a value that is constrained to the range of the slider. */
  virtual double constrainValue(double value) const throw();

  /** Returns a value that is constrained to the range of the slider and quantized to the interval 
  as set by setRange(). */
  virtual double constrainAndQuantizeValue(double value) const throw();

  /** Constrains the members currentValue and defaultValue via constrainedValue(). */
  virtual void valueSanityCheck();

  //virtual void enablementChanged();
  virtual void handleAsyncUpdate();
  virtual void triggerChangeMessage(const bool synchronous);

  // our listeners:
  juce::Array<RSliderListener*> listeners;  // rename into sliderObservers

  // internal state variables:
  double currentValue, minValue, maxValue, defaultValue, interval;
  juce::Array<double> defaultValues;
  int    scaling;
  int    layout;

  Rectangle<int> handleRectangle;
  juce::String   sliderName;
  Component      *nameRectangle; // just a dummy in order to not receive mouse-events when the user clicks on the name-field

private:

  // make copy constructor and assignment operator unavailable:
  RSlider(const RSlider&);
  const RSlider& operator= (const RSlider&);

  juce_UseDebuggingNewOperator;
};

#endif   
