#ifndef jura_RSlider_h
#define jura_RSlider_h

class RSlider;

//=================================================================================================

class JUCE_API RSliderListener // rename to RSliderObserver, move into RSlider
{
public:
  RSliderListener() {}
  virtual ~RSliderListener() {}
  virtual void rSliderValueChanged(RSlider* rSlider) = 0;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RSliderListener)
};

//=================================================================================================

/** Baseclass for objects that paint a slider. */

class JUCE_API RSliderPainter
{
public:
  RSliderPainter() {}
  virtual ~RSliderPainter() {}

  /** You need to override this method to actually paint the slider using the passed Graphics 
  object. */
  virtual void paint(Graphics& g, RSlider* slider) = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RSliderPainter)
};

//=================================================================================================

/** This is a class for horizontal sliders....

\todo: if needed, make a subclass VerticalSlider that overrides paint and the mouse-callbacks

*/

//class JUCE_API RSlider : public RWidget, public RPopUpMenuObserver
//class JUCE_API RSlider : public rsAutomatableWidget
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
  RSlider(const juce::String& sliderName = "");
  //RSlider(const juce::String& componentName = "");

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
  // maybe rename to setStringConversion - the current name is too long and that it's a function 
  // pointer is clear from the parameter type anyway and in MeteringDisplayWithText we also use
  // the shorter version already so it would be good to make that consistent

  /** Switches the slider into vertical mode of operation - the default mode is horizontal. */
  //virtual void setVertical(bool shouldBeVertical);

  /** Sets up all the slider's settings all at once. This is recommended way of setting it up. 
  Using a sequence of single calls to setRange, setValue etc. is not recommended because within 
  the sequence of calls, invalid states may occur such as: in setRange, the old/current value may 
  be outside the new range, in setValue, the old/current range may be to small to accept the new 
  value etc. Such conditions cause our valueSanityCheck() to raise an assertion. If you can, set 
  it up all at once. It's also more efficient. */
  virtual void setup(double newMinimum, double newMaximum, double newInterval,
    double newDefaultValue, jura::Parameter::Scaling newScaling, double newValue);

  /** Sets up the range of the slider, a quantization interval for the values, a default value and
  optionally initializes the current value to the default value. */
  virtual void setRange(double newMinimum, double newMaximum, double newInterval,
    double newDefaultValue, bool initToDefault = true);

  /** Sets the scaling behavior of the slider (linear, exponential, etc.)
  @see Parameter::scalings. */
  virtual void setScaling(int newScaling);

  /** Sets the current value of the slider. 
  The optional arguments are for triggering an (asynchrnous) update of the lsider itself
  ...i think, this is to make sure, that repaint gets called on the message thread when a parameter
  is changed form the audio thread...figure out.. OK - we use repaintOnMessageThread now and may 
  now get rid of these boolean parameters - but test it before we do that clean up, so we may
  roll back, if necessary */
  virtual void setValue(double newValue, 
    const bool sendUpdateMessage = true /*, const bool sendMessageSynchronously = false*/);

  /** Overriden from RWidget - sets the current value of the slider from a string and optionally
  sends out a callback message. */
  void setStateFromString(const juce::String &stateString, bool sendChangeMessage = true) override;

  /** Sets the value of the slider expressed as proportion of the slider's range, taking into
  account the scaling behaviour. The value is thus between 0...1. */
  virtual void setNormalizedValue(double newValue, const bool sendUpdateMessage = true
    /*, const bool sendMessageSynchronously = false*/);
    // todo: rename to setNormalizedValue

  /** Sets the single default value that will be used on ctrl-click. */
  virtual void setDefaultValue(double newDefaultValue);

  /** Sets the default values that are accessible via right-click menu */
  virtual void setDefaultValues(std::vector<double> newDefaultValue);

  /** Sets the current value of the slider to the default value and optionally sends out a callback
  message.. */
  virtual void setToDefaultValue(const bool sendUpdateMessage = true
    /*, const bool sendMessageSynchronously = false*/);

  /** Overrides a RWidget::assignParameter in order to retrieve some infos from the Parameter (such
  as range, default-value, etc.) and sets up the slider accordingly. It will also set the slider's 
  name equal to the parameter's name. You can change that slider name later, if you want to display
  a name different from underlying Parameter's name. But in most cases, the slider reflecting the
  parameter's name is what is desired, so copying it over may save a lot of boilerplate code 
  typing. */
  void assignParameter(Parameter* parameterToAssign) override;

  //void parameterChanged(Parameter* p) override; 
    // overrided only for debugging purposes - comment when done

  /** The callback method that will get called when one of our observed parameters has changed its
  range. */
  void parameterRangeChanged(Parameter* parameterThatHasChangedRange) override;

  /** Copies the settings (such as the range, scaling etc. but not the assignedParameter) from
  another RSlider into this one.*/
  virtual void copySettingsFrom(const RSlider* otherSlider);

  /** Sets up a painter object that can be used for custom painting. By default, this is a nullptr
  here and we will paint directly in our paint implementation but in case it is a non-nullptr, the
  actual painting will be delegated to the painter object. */
  virtual void setPainter(RSliderPainter* painterToUse) { painter = painterToUse; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the name of the slider. */
  virtual const juce::String& getSliderName() const;

  /** Returns the current value of the slider. */
  virtual double getValue() const;

  /** Returns the value of the slider expressed as proportion of the slider's range, taking into
  account the scaling behaviour. The value is thus between 0...1. */
  virtual double getNormalizedValue() const;

  /** Returns the value of the slider expressed as proportion of the slider's range. */
  virtual double getNormalizedDefaultValue() const;

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
  juce::String getStateAsString() const override;

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

  /** Notifies our slider-listeners about value-change. */
  void notifyListeners();

  /** Overrides the method inherited from RWidget */
  void updateWidgetFromAssignedParameter(bool sendChangeMessage = false) override;


protected:

  /** Returns true when this slider has an underlying parameter that can handle only calls to 
  setNormalizedValue (not setValue) due to possibly nonmonotonic mapping. */
  //virtual bool needsToSetNormalizedParameter();

  /** Returns a value that is constrained to the range of the slider. */
  virtual double constrainValue(double value) const throw();

  /** Returns a value that is constrained to the range of the slider and quantized to the interval
  as set by setRange(). */
  virtual double constrainAndQuantizeValue(double value) const throw();

  /** Constrains the members currentValue and defaultValue via constrainedValue(). */
  virtual void valueSanityCheck();

  // our listeners:
  juce::Array<RSliderListener*> listeners;  // rename into sliderObservers

  // internal state variables:
  double currentValue, minValue, maxValue, defaultValue, interval;

  double normalizedValueOnMouseDown = 0;

  int oldDragDistance;
  double dragValue = 0;

  std::vector<double> defaultValues;
  int    scaling;
  int    layout;

  // bool wrapAround = false;
  // double bipolarCenter = RS_NAN(double);

  juce::Rectangle<int> handleRectangle;
  juce::String   sliderName;
  Component      *nameRectangle; // just a dummy in order to not receive mouse-events when the user clicks on the name-field

  RSliderPainter *painter = nullptr;

private:

  // make copy constructor and assignment operator unavailable:
  RSlider(const RSlider&);
  const RSlider& operator= (const RSlider&);

  juce_UseDebuggingNewOperator;
};

#endif
