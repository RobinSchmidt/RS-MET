#ifndef jura_Parameter_h
#define jura_Parameter_h

#include <functional>

class Parameter; // forward declaration

//=================================================================================================
// class ParameterObserver:

/** This class is the baseclass for objects that are interested in callbacks from Parameter objects
via parameterChanged. Subclasses must implement this purely virtual callback method which will be
called back whenever a Parameter's value was changed. Furthermore, they must implement and the
parameterIsGoingToBeDeleted method which will be called from the destructor of Parameter objects.
If your subclass maintains pointers to the observed parameters, it should invalidate them inside
this callback and not use them anymore thereafter. */

class JUCE_API ParameterObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ParameterObserver();

  /** Destructor. */
  virtual ~ParameterObserver();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** The callback method that will get called when one of our observed parameters was changed. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) = 0;

  /** A callback that will be called from the destructor of the Parameter class - it is intended to
  give ParameterObserver objects an opportunity to invalidate any pointers to a particular
  Parameter object that they may hold. */
  virtual void parameterIsGoingToBeDeleted(Parameter* /*parameterThatWillBeDeleted*/) {};
    // rename to parameterWillBeDeleted

  /** The callback method that will get called when one of our observed parameters has changed its
  range. */
  virtual void parameterRangeChanged(Parameter* /*parameterThatHasChangedRange*/) {}

  // maybe, we should get rid of the stuff below:

  /** Informs whether this instance wants automation notifications */
  virtual bool wantsAutomationNotification();

  //-----------------------------------------------------------------------------------------------
  // public data members:

  /** A flag which can turn off automation listening globally for all instances of
  ParameterObserver that have the isGuiElement flag set to true. */
  static bool guiAutomationSwitch;

  /** A flag which can turn off automation listening globally for all instances of
  ParameterObserver. */
  static bool globalAutomationSwitch;

  /** A flag which indicates whether or not this ParameterObserver want to receive notifications.
  It can make sense to switch that off temporarily when this listener also manipulates the
  parameter. */
  bool localAutomationSwitch;

  /** A flag to indicate that this ParameterObserver is a GUI element - set this to true in your
  subclass when it is a GUI element. */
  bool isGuiElement;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterObserver)
};

//=================================================================================================

/** Objects of (subclasses of) this class are used in Parameter to do the mapping between the 
normalized range 0..1 and the actual parameter range. Subclasses can implement different mapping 
functions (shapes) to go from some minimum to some maximum parameter value. */

class JUCE_API rsParameterMapper
{

public:

  /** Constructor. */
  rsParameterMapper() {}

  /** Destructor. */
  virtual ~rsParameterMapper() {}

  /** Override this function in your subclass to map a normalized input value (in the range 0..1) 
  to the corresponding actual parameter value (in the range min..max).  */
  virtual double map(double normalizedValue) const = 0;

  /** Override this function in your subclass to map the actual parameter value (in 
  the range min..max) to the corresponding normalized value (in the range 0..1). It should be the
  inverse function of map. For example, if map is exp then unmap should be log. */
  virtual double unmap(double value) const = 0;

  /** Sets up the range for the (mapped, actual) parameter value. */
  virtual void setRange(double newMin, double newMax)
  {
    jassert(newMin <= newMax);
    min = newMin;
    max = newMax;
  }

  /** Returns the minimum value. */
  inline double getMin() const { return min; }

  /** Returns the maximum value. */
  inline double getMax() const { return max; }

protected:

  double min = 0, max = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapper)
};

/** Subclass of rsParameterMapper for linear mapping. This is appropriate for parameters that are
either intrinsically perceived on a linear scale (such as a phase between -180..+180) or a 
perceptually linearized measure of some quantity (such as decibels or semitones). */

class JUCE_API rsParameterMapperLinear : public rsParameterMapper
{
public:
  rsParameterMapperLinear() {}
  double   map(double x) const override { return min + (max-min) * x; }
  double unmap(double y) const override { return (y-min) / (max-min); }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperLinear)
};

/** Subclass of rsParameterMapper for exponential mapping. This is appropriate for parameters that 
are perceived linearly on a logarithmic scale, but nevertheless are entered as raw values by the 
user. An example would be a frequency parameter with a range of 20..20000. With exponential mapping, 
equal differences in slider value translate to equal multiplication factors for the parameter 
value. */

class JUCE_API rsParameterMapperExponential : public rsParameterMapper
{
public:
  rsParameterMapperExponential() {}
  double   map(double x) const override { return min * exp(x*(log(max/min))); }
  double unmap(double y) const override 
  {
    return jlimit(0.0, 1.0, log(y/min) / (log(max/min)));
    //return log(y/min) / (log(max/min)); 
  }

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperExponential)
};

/** A mapper based on the hyperbolic sine function, using y = a * sinh(b*x) where x is a value
between -1 and +1 (derived from the normalized 0..1 parameter p as x = 2*p-1). This function is 
suitable for parameters that should be mapped exponentially but should nevertheless by bipolar.
An example would be a frequency between -20000 and +20000 Hz. You can set up a shape parameter 
which controls the trade-off between precision around zero and high frequency precision (in the 
case of the freq-example). This shape parameter is actually the "b" in the formula. The "a" will 
then be determined by "b" and the max value (i.e. 20000). It currently supports only ranges that 
are symmetric around zero, i.e. min should be -max. */

class JUCE_API rsParameterMapperSinh : public rsParameterMapper
{
public:

  rsParameterMapperSinh(double minValue, double maxValue, double shape)
  {
    b = shape;
    setRange(minValue, maxValue); // updates the a-coeff
  }

  double map(double x) const override
  { 
    x = 2*x - 1;           // 0..1 to -1..+1
    return a * sinh(b*x);  // -max..max
    // maybe generalize to y = a * sinh(b*(x+c)) + d for unsymmetric ranges, etc. 
  } 

  double unmap(double y) const override
  { 
    y = asinh(y/a) / b;    // -1..+1
    return 0.5 * (y+1);    //  0..1
  }

  /** The range must be symmetrical around 0: newMin == -newMax. */
  void setRange(double newMin, double newMax) override
  {
    jassert(newMin == -newMax); // supports currently only 0-centered symmetric mapping
    rsParameterMapper::setRange(newMin, newMax);
    updateCoeffs();
  }

  void setShape(double newShape)
  {
    jassert(newShape > 0);
    b = newShape;
    updateCoeffs();
  }

protected:

  void updateCoeffs()
  {
    a = max / sinh(b);
  }

  double a = 1, b = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperSinh)
};

/** Parameter mapper based on the hyperbolic tangent. ...the code actually exactly parallels the 
sinh mapper - maybe we can avoid the duplication by refactoring? ...maybe it needs to use function 
pointers to sinh/asinh and tanh/atanh respectively ...and maybe that can then be generalized.
maybe factor out a class rsParameterMapperBipolar that is defined via the function:

y = a * f(b*x)

for some function f wich can be sinh, tanh, etc. and that is defined by a function pointer. 
subclasses then just set the function-pointer in the constructor (or actually 2 function pointers, 
for forward (sinh/tanh) and backward (asinh/atanh) mapping */
class JUCE_API rsParameterMapperTanh : public rsParameterMapper
{
public:

  rsParameterMapperTanh(double minValue, double maxValue, double shape)
  {
    b = shape;
    setRange(minValue, maxValue); // updates the a-coeff
  }

  double map(double x) const override
  { 
    x = 2*x - 1;           // 0..1 to -1..+1
    return a * tanh(b*x);  // -max..max
                           // maybe generalize to y = a * tanh(b*(x+c)) + d for unsymmetric ranges, etc. 
  } 

  double unmap(double y) const override
  { 
    y = atanh(y/a) / b;    // -1..+1
    return 0.5 * (y+1);    //  0..1
  }

  /** The range must be symmetrical around 0: newMin == -newMax. */
  void setRange(double newMin, double newMax) override
  {
    jassert(newMin == -newMax); // supports currently only 0-centered symmetric mapping
    rsParameterMapper::setRange(newMin, newMax);
    updateCoeffs();
  }

  void setShape(double newShape)
  {
    jassert(newShape > 0);
    b = newShape;
    updateCoeffs();
  }

protected:

  void updateCoeffs()
  {
    a = max / tanh(b);
  }

  double a = 1, b = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperTanh)
};

class JUCE_API rsParameterMapperRational : public rsParameterMapper
{
public:

  rsParameterMapperRational(double minValue, double maxValue, double shape)
  {
    tension = shape;
    rsParameterMapper::setRange(minValue, maxValue);
  }

  double curve(double t, double v) const
  {
    double tv = t*v;
    return (tv-v)/(2*tv-t-1);
  }

  /** Shape must be between -1 and +1, negative for log, positive for exp */
  void setShape(double newShape) { tension = newShape;}

  double map(double x) const override 
  {    
    return jmap(curve(tension, x), min, max);
  }
  double unmap(double y) const override 
  { 
    return curve(-tension, jmap(y, min, max, 0.0, 1.0));
  }

protected:

  double tension;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperRational)
};

class JUCE_API rsParameterMapperRationalBipolar : public rsParameterMapper
{
public:

  rsParameterMapperRationalBipolar(double minValue, double maxValue, double shape)
  {
    tension = shape;
    rsParameterMapper::setRange(minValue, maxValue);
  }

  double s_curve(double t, double v) const
  {
    double tv = t*v;

    if (v < 0.5)
      return (tv-v) / (4*tv-t-1);

    v += -0.5;
    t *= -0.5;
    tv = t*v;
    return (tv-v*0.5) / (4*tv-t-0.5) + 0.5;
  }

  /** Shape must be between -1 and +1, negative for log, positive for exp */
  void setShape(double newShape) { tension = newShape; }

  double map(double x) const override
  {
    return jmap(s_curve(tension, x), min, max);
  }
  double unmap(double y) const override
  {
    return s_curve(-tension, jmap(y, min, max, 0.0, 1.0));
  }

protected:

  double tension;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsParameterMapperRationalBipolar)
};

/*
these are elan's rational mapping functions, they could be useful for a mapper, too:

https://www.desmos.com/calculator/xaklfkriac

...also a tanh based mapper could be useful (maybe unipolar and bipolar) - a saturating curve
*/

//=================================================================================================
// the actual Parameter class:

/** This class represents a parameter that may notify interested objects whenever its value was
changed. For such notifications, two complementary callback meachnisms are provided.

The first callback mechanism uses the observer pattern and can be used notify all objects that are
subclasses of ParameterObserver. To put this to work, the interested object must register itself
via registerParameterObserver. When this is done, it will be called back through its
parameterChanged callback method whenever the value has changed. It may also de-register itself
later via deRegisterParameterObserver. Such observer objects often like to maintain pointers to the
observees, mainly to compare the stored pointers to the argument of the callback function in order
to trigger particular actions depending on which particular Parameter was changed. In order to
avoid dangling pointers inside ParameterObserver objects when Parameter obejcts get deleted, the
Parameter class calls another callback parameterIsGoingToBeDeleted in its destructor to give the
obervers an opportunity to invalidate their pointers.

The second callback mechanism relies on pointers to member-functions of the interested object. The
object that wants to be called back does not need to derive from ParameterObserver, but it must
provide a member function with one of the following signatures:
void setSomeValue(double newValue),
void setSomeValue(int newValue) or
void setSomeValue(bool newValue). Such member-functions can be
connected to a Paramter object via registerValueChangeCallback and will be called back with the
new value of the Parameter object as argument whenever the value has changed.

The second callback mechanism has the advantage that the observing object does not need to derive
from ParameterObserver and does not need to know about the existence of class Parameter at all.
Moreover, it may provide better performance than the 1st meachanism due to not having to figure
out which Parameter has been changed on callback invocation. Disadvantages are, that it dictates a
signature for the to-be-called-back member functions and potential code-bloat due to template
instantiation (it's based on template functors). It's recommended to be used mainly for performance
critical callback applications and/or when it's not acceptable to derive from ParameterObserver.

NOTE: currently there's no mechanism to de-register a member-function based callback -> the object
which member function is connected to a parameter is assumed to to be valid for the whole lifetime
of the Parameter object that calls it back (\todo provide a mechanism for de-registering member
functions)

\todo: getStorageString, setFromStorageString, getDisplayString
-allows for storing string/enum parameters with different stings than what is displayed in the
 ComboBox - example:
 displayString == 24 dB/oct (if we are already in a submenu 'Lowpass')
 storageString == Lowpass24 or Lowpass>24
-include a description-string ...will be copied into the corresponding widget's description on
 assignment

\todo
-for performance, provide two versions of setValue - the regular one and a version
 setValueWithoutLocking - where the latter should be called only under the premise that the
 mutex-lock is already held by the caller

*/

class JUCE_API Parameter
{

public:

  /** An enumeration of the available scalings of the parameter. BOOLEAN will map all
  values >= 0.5 to 1.0, and all values < 0.5 to 0.0, INTEGER will linearly map to the range between
  minValue and maxValue and round to the nearest integer, LINEAR will provide a linear mapping
  without rounding and EXPONENTIAL will provide exponential mapping. LINEAR_BIPOLAR is mapping 
  function wise the same as LINEAR, it is just distinguished to serve as a hint for sliders to 
  draw themselves differently for bipolar parameters */
  enum scalings
  {
    BOOLEAN = 0,
    STRING,         // rename to CHOICE
    INTEGER,        // maybe remove, linear should be just as good
    LINEAR,
    EXPONENTIAL,
    LINEAR_BIPOLAR,
    CUSTOM,         // requires a rsParameterMapper object to be passed

    NUM_SCALINGS
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. You need to pass a name and optionally min/max/deafult values, a type of
  scaling @see secalings and a quantization interval. */
  Parameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);

  /** DEPRECATED constructor. Use the other constructor. */
  Parameter(CriticalSection *criticalSectionToUse, const juce::String& newName,
    double newMinValue = 0.0, double newMaxValue = 1.0, double newInterval = 0.0,
    double newDefaultValue = 0.0, int newScaling = LINEAR);

  /** Destructor. */
  virtual ~Parameter();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the value of the parameter and (optionally) notifies all the listeners and calls the
  assigned callback function. */
  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks);

  /** Sets a new range and value as single operation to avoid inconsistencies that may occur when
  setting these things one after another. */
  virtual void setRangeAndValue(double newMin, double newMax, double newValue, 
    bool sendNotification, bool callCallbacks);

  /** Sets the value of the parameter where the input argument is assumed to be normalized to the
  range 0...1  .... */
  virtual void setProportionalValue(double newProportionalValue, bool sendNotification,
    bool callCallbacks);

  /** Resets the value of the parameter to its default value and (optionally) notifies all the
  listeners and calls the assigned callback function. */
  virtual void resetToDefaultValue(bool sendNotification, bool callCallbacks);

  /** Sets up string-based parameters to one of the stringValues. If the passed string is found in
  the array stringValues, it will set up the value to the index of the string. Otherwise it will
  set it to the default value/index. */
  virtual void setStringValue(const juce::String& newString, bool sendNotification,
    bool callCallbacks);

  /** Sets the lower and upper limit for the parameter's value. */
  virtual void setRange(double newMinValue, double newMaxValue);

  /** Sets a lower limit for the value of the parameter. */
  virtual void setMinValue(double newLimit);

  /** Sets an upper limit for the value of the parameter. */
  virtual void setMaxValue(double newLimit);

  /** Sets the main default value of the parameter and optionally sets the current value to this
  new default value. */
  virtual void setDefaultValue(double newDefaultValue, bool setToDefault = false);

  /** Sets multiple default values that can be easily accessed via a popup menu. */
  virtual void setDefaultValues(std::vector<double> newDefaultValues)
  {
    ScopedPointerLock spl(mutex);
    defaultValues = newDefaultValues;
  }

  /** Chooses one of the scaling methods for this parameter. @see: scalings */
  virtual void setScaling(int newScaling);

  /** Chooses one of the scaling methods for this parameter from a string. @see: scalings */
  virtual void setScalingFromString(juce::String newScalingString);

  /** Sets up a custom parameter mapper object to be used for mapping back and forth between 
  normalized (0..1) and actual (min..max) values. You can use this, whenever you need a mapping
  that is not listed in the scalings enum. This Parameter will take over ownership of the passed
  object, i.e. delete it on destruction. */
  virtual void setMapper(rsParameterMapper* newMapper);

  /** Assigns a name to this parameter. */
  virtual void setName(const juce::String& newName)
  {
    ScopedPointerLock spl(mutex);
    name = newName;
  }

  /** Sets a flag to indicate whether this parameter should be automatically saved into a preset
  (as a parameter of AutamatableModule) - it may be desirable to turn this off, when the parameter
  is saved by other means in order to avoid redundancies in the preset files. */
  virtual void setSaveAndRecall(bool newSaveAndRecall)
  {
    ScopedPointerLock spl(mutex);
    saveAndRecall = newSaveAndRecall;
  }

  /** Adds a string value to the array - relevant only for string-based parameters. */
  virtual void addStringValue(const juce::String& valueToAdd);

  /** Sets the mutex lock to use for parameter changes. */
  virtual void setMutexToUse(CriticalSection* cs) { mutex = cs; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the current value of the parameter. */
  virtual double getValue() const { ScopedPointerLock spl(mutex); return value; } // remove mutex, inline

  /** Returns the normalized value in the range 0..1. */
  inline double getProportionalValue() { return valueToProportion(value); }

  /** Converts the clear text value to a proportional value in the range 0..1 according to our
  scaling/mapping function. */
  virtual double valueToProportion(double value);

  /** Converts a proportional value in the range 0..1 to a clear text value according to our
  scaling/mapping function. */
  virtual double proportionToValue(double proportion);

  /** Returns the currently chosen string-value for string based parameters. When this parameter is
  not actually a string based parameter, it will return String::empty. */
  virtual juce::String getStringValue() const;

  /** Returns one of the option strings with the given index (or String::empty when the index is
  out of range or this is not a string based parameter). */
  virtual juce::String getOptionStringAtIndex(int index) const;

  /** Returns the number of available options for string-based parameters. */
  virtual int getNumStringValues() const
  {
    ScopedPointerLock spl(mutex); return stringValues.size();
  }

  /** Returns the lower limit of the parameter. */
  virtual double getMinValue() const
  {
    ScopedPointerLock spl(mutex); return mapper->getMin();
  }

  /** Returns the upper limit of the parameter. */
  virtual double getMaxValue() const 
  { 
    ScopedPointerLock spl(mutex); return mapper->getMax(); 
  }

  /** Returns the interval in which the parameter is adjusted. */
  virtual double getInterval() const { ScopedPointerLock spl(mutex); return interval; }

  /** Returns the default value of the parameter. */
  virtual double getDefaultValue() const { ScopedPointerLock spl(mutex); return defaultValue; }

  /** Returns a pointer to an juce:Array of double values that are the default values for this
  parameter. */
  virtual std::vector<double> getDefaultValues() const
  {
    ScopedPointerLock spl(mutex);
    return defaultValues;
  }

  /** Returns the default value for string-based parameters. */
  virtual juce::String getDefaultStringValue() const;

  /** Informs, whether or not the current value of the parameter is equal to the default value. */
  virtual bool isCurrentValueDefaultValue() const
  {
    ScopedPointerLock spl(mutex); return defaultValue == value;
  }

  /** Returns the scaling methods for this parameter. @see: scalings */
  virtual int getScaling() const { ScopedPointerLock spl(mutex); return scaling; }

  /** Returns the scaling methods for this parameter as a String. @see: scalings */
  virtual juce::String getScalingString() const;

  /** Returns the name of the parameter. */
  virtual const juce::String& getName() const { ScopedPointerLock spl(mutex); return name; }

  /** Retruns true if 'nameToCompareTo' equals the name of this parameter. */
  virtual bool hasName(const juce::String& nameToCompareTo) const
  {
    ScopedPointerLock spl(mutex); return (name == nameToCompareTo);
  }

  /** Informs whether or not this parameter is in default state, that is: are the assigned lower
  and upper limits and the assigned controller still in their default settings. */
  //bool isInDefaultState();

  /** Returns true when this parameter should be saved into preset files (as a parameter of
  AutamatableModule) - it may be desirable to turn this off, when the parameter is save by other
  means in order to avoid redundancies in the preset files. mmmhhh...this seems kinda bad design */
  virtual bool shouldBeSavedAndRecalled() { ScopedPointerLock spl(mutex); return saveAndRecall; }

  /** Returns, whether or not this is a string-based parameter, that is: a parameter which can take
  on a number of named values such as used by comboboxes. */
  virtual bool isStringParameter() const
  {
    ScopedPointerLock spl(mutex);
    return scaling == STRING && stringValues.size() != 0;
    //return stringValues.size() != 0;
  }

  //-----------------------------------------------------------------------------------------------
  // functions for the callback-mechanisms:

  /** Adds an ParameterObserver to this parameter - the observer will get notified about parameter
  changes. */
  virtual void registerParameterObserver(ParameterObserver *observerToAdd);

  /** Removes an ParameterObserver from this parameter - this should be called for each observer
  before it is destroid (like, for example, sliders on a GUI which is closed). */
  virtual void deRegisterParameterObserver(ParameterObserver *observerToRemove);

  /** Sets up a member function of some class as callback that will be called whenever the value of
  the Parameter changes. The member function must have the following signature:
  'void myMemberFunction(double newValue)' and will be called with the new value of the Parameter
  as its argument for 'newValue'. To register the member-function for a particular instance, you
  should call it like: setValueChangeCallback(&myInstance, &MyClass::myMemberFunction); where
  myInstance is the instance (and &myInstance therefore is a pointer to it). Calling it with NULL
  as first argument will reset the callback to call nothing (2nd argument does not matter in this
  case, but you still must pass something suitable). */
  template<class CalleeClass>
  void setValueChangeCallback(CalleeClass *calleeObject,
    void (CalleeClass::*memberToCall) (double))
  {
    ScopedPointerLock spl(mutex);
    if(valueChangeCallbackDouble != nullptr)
    {
      delete valueChangeCallbackDouble;
      valueChangeCallbackDouble = nullptr;
    }
    if(calleeObject != nullptr)
      valueChangeCallbackDouble =
      new SpecificMemberFunctionCallback1<CalleeClass, void, double>(calleeObject, memberToCall);
    callValueChangeCallbacks(); // enforce consistency of parameter with core object after wiring
                                // up the callback
  }

  /** @see registerValueChangeCallback(CalleeClass *calleeObject, void (CalleeClass::*memberToCall) (double)) */
  template<class CalleeClass>
  void setValueChangeCallback(CalleeClass *calleeObject, void (CalleeClass::*memberToCall) (int))
  {
    ScopedPointerLock spl(mutex);
    if(valueChangeCallbackInt != nullptr)
    {
      delete valueChangeCallbackInt;
      valueChangeCallbackInt = nullptr;
    }
    if(calleeObject != nullptr)
      valueChangeCallbackInt
      = new SpecificMemberFunctionCallback1<CalleeClass, void, int>(calleeObject, memberToCall);
    callValueChangeCallbacks();
  }

  /** @see registerValueChangeCallback(CalleeClass *calleeObject, void (CalleeClass::*memberToCall) (double)) */
  template<class CalleeClass>
  void setValueChangeCallback(CalleeClass *calleeObject, void (CalleeClass::*memberToCall) (bool))
  {
    ScopedPointerLock spl(mutex);
    if(valueChangeCallbackBool != nullptr)
    {
      delete valueChangeCallbackBool;
      valueChangeCallbackBool = nullptr;
    }
    if(calleeObject != nullptr)
      valueChangeCallbackBool =
      new SpecificMemberFunctionCallback1<CalleeClass, void, bool>(calleeObject, memberToCall);
    callValueChangeCallbacks();
  }

  void setValueChangeCallback(std::function<void(double)> cb)
  {
	  ScopedPointerLock spl(mutex);
	  valueChangeCallbackFunction = cb;
    callValueChangeCallbacks();
  }

  /** Clears our valueChangeCallbacks, so they will call back nothing. */
  void clearValueChangeCallbacks()
  {
    ScopedPointerLock spl(mutex);

	  valueChangeCallbackFunction = [](double) {};

    delete valueChangeCallbackDouble;
    delete valueChangeCallbackInt;
    delete valueChangeCallbackBool;

    valueChangeCallbackDouble = nullptr;
    valueChangeCallbackInt    = nullptr;
    valueChangeCallbackBool   = nullptr;
  }

  /** Notifies the observers that this parameter has been changed. */
  virtual void notifyObservers();

  /** Calls all currently registered callback functions for a value change. */
  virtual void callValueChangeCallbacks();

protected:


  /** Restrics a value to the permitted parameter range as defined by minValue and maxValue. */
  virtual double restrictValueToParameterRange(double valueToRestrict);

  /** Restricts the value, defaultValue, and automation limits to the range of the parameter, makes sure that max = min, makes sure that
  values > 0 for exponential scaling, etc.. */
  virtual void valueSanityCheck();

  juce::String name;                 // string for the parameter name
  double       value;                // actual value of the parameter
  double       interval;             // interval for adjustments ...rename to stepSize
  double       defaultValue;         // default value of this parameter ...maybe rename to resetValue
  int          scaling;              // index to the scaling/mapping to be used
  bool         saveAndRecall = true; // flag, to switch automatic saving on/off - why?

  // array of some more default values, meant to be used for easy access via popup menu:
  std::vector<double> defaultValues;

  // array of strings to be used for enum-based parameters (for comboboxes, etc.):
  StringArray stringValues;

  // mutex lock to acquire when we notify our observers and call our callbacks:
  CriticalSection* mutex = nullptr;

  // mapper object for mapping back and forth between normalized and actual value:
  rsParameterMapper* mapper = nullptr;

  // the callback objects
  GenericMemberFunctionCallback1<void, double> *valueChangeCallbackDouble = nullptr;
  GenericMemberFunctionCallback1<void, int>    *valueChangeCallbackInt    = nullptr;
  GenericMemberFunctionCallback1<void, bool>   *valueChangeCallbackBool   = nullptr;
  std::function<void(double)> valueChangeCallbackFunction;
  std::vector<ParameterObserver*> parameterObservers;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Parameter)
};

//=================================================================================================
// classes for handling the management of a variable number of parameters:

class ParameterSetHolder;

/** This class is the baseclass for objects that are interested in callbacks from
ParameterSetHolder objects that can hold a variable number of parameters. The idea is that one
object holds the Parameters as members (possibly in arrays or other datastructures) and another
object may maintain pointers to some of these Parameter-members - for example, an equalizer may
have a variable number of bands (each with a set of per-band parameters) and a GUI object such as
a frequency-response editor maintains pointers to the parameters to the currently selected band.
Now, whenever a parameter is added or removed, memory re-allocations may occur which would
invalidate these pointers. Therefore, the object that holds the parameters should send out a
callback to all attached observers such that they can re-assign or NULL their pointers. */

class ParameterSetObserver
{

public:


  ParameterSetObserver() {}
  virtual ~ParameterSetObserver() {}

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** The callback method that will is to be called when parameters (possibly) were
  de-/re-allocated. */
  virtual void parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged) = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterSetObserver)
};


class ParameterSetHolder
{

public:

  ParameterSetHolder() {}
  virtual ~ParameterSetHolder() {}

  //-----------------------------------------------------------------------------------------------
  // de-/registering and notification:

  /** Registers a ParameterSetObserver that will be called back whenever parameters are de- or
  re-allocated. */
  virtual void registerParameterSetObserver(ParameterSetObserver *observerToRegister);

  /** De-registers a previously registered ParameterSetObserver. */
  virtual void deRegisterParameterSetObserver(ParameterSetObserver *observerToDeRegister);

  /** Calls parameterSetChanged with the passed ParameterSetHolder as argument on all our
  registered ParameterSetObservers. */
  virtual void sendParameterSetChangeNotification(ParameterSetHolder*
    parameterSetHolderThatHasChanged);

protected:

  std::vector<ParameterSetObserver*> parameterSetObservers; // use std::vector

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterSetHolder)
};


#endif
