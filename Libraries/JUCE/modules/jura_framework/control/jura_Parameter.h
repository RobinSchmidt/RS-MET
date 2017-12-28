#ifndef jura_Parameter_h
#define jura_Parameter_h

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
  /** \name Construction/Destruction */

  /** Constructor. */
  ParameterObserver();

  /** Destructor. */
  virtual ~ParameterObserver();

  //-----------------------------------------------------------------------------------------------
  /** \name Callbacks */

  /** The callback method that will get called when one of our observed parameters was changed. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) = 0;

  /** A callback that will be called from the destructor of the Parameter class - it is intended to
  give ParameterObserver objects an opportunity to invalidate any pointers to a particular
  Parameter object that they may hold. */
  virtual void parameterWillBeDeleted(Parameter* /*parameterThatWillBeDeleted*/) {};

  /** The callback method that will get called when one of our observed parameters has changed its
  range. */
  virtual void parameterRangeChanged(Parameter* /*parameterThatHasChangedRange*/) {}

  /** Informs whether this instance wants to receive parameterChanged notifications */
  virtual bool wantsAutomationNotification();
    // rename to wantsNotification

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Turns on/off automation listening globally for all instances of ParameterObserver that have 
  the isGuiElement flag set to true. */
  static void setGuiAutomationSwitch(bool newSwitch) { guiAutomationSwitch = newSwitch; }

  /** Turns on/off automation listening globally for all instances of ParameterObserver. */
  static void setGlobalAutomationSwitch(bool newSwitch) { globalAutomationSwitch = newSwitch; }

  /** Sets a flag which indicates whether or not this ParameterObserver want to receive 
  notifications.It can make sense to switch that off temporarily when this listener also 
  manipulates the parameter. */
  void setLocalAutomationSwitch(bool newSwitch) { localAutomationSwitch = newSwitch; }

  /** Sets a flag to indicate that this ParameterObserver is a GUI element - set this to true in 
  your subclass when it is a GUI element. */
  void setIsGuiElement(bool newIsGui) { guiElement = newIsGui; }

  /** Relevant only for smoothable parameters. Sets up whether or not this observer want to receive
  notifications before smoothing starts. */
  void notifyPreSmoothing(bool shouldNotify) { preSmoothNotify = shouldNotify; }

  /** Relevant only for smoothable parameters. Sets up whether or not this observer want to receive
  notifications after smoothing has ended. */
  void notifyPostSmoothing(bool shouldNotify) { postSmoothNotify = shouldNotify; }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, if this is a GUI element (assuming that you correctly called setIsGuiElement in
  your subclass constructor). */
  bool isGuiElement() const { return guiElement; }
   // todo: maybe enforce correctness by making isGuiElement() purely virtual

  bool wantsPreSmoothingNotification()  const { return preSmoothNotify; }

  bool wantsPostSmoothingNotification() const { return postSmoothNotify; }




private:

  //-----------------------------------------------------------------------------------------------
  // data 

  static bool guiAutomationSwitch;
  static bool globalAutomationSwitch;
  bool localAutomationSwitch = true;
  bool guiElement = false;
  bool preSmoothNotify  = true;
  bool postSmoothNotify = true;

  // data is made private to enforce to set them by function calls in subclasses, so we may hook
  // into changes to them with the debugger (changing them temporarily is a source for subtle 
  // bugs)

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterObserver)
};

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

The callbacks are called before the observer notifications - this may be relevant in situations 
where a parameter change triggers a gui update that, for example, re-draws a plot. In this case, it 
may be important that a underlying dsp object (connected via callback) has already a new parameter
value when the gui repaints itself in the observer notification.

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
    IDENTITY = 0,
    BOOLEAN,
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

  /** Sometimes you want to use numeric values but represent them as strings, for example in 
  dropdown menus. in such cases, you can conveniently add a range of numeric string values via 
  this function, for example addNumericStringValues(2, 8, 2) would add "2", "4", "6", "8". */
  virtual void addNumericStringValues(int min, int max, int step = 1);

  /** Sets the mutex lock to use for parameter changes. */
  virtual void setMutexToUse(CriticalSection* cs) { mutex = cs; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the current value of the parameter. Normally, this is the value of our "value" 
  member, but subclasses may override it to wrangle the returned value before, as can be seen in 
  @see ParameterGridInterval. */
  virtual double getValue() const { ScopedPointerLock spl(mutex); return value; } // remove mutex, inline

  /** Returns the raw value of our "value" member without the possibility of any wrangling by 
  subclass overrides. */
  double getRawValue() const { ScopedPointerLock spl(mutex); return value; }

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

  /** Notifies only observers that are GUI elements. */
  //virtual void notifyGuiObservers();

  /** Notifies only observers that are not GUI elements. */
  //virtual void notifyNonGuiObservers();

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
