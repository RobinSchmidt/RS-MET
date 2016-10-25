#ifndef jura_Parameter_h
#define jura_Parameter_h

class JUCE_API Parameter; // forward declaration

//=================================================================================================
// class ParameterObserver:

/**

This class is the baseclass for objects that are interested in callbacks from Parameter objects via
parameterChanged. Subclasses must implement this purely virtual callback method which will be 
called back whenever a Parameter's value was changed. Furthermore, they must implement and the 
parameterIsGoingToBeDeleted method which will be called from the destructor of Parameter objects.
If your subclass maintains pointers to the observed parameters, it should invalidate them inside 
this callback and not use them anymore thereafter.

*/

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

  /** A callback that will be called from the destructor of the Parameter class - it is intended to 
  give ParameterObserver objects an opportunity to invalidate any pointers to a particular 
  Parameter object that they may hold. */
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted) = 0;

  /** The callback method that will get called when one of our observed parameters was changed. */
  virtual void parameterChanged(Parameter* parameterThatHasChanged) = 0;

  /** The callback method that will get called when one of our observed parameters has changed its 
  range. */
  virtual void parameterRangeChanged(Parameter* parameterThatHasChangedRange) { }

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

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// the actual Parameter class:

/**

This class represents a parameter that may notify interested objects whenever its value was 
changed. To this end, it uses the observer pattern and can be used notify all objects that are 
subclasses of ParameterObserver. To put this to work, the interested object must register itself 
via registerParameterObserver. When this is done, it will be called back through its 
parameterChanged callback method whenever the value has changed. It may also de-register itself 
later via deRegisterParameterObserver. Such observer objects often like to maintain pointers to the 
observees, mainly to compare the stored pointers to the argument of the callback function in order 
to trigger particular actions depending on which particular Parameter was changed. In order to 
avoid dangling pointers inside ParameterObserver objects when Parameter obejcts get deleted, the 
Parameter class calls another callback parameterIsGoingToBeDeleted in its destructor to give the 
obervers an opportunity to invalidate their pointers.

\todo: getStorageString, setFromStorageString, getDisplayString
-allows for storing string/enum parameters with different stings than what is displayed in the 
 ComboBox:
 example:
 displayString == 24 dB/oct (if we are already in a submenu 'Lowpass')
 storageString == Lowpass24 or Lowpass>24
-include a description-string ...will be copied into the corresponding widget's description on 
 assignment

\todo
-for performance, provide two versions of setValue - the regular one and a version 
 setValueWithoutLocking - where the latter should be called only under the premise that the 
 mutex-lock is already held by the caller

\todo factor out a ParameterMapper class that does all the mapping - maybe it should be more 
generally a ValueMapper or something

\todo: maybe use a std::vector for the observers, a std::mutex for the critical section, and a
std::string for the strings, such that the class becomes independent from JUCE, maybe the parameter
class should not even deal with mutex-stuff - that should better be dealt with at the level of
the ParameterObserver


*/

class JUCE_API Parameter
{

public:

  /** An enumeration of the available scalings of the parameter. BOOLEAN will map all values >= 0.5 
  to 1.0, and all values < 0.5 to 0.0, INTEGER will linearly map to the range between minValue and 
  maxValue and round to the nearest integer, LINEAR will provide a linear mapping without rounding 
  and EXPONENTIAL will provide exponential mapping. */
  enum scalings
  {
    BOOLEAN = 0,
    STRING,
    INTEGER,
    LINEAR,
    EXPONENTIAL,
    LINEAR_BIPOLAR
    // EXPONENTIAL_WITH_OFFSET
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  Parameter(CriticalSection *criticalSectionToUse, const juce::String& newName, 
    double newMinValue = 0.0, double newMaxValue = 1.0, double newInterval = 0.0, 
    double newDefaultValue = 0.0, int newScaling = LINEAR);

  /** Destructor. */
  virtual ~Parameter();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the value of the parameter and (optionally) notifies all the observers. */
  virtual void setValue(double newValue, bool sendNotification);

  /** Resets the value of the parameter to its default value and (optionally) notifies all the 
  observers. */
  virtual void resetToDefaultValue(bool sendNotification);

  /** Sets up string-based parameters to one of the stringValues. If the passed string is found in 
  the array stringValues, it will set up the value to the index of the string. Otherwise it will 
  set it to the default value/index. */
  virtual void setStringValue(const juce::String& newString, bool sendNotification);

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
  virtual void setDefaultValues(juce::Array<double> newDefaultValues)
  {
    ScopedPointerLock spl(mutex);
    defaultValues = newDefaultValues;
  }

  /** Chooses one of the scaling methods for this parameter. @see: scalings */
  virtual void setScaling(int newScaling);

  /** Chooses one of the scaling methods for this parameter from a string. @see: scalings */
  virtual void setScalingFromString(juce::String newScalingString);

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

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the current value of the parameter. */
  virtual double getValue() const { ScopedPointerLock spl(mutex); return value; }

  /** Returns the currently chosen string-value for string based parameters. When this parameter is 
  not actually a string based parameter, it will return String::empty. */
  virtual juce::String getStringValue() const;

  /** Returns one of the option strings with the given index (or String::empty when the index is 
  out of range or this is not a string based parameter). */
  virtual juce::String getOptionStringAtIndex(int index) const;

  /** Returns the number of available options for string-based parameters. */
  virtual int getNumStringValues() const 
  { 
    ScopedPointerLock spl(mutex); 
    return stringValues.size(); 
  }

  /** Returns the lower limit of the parameter. */
  virtual double getMinValue() const { ScopedPointerLock spl(mutex); return minValue; }

  /** Returns the upper limit of the parameter. */
  virtual double getMaxValue() const { ScopedPointerLock spl(mutex); return maxValue; }

  /** Returns the interval in which the parameter is adjusted. */
  virtual double getInterval() const { ScopedPointerLock spl(mutex); return interval; }

  /** Returns the default value of the parameter. */
  virtual double getDefaultValue() const { ScopedPointerLock spl(mutex); return defaultValue; }

  /** Returns a pointer to an juce:Array of double values that are the default values for this 
  parameter. */
  virtual juce::Array<double> getDefaultValues() const 
  { 
    ScopedPointerLock spl(mutex); 
    return defaultValues; 
  }

  /** Returns the default value for string-based parameters. */
  virtual juce::String getDefaultStringValue() const;

  /** Informs, whether or not the current value of the parameter is equal to the default value. */
  virtual bool isCurrentValueDefaultValue() const 
  { 
    ScopedPointerLock spl(mutex); 
    return defaultValue == value; 
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
  AutamatableModule) - it may be desirable to turn this off, when the parameter is saved by other 
  means in order to avoid redundancies in the preset files. mmmhhh...this seems kinda bad 
  design.... */
  virtual bool shouldBeSavedAndRecalled() { ScopedPointerLock spl(mutex); return saveAndRecall; }

  /** Returns, whether or not this is a string-based parameter, that is: a parameter which can take 
  on a number of named values such as used by comboboxes. */
  virtual bool isStringParameter() const 
  { 
    ScopedPointerLock spl(mutex); 
    return stringValues.size() != 0; 
  }

  //-----------------------------------------------------------------------------------------------
  // functions for the callback-mechanisms:

  /** Adds an ParameterObserver to this parameter - the observer will get notified about parameter 
  changes. */
  virtual void registerParameterObserver(ParameterObserver *observerToAdd);

  /** Removes an ParameterObserver from this parameter - this should be called for each observer 
  before it is destroyed (like, for example, sliders on a GUI which is closed). */
  virtual void deRegisterParameterObserver(ParameterObserver *observerToRemove);


  juce_UseDebuggingNewOperator;

protected:

  /** Notifies the observers that this parameter has been changed. */
  virtual void notifyObservers();

  /** Restrics a value to the permitted parameter range as defined by minValue and maxValue. */
  virtual double restrictValueToParameterRange(double valueToRestrict);

  /** Restricts the value, defaultValue, and automation limits to the range of the parameter, makes 
  sure that max = min, makes sure that values > 0 for exponential scaling, etc.. */
  virtual void valueSanityCheck();

  juce::String name;                 // string for the parameter name
  double       value;                // actual value of the parameter
  double       minValue;             // lower limit for the parameter
  double       maxValue;             // upper limit for the parameter
  double       interval;             // interval for adjustments
  double       defaultValue;         // default value of this parameter
  int          scaling;              // index to the scaling/mapping to be used
  bool         saveAndRecall;        // flag, to switch automatic saving on/off - why?

  // array of some more default values, meant to be used for easy access via popup menu:
  juce::Array<double> defaultValues;

  // array of strings to be used for enum-based parameters (for comboboxes, etc.):
  StringArray stringValues;

  CriticalSection *mutex;


  //juce::Array<const ParameterObserver*> parameterObservers;
  juce::Array<ParameterObserver*> parameterObservers;
  //std::vector<ParameterObserver*> parameterObservers;
   // \todo switch back to juceArray - std::vector was only for debugging purposes

};

//=================================================================================================
// classes for handling the management of a variable number of parameters:

class JUCE_API ParameterSetHolder;

/**

This class is the baseclass for objects that are interested in callbacks from ParameterSetHolder 
objects that can hold a variable number of parameters. The idea is that one object holds the 
Parameters as members (possibly in arrays or other datastructures) and another object may maintain
pointers to some of these Parameter-members - for example, an equalizer may have a variable number 
of bands (each with a set of per-band parameters) and a GUI object such as a frequency-response 
editor maintains pointers to the parameters to the currently selected band. Now, whenever a 
parameter is added or removed, memory re-allocations may occur which would invalidate these 
pointers. Therefore, the object that holds the parameters should send out a callback to all 
attached observers such that they can re-assign or NULL their pointers.

*/

class JUCE_API ParameterSetObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** The callback method that will is to be called when parameters (possibly) were 
  de-/re-allocated. */
  virtual void parameterSetChanged(ParameterSetHolder* parameterSetHolderThatHasChanged) = 0;

  //===============================================================================================
  juce_UseDebuggingNewOperator;

};


class JUCE_API ParameterSetHolder
{

public:

  //-----------------------------------------------------------------------------------------------
  // de-/registering and notification:

  /** Registers a ParameterSetObserver that will be called back whenever parameters are de- or 
  re-allocated. */
  virtual void registerParameterSetObserver(ParameterSetObserver *observerToRegister);

  /** De-registers a previously registered ParameterSetObserver. */
  virtual void deRegisterParameterSetObserver(ParameterSetObserver *observerToDeRegister);

  /** Calls parameterSetChanged with the passed ParameterSetHolder as argument on all our 
  registered ParameterSetObservers. */
  virtual void sendParameterSetChangeNotification(
    ParameterSetHolder* parameterSetHolderThatHasChanged);

  juce_UseDebuggingNewOperator;

protected:

  juce::Array<ParameterSetObserver*> parameterSetObservers;

};

#endif 
