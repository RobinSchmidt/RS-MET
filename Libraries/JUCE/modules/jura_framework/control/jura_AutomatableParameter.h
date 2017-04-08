#ifndef jura_AutomatableParameter_h
#define jura_AutomatableParameter_h


/** This class represents a parameter which can be automated via a MIDI controller.


\todo:
-rename to MidiControlledParameter
-rename get/setAutomationValue in get/setNormalizedValue
-rename get/setMidi... in get/setMeta  */

class JUCE_API AutomatableParameter : public Parameter
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AutomatableParameter(CriticalSection *criticalSectionToUse, const juce::String& newName, 
    double newMinValue = 0.0, double newMaxValue = 1.0, double newInterval = 0.0,
    double newDefaultValue = 0.0, int newScaling = LINEAR, int newDefaultMidiController = -1);

  /** Destructor. */
  virtual ~AutomatableParameter();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the value of the parameter and (optionally) notifies all the listeners. */
  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks)
  {
    ScopedPointerLock spl(mutex);
    Parameter::setValue(newValue, sendNotification, callCallbacks);
    updateAutomationFromParameter();
  }

  /** Sets the lower and upper limit for the parameter's value and also updates the 
  automation-limits according to the new range. */
  virtual void setRange(double newMinValue, double newMaxValue)
  {
    ScopedPointerLock spl(mutex);
    Parameter::setRange(newMinValue, newMaxValue);
    lowerAutomationLimit = newMinValue;
    upperAutomationLimit = newMaxValue;
    updateAutomationFromParameter();
  }

  /** Sets a lower limit for the value of the parameter and also sets the lower limit for the 
  automation to this value. */
  virtual void setMinValue(double newMinValue)
  {
    ScopedPointerLock spl(mutex);
    Parameter::setMinValue(newMinValue);
    lowerAutomationLimit = newMinValue;
    updateAutomationFromParameter();
  }

  /** Sets a lower limit for the value of the parameter and also sets the upper limit for the 
  automation to this value. */
  virtual void setMaxValue(double newMaxValue)
  {
    ScopedPointerLock spl(mutex);
    Parameter::setMaxValue(newMaxValue);
    upperAutomationLimit = newMaxValue;
    updateAutomationFromParameter();
  }

  /** Chooses one of the scaling methods for this parameter. @see: scalings */
  virtual void setScaling(int newScaling)
  {
    ScopedPointerLock spl(mutex);
    Parameter::setScaling(newScaling);
    updateAutomationFromParameter();
  }

  /** Sets the value of the parameter where the input argument is assumed to be normalized to the 
  range 0...1 - arguments in this range will be mapped to the range between lowerAutomationLimit 
  and upperAutomationLimit and the actual parameter will be adjusted accordingly. */
  virtual void setAutomationValue(double newAutomationValue, bool sendNotification, 
    bool callCallbacks);

  /** Sets a lower limit for the automatable range of the parameter and informs if this was 
  successful - if not, it means that the desired (new) lower limit was higher than the 
  current (old) upper limit. */
  virtual void setLowerAutomationLimit(double newLimit);

  /** Sets an upper limit for the automatable range of the parameter and informs if this was 
  successful - if not, it means that the desired (new) upper limit was lower than the 
  current (old) lower limit. */
  virtual void setUpperAutomationLimit(double newLimit);

  /** Sets up the default MIDI controller for this parameter. */
  virtual void setDefaultMidiController(int newDefaultControllerNumber);

  /** Assigns a MIDI controller to this parameter. */
  virtual void assignMidiController(int controllerNumberToAssign);

  /** Sets the assigned MIDI controller back to the default controller number. */
  virtual void revertToDefaultMidiController();

  /** Checks whether the controllerNumber matches the controller which is assigned to this 
  parameter, and if so, it sets up the automation value to controllerValue / 127.0. */
  virtual void setMidiController(int controllerNumber, float controllerValue);
  //virtual void setMidiController(int controllerNumber, int controllerValue);

  /** Switches this parameter into midi-learn mode - the very next received MIDI controller 
  will be assigned to this parameter then. */
  virtual void switchIntoMidiLearnMode(bool shouldLearn = true);

  /** Reverts all automation settings to their defaults and optionally also resets the current 
  value - if this second is done (first parameter true), there is a further option to do this 
  with or without sending a notification to the listeners. */
  //virtual void revertToDefaults(bool resetValueAlso = false, bool sendNotification = false);
  virtual void revertToDefaults(bool resetValueAlso, bool sendNotification, bool callCallbacks);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the automation value of this parameter (in the range 0...1) - when the actual 
  parameter is below the lower automation limit, it will return 0.0, when it is above the upper 
  automation limit, it will return 1.0. */
  virtual double getAutomationValue() const 
  { ScopedPointerLock spl(mutex); return automationValue; }

  /** Returns the lower limit of the automatable range of the parameter. */
  virtual double getLowerAutomationLimit() const 
  { ScopedPointerLock spl(mutex); return lowerAutomationLimit; }

  /** Returns the upper limit of the automatable range of the parameter. */
  virtual double getUpperAutomationLimit() const 
  { ScopedPointerLock spl(mutex); return upperAutomationLimit; }

  /** Returns the default MIDI controller assigned to this parameter, default: -1 (none). */
  virtual int getDefaultMidiController() const 
  { ScopedPointerLock spl(mutex); return defaultMidiController; }

  /** Returns the MIDI controller assigned to this parameter, default: -1 (none). */
  virtual int getAssignedMidiController() const 
  { ScopedPointerLock spl(mutex); return assignedMidiController; }

  /** Informs whether or not the assigned MIDI controller is the default setting. */
  virtual bool isAssignedControllerDefault() const
  {
    ScopedPointerLock spl(mutex); return (assignedMidiController == defaultMidiController);
  }

/** Converts the current automation value into a MIDI control value and returns it. */
  virtual int getMidiControllerValue() const
  {
    ScopedPointerLock spl(mutex); return (int)round(127.0*automationValue);
  }

  /** Informs whether or not this parameter is in default state, that is: are the assigned 
  lower and upper limits and the assigned Controller still in their default settings. */
  virtual bool isInDefaultState() const;

protected:

  /** Restrics a value to the permitted parameter range as defined by minValue and
  maxValue. */
  //virtual double restrictValueToParameterRange(double valueToRestrict);

  /** Restrics a value to the permitted automatable range as defined by lowerAutomationLimit and 
  upperAutomationLimit. */
  virtual double restrictValueToAutomatableRange(double valueToRestrict);

  /** Restricts the value, defaultValue, and automation limits to the range of the parameter, 
  makes sure that max = min, makes sure that values > 0 for exponential scaling, etc.. */
  virtual void valueSanityCheck();

  /** Updates the current value of the parameter from a (possibly new) automation value. */
  virtual void updateParameterFromAutomation();

  /** Updates the current automaztion-value from the parameter itself (which may have a new value). */
  virtual void updateAutomationFromParameter();

  double automationValue;        // automation setting, in the range 0...1
  double lowerAutomationLimit;   // lower limit for automation
  double upperAutomationLimit;   // upper limit for automation
  int    defaultMidiController;  // default MIDI controller
  int    assignedMidiController; // currently assigned MIDI controller
  bool   learnMode;              // flag, to indicate midi learn mode - in learn mode, the 
                                 // parameter will use the next incoming cotroller-event to 
                                 // re-assign the midi controller

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableParameter)
};



//=================================================================================================


class MetaParameterManager;


/** A subcclass of Parameter that can be controlled via a MetaParameter. To do so, it maintains a 
pointer to a MetaParameterManager, where it can attach itself to one of the managed MetaParameters
there. To make this work, you will have to call setMetaParameterManager to pass in the 
MetaParameterManager to use. once this is done, you can attach this parameter to one of the 
MetaParameters by calling attachToMetaParameter

\todo maybe factor the handling of proportionalValue out into parameter baseclass
*/

class JUCE_API MetaControlledParameter : public Parameter
{

public:

  /** Constructor */
  MetaControlledParameter(const juce::String& name, double min = 0.0, double max = 1.0, 
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);

  /** Sets the value of the parameter where the input argument is assumed to be normalized to the 
  range 0...1 - arguments in this range will be mapped to the range between lowerAutomationLimit 
  and upperAutomationLimit and the actual parameter will be adjusted accordingly. */
  virtual void setProportionalValue(double newProportionalValue, bool sendNotification, 
    bool callCallbacks);

  /** Converts the clear text value to a proportional value in the range 0..1 according to our 
  scaling/mapping function. */
  virtual double valueToProportion(double value);

  /** Converts a proportional value in the range 0..1 to a clear text value according to our 
  scaling/mapping function. */
  virtual double proportionToValue(double proportion);

  /** Sets up the MetaParameterManager to use. */
  virtual void setMetaParameterManager(MetaParameterManager *newManager);

  /** Attaches this parameter to the MetaParameter with the given index (in the 
  MetaParameterManager). */
  virtual void attachToMetaParameter(int metaParameterIndex);

  /** Returns the normalized value in the range 0..1. */
  inline double getProportionalValue() { return proportionalValue; }

protected:

  double proportionalValue;  // maybe factor out to Parameter

  int metaIndex = -1;
  MetaParameterManager* metaParaManager = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaControlledParameter)
};

//=================================================================================================

/** A class to represent meta-parameters, i.e. parameters that control other (lower level) 
parameters. It derives from ParameterObserver in order to also provide a means of cross-coupling
between the dependent parameters - whenever one of them changes, we get notified here and also
set up all other dependent parameters accordingly 
\todo: maybe factor out coupling functionality into subclass CouplingMetaParameter */

class JUCE_API MetaParameter : public ParameterObserver
{

public:

  MetaParameter();

  /** Sets this MetaParameter to the given value and updates all dependent MetaControlledParameters
  accordingly. MetaParameter values are always in the normalized range 0..1. */
  void setValue(double newValue);

  /** Attaches the given MetaControlledParameter to this MetaParameter. */
  void attachParameter(MetaControlledParameter* p);

  /** Detaches the given MetaControlledParameter from this MetaParameter. */
  void detachParameter(MetaControlledParameter* p);

  // callbacks:
  virtual void parameterChanged(Parameter* p) override;
  virtual void parameterIsGoingToBeDeleted(Parameter* p) override;


protected:

  double value = 0.0;

  std::vector<MetaControlledParameter*> params; // list of pointers to the dependent parameters

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameter)
};

//=================================================================================================

/** A class to manage a bunch of MetaParameters, allowing MetaControlledParameter objects to attach
themselves to any of our manages MetaParameters. */

class JUCE_API MetaParameterManager /*: public ParameterObserver*/
{

public:

  /** Attaches the passed MetaControlledParameter to the MetaParameter with given index and 
  returns if this was successful (it may fail, if you pass an out-of-range index). */
  virtual bool attachParameter(MetaControlledParameter* param, int index);

protected:

  std::vector<MetaParameter*> metaParams;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameterManager)
};


/* 
Idea: 3 classes

MetaControlledParameter: 
-has a pointer to a MetaParameterManager where it can register itself to listen to one of the
 MetaParameters

MetaParameter:
-maintains a list of dependent MetaControlledParameters and sets all of of them when its 
 setAutomationValue gets called
-listens to each of its dependent parameters to provide cross-coupling: when one of the dependent
 parameters changes, all others should change as well

MetaParameterManager:
-maintains a list of MetaParameters to allow MetaControlledParameter objects set themselves up
 to listen to one specific MetaParameter


todo:
-don't use class AutomatableParameter anymore
-instead, make it possible to assign MIDI-controllers to a MetaParameter (maybe through a subclass
 MidiControlledMetaParameter)

*/

#endif 
