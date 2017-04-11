#ifndef jura_MetaParameter_h
#define jura_MetaParameter_h

/* 
Idea: 3 classes

MetaControlledParameter: 
-subclass of Parameter that can be controlled by a MetaParameter
-has a pointer to a MetaParameterManager where it can register itself to listen to one of the
MetaParameters
-todo: provide custom mapping curves (from input 0..1 to output 0..1 - before it goes into the
regular mapping provided by the Parameter baseclass)
-todo: provide Smoothing (maybe by an intermediate subclass SmoothedParameter)

MetaParameter:
-maintains a list of dependent MetaControlledParameters and sets all of of them when its 
setValue method gets called
-watches each of its dependent parameters for changes (by means of being subclass of 
ParameterObserver) to provide cross-coupling: when one of the dependent parameters changes, all 
others are updated as well

MetaParameterManager:
-maintains a list of MetaParameters to allow MetaControlledParameter objects set themselves up
to listen to one specific MetaParameter

todo:
-make it possible to assign MIDI-controllers to a MetaParameter (maybe through a subclass
MidiControlledMetaParameter)
-don't use class AutomatableParameter anymore
*/

class MetaParameterManager;

/** A subclass of Parameter that can be controlled via a MetaParameter. To do so, it maintains a 
pointer to a MetaParameterManager, where it can attach itself to one of the managed MetaParameters
there. To make this work, you will have to call setMetaParameterManager to pass in the 
MetaParameterManager to use. Once this is done, you can attach this parameter to one of the 
MetaParameters by calling attachToMetaParameter.

\todo 
-maybe factor the handling of proportionalValue out into parameter baseclass
-maybe we should here have a setMetaValue function to replace the setProportionalValue that
 1st maps the meta-range 0..1 arbitrarily to itself (via some kind of ParameterMapper object) and
 the calls Parameter::setProportionalValue with the mapped value. the mapper should be optional
 and default to the identity-mapping.
*/

class JUCE_API MetaControlledParameter : public Parameter
{

public:

  /** Constructor */
  MetaControlledParameter(const juce::String& name, double min = 0.0, double max = 1.0, 
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);


  //--------------------
  // move to Parameter:
  /** Sets the value of the parameter where the input argument is assumed to be normalized to the 
  range 0...1  .... */
  virtual void setProportionalValue(double newProportionalValue, bool sendNotification, 
    bool callCallbacks);

  /** Converts the clear text value to a proportional value in the range 0..1 according to our 
  scaling/mapping function. */
  virtual double valueToProportion(double value);

  /** Converts a proportional value in the range 0..1 to a clear text value according to our 
  scaling/mapping function. */
  virtual double proportionToValue(double proportion);

  /** Returns the normalized value in the range 0..1. */
  inline double getProportionalValue() { return proportionalValue; }
  //----------------------


  /** Sets up the MetaParameterManager to use. */
  virtual void setMetaParameterManager(MetaParameterManager *newManager);

  /** Attaches this parameter to the MetaParameter with the given index (in the 
  MetaParameterManager). */
  virtual void attachToMetaParameter(int metaParameterIndex);

  /** Detaches this parameter from any MetaParameter, it may be attched to. */
  virtual void detachFromMetaParameter();

  /** Returns the index of the MetaParameter that this parameter is attached to. If it's not 
  attached to any MetaParameter, it returns -1. */
  inline int getMetaParameterIndex() { return metaIndex; }

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
\todo: maybe factor out coupling functionality into subclass CrossCouplingMetaParameter */

class JUCE_API MetaParameter : public ParameterObserver
{

public:

  MetaParameter();

  /** Sets this MetaParameter to the given value and updates all dependent MetaControlledParameters
  accordingly. MetaParameter values are always in the normalized range 0..1. */
  void setMetaValue(double newValue);

  /** Attaches the given MetaControlledParameter to this MetaParameter. */
  void attachParameter(MetaControlledParameter* p);

  /** Detaches the given MetaControlledParameter from this MetaParameter. */
  void detachParameter(MetaControlledParameter* p);

  // callbacks:
  virtual void parameterChanged(Parameter* p) override;
  virtual void parameterIsGoingToBeDeleted(Parameter* p) override;


protected:

  double metaValue = 0.0;

  std::vector<MetaControlledParameter*> params; // list of pointers to the dependent parameters

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameter)
};

//=================================================================================================

/** A class to manage a bunch of MetaParameters, allowing MetaControlledParameter objects to attach
themselves to any of our managed MetaParameters. */

class JUCE_API MetaParameterManager
{

public:

  MetaParameterManager() {};

  /** Adds the passed MetaParameter to our list of managed MetaParameters. */
  void addMetaParamater(MetaParameter* metaParameterToAdd);

  /** Attaches the passed MetaControlledParameter to the MetaParameter with given index and 
  returns if this was successful (it may fail, if you pass an out-of-range index). */
  bool attachParameter(MetaControlledParameter* param, int metaIndex);

  /** If the passed parameter is attached to any of our managed MetaParameters, this function
  will detach it (otherwise it will have no effect). */
  void detachParameter(MetaControlledParameter* param);


protected:

  std::vector<MetaParameter*> metaParams;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameterManager)
};

#endif 
