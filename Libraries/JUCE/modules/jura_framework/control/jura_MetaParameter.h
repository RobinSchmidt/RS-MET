#ifndef jura_MetaParameter_h
#define jura_MetaParameter_h

/*
The meta-parameter handling involves 3 classes:

MetaControlledParameter:
-subclass of Parameter that can be controlled by a MetaParameter
-has a pointer to a MetaParameterManager where it can register itself to listen to one of the
 MetaParameters (which are manitained as a list there)
-todo: provide custom mapping curves (from input 0..1 to output 0..1 - before it goes into the
 regular mapping provided by the Parameter baseclass, so we have a two level mapping)
-todo: provide smoothing (maybe by an intermediate subclass SmoothedParameter:
 i.e. Parameter <- SmoothedParameter <- MetaControlledParameter)...or maybe just integrate it
 directly .. and at some stage, we may need non-smoothed meta-controlled or non-meta-controlled
 smoothed parameters - then we can factor out the more basic class 
 or: use the decorator pattern: make an abstract ParameterBase baseclass with an interface like
 getValue, setValue, etc - then Parameter, MetaControlledParameter and SmoothedParameter can all 
 derive from ParameterBase, and SmoothedParameter and MetaControlledParameter do not derive from 
 Parameter but instead keep a reference to a ParameterBase object and forward requests to it....or 
 something
  
MetaParameter:
-maintains a list of dependent MetaControlledParameters and updates all of them when its
 setMetaValue method gets called
-watches each of its dependent MetaControlledParameters for changes (by means of being subclass of
 ParameterObserver) to provide cross-coupling: when one of the dependent parameters changes, all
 others are updated as well

MetaParameterManager:
-maintains a list of MetaParameters to allow MetaControlledParameter objects register themselves
 to listen to one specific MetaParameter

todo:
-make it possible to assign MIDI-controllers to a MetaParameter (maybe through a subclass
 MidiControlledMetaParameter)
-don't use class AutomatableParameter anymore - consider it deprecated, but maybe it should still
 be renamed into MidiControlledParameter...and perhaps be moved inot the jura_processors module
 we may still need it for legacy code compatibility
*/

/** This is an abstract baseclass for classes that are supposed to map an input value (which is 
normalized to the range 0..1) to an output value in the same range. It may (optionally) be used
in a MetaControlledParameter to realize arbitrary mappings from the MetaParameter to the target 
Parameter. To this end, you would have to define a subclass of NormalizedParameterMapper and 
configure the respective MetaControlledParameter object with an object of your subclass
(using MetaControlledParameter::setParameterMapper). */

class JUCE_API NormalizedParameterMapper
{
public:
  virtual ~NormalizedParameterMapper(){};

  /** Override this function in your subclass to map a normalized input value (in the range 0..1) 
  to the corresponding output value (in the same range). You may map different input values to the 
  same output value, i.e. your mapping function does not need to be injective/invertible, but it 
  should have the property that any value between 0..1 (ends inclusive) is a valid input and the 
  output values should be restricted to that same range. */
  virtual double map(double input) = 0;

};

//=================================================================================================

class MetaParameterManager;

/** A subclass of Parameter that can be controlled via a MetaParameter. To do so, it maintains a
pointer to a MetaParameterManager, where it can attach itself to one of the managed MetaParameters
there. To make this work, you will have to call setMetaParameterManager to pass in the
MetaParameterManager to use. Once this is done, you can attach this parameter to one of the
MetaParameters by calling attachToMetaParameter. */

class JUCE_API MetaControlledParameter : public Parameter
{

public:

  /** Constructor */
  MetaControlledParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0);

  /** Sets this parameter up according to a given MetaParameter value and optionally notifies
  observers and/or calls the callbacks. */
  virtual void setFromMetaValue(double newMetaValue, bool sendNotification,
    bool callCallbacks);

  /** Sets up the MetaParameterManager to use. This function should be called once shortly after
  this MetaControlledParameter object has been created and the passed manager object should remain
  valid for the whole lifetime of this object. */
  virtual void setMetaParameterManager(MetaParameterManager* newManager);

  /** Sets the mapper object that maps between the MetaParameter value at the input side and the
  proportional value of the target parameter on the output side. Such a mapper is optional, if you 
  pass none, an identity mapping will be used by default. The mapper object should be valid for the
  entire lifetime of this MetaControlledParameter. To reset it, you can use this function with a 
  nullptr argument. */
  virtual void setParameterMapper(NormalizedParameterMapper* newMapper);
    // feature is not yet tested

  /** Attaches this parameter to the MetaParameter with the given index (in the
  MetaParameterManager). */
  virtual void attachToMetaParameter(int metaParameterIndex);

  /** Detaches this parameter from any MetaParameter, it may be attched to. */
  virtual void detachFromMetaParameter();

  /** Return a pointer to the MetaParameterManager object that is used here. */
  inline MetaParameterManager* getMetaParameterManager() { return metaParaManager; }

  /** Returns the index of the MetaParameter that this parameter is attached to. If it's not
  attached to any MetaParameter, it returns -1. */
  inline int getMetaParameterIndex() { return metaIndex; }

  /** Returns the name of the MetaParameter which this Parameter is attached to. */
  String getMetaParameterName();

protected:

  int metaIndex = -1;
  MetaParameterManager* metaParaManager = nullptr; // use a Null Object by default
  NormalizedParameterMapper* mapper = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaControlledParameter)
};

//=================================================================================================

/** A class to represent meta-parameters, i.e. parameters that control other (lower level)
parameters. It derives from ParameterObserver in order to also provide a means of cross-coupling
between the dependent parameters - whenever one of them changes, we get notified here and also
update all other dependent parameters accordingly.
\todo: maybe factor out coupling functionality into subclass CrossCouplingMetaParameter */

class JUCE_API MetaParameter : public ParameterObserver
{

public:

  MetaParameter();

  /** Sets this MetaParameter to the given value and updates all dependent MetaControlledParameters
  accordingly. MetaParameter values are always in the normalized range 0..1. */
  void setMetaValue(double newValue);

  /** Returns the current (meta) value of the MetaParameter. */
  inline double getMetaValue() { return metaValue; }

  /** Attaches the given MetaControlledParameter to this MetaParameter. */
  void attachParameter(MetaControlledParameter* p);

  /** Detaches the given MetaControlledParameter from this MetaParameter. */
  void detachParameter(MetaControlledParameter* p);

  /** Resets this MetaParameter to its default value of 0.5 (causing callbacks and
  notifications). */
  inline void resetToDefaultValue() { setMetaValue(0.5); }

  /** Gives this MetaParameter a name. */
  void setName(const String& newName) { name = newName; }

  /** Returns the name of this MetaParameter. */
  String getName() const { return name; }

  // callbacks:
  virtual void parameterChanged(Parameter* p) override;
  virtual void parameterIsGoingToBeDeleted(Parameter* p) override;


protected:

  /** The value of the meta parameter. It is initialized to 0.5 which is used as the default value.
  The center value seems a reasonable choice for most continuous parameters (which we assume to be
  the majority). Note that changing this may break presets and states because meta values are only
  stored when they differ from the default value - so don't change that. */
  double metaValue = 0.5; // NEVER change the 0.5 initialization, state save/recall relies on that

  /** Name of this MetaParameter. */
  juce::String name;

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

  /** Returns the number of MetaParameters that are managed by this object. */
  inline int getNumMetaParameters() { return (int) size(metaParams); }

  /** Returns a pointer to the MetaParameter with given index. If the index is out of range, it
  will be a nullptr. */
  MetaParameter* getMetaParameter(int index);

  /** Returns the name of the MetaParameter with given index (empty string, if index is out of
  range). */
  String getMetaParameterName(int index);

  /** Resets all MetaParameters in our array to their default value of 0.5. */
  void resetAllToDefaults();

  /** Tries to set the MetaParameter with given index to the passed newValue and reports if
  this was successful (it will fail, if index is out of range). */
  bool setMetaValue(int index, double newValue);

  /** Tries to give a new name to the MetaParameter with given index and reports if this was
  successful (it will fail, if index is out of range). */
  bool setMetaName(int index, const String& newName);



protected:

  std::vector<MetaParameter*> metaParams;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameterManager)
};

#endif
