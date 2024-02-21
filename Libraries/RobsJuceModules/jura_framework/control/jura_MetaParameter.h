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


 meta-mapping use case:
 -cutoff and reso are attached to the same meta
 -reso has nonmonotonic meta-mapping, say, with a bump in the middle (i.e. high-reso at mid-cutoff)
 -user moves reso-slider
 -reso-slider calls setValue on the reso-parameter
 -reso-parameter notifies meta about the value-change
 -meta sets up cutoff
 -how does the the meta know which is the right cutoff? we can't map back from the actual 
  reso-value - somehow we must keep track of reso-slider's normalized value
 -the reso-slider's position should reflect that unmapped narmalized value but the numeric readout
  should reflect the actual resnance
 ->the presence of the meta-attachment would change the behavior (mapping) of the reso-slider
   ...maybe we must override the mapping function in MetaControlledParameter

*/

//=================================================================================================

/** This class is used to map the normalized parameter range 0..1 that is used by meta-parameters
to itself via a user definable mapping function, such that several parameters that are controlled
via the same MetaParameter can map this input-range of the meta parameter to their own normalized
range in arbitrary ways. */

class JUCE_API rsMetaParameterMapper : public RAPT::rsNodeBasedFunction<double>
{

public:

  rsMetaParameterMapper();
  
  virtual ~rsMetaParameterMapper() {}

  // overriden from rsNodeBasedFunction to apply some restrictions:
  //size_t addNode(double x, double y) override;
  //bool removeNode(size_t index) override;
  size_t moveNode(size_t index, double newX, double newY) override;

  virtual bool isNodeRemovable(size_t index) const override;
  virtual size_t constrainNode(size_t index) override;


  /** Maps an incoming normalized parameter in the range 0..1 to mapped 0..1 range using the 
  function defined by the nodes in our rsNodeBasedFunction baseclass. */
  double map(double x) { return rsNodeBasedFunction<double>::getValue(x); } // rename to applyFunction

  /** Tries to find a normalized x-value that corresponds to the given normalized y-value... */
  double unmap(double y) { return rsNodeBasedFunction<double>::applyInverseFunction(y); } 

  /** Returns true, if this map is the default identity map with two nodes at (0,0) and (1,1). */
  bool isDefaultMap() const;
    // maybe factor out into a baseclass that handles all kinds of maps and
    // make virtual, so subclasses can define their own default maps

  /** Initializes this maps to its default state which is the identity map. */
  void initToIdentity();

  /** Initializes this maps to a flat line at the given value. */
  void initToFlat(double value);

  /** Returns the state of this mapper as xml element with given tag name. */
  XmlElement* getStateAsXml(const juce::String& tagName) const;

  /** Sets up the state of this mapper from the given XmlElement. */
  void setStateFromXml(const XmlElement& xmlState);


  static String shapeIndexToString(int index);

  static int stringToShapeIndex(const String& shapeString);

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsMetaParameterMapper)
};

//=================================================================================================

class MetaParameter;
class MetaParameterManager;

/** A subclass of Parameter that can be controlled via a MetaParameter. To do so, it maintains a
pointer to a MetaParameterManager, where it can attach itself to one of the managed MetaParameters
there. To make this work, you will have to call setMetaParameterManager to pass in the
MetaParameterManager to use. Once this is done, you can attach this parameter to one of the
MetaParameters by calling attachToMetaParameter. */

class JUCE_API MetaControlledParameter : public rsSmoothableParameter
{

public:

  /** Constructor */
  MetaControlledParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, Parameter::Scaling scaling = Parameter::Scaling::LINEAR, 
    double interval = 0.0);

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets this parameter up according to a given MetaParameter value and optionally notifies
  observers and/or calls the callbacks. */
  virtual void setFromMetaValue(double newMetaValue, bool sendNotification,
    bool callCallbacks);

  /** Overriden to find the normalized value that corresponds to the passed desired value, taking
  both mapping functions into account. If there are several normalized values possible that map to
  the desired value, the first one will be used. If there is no possible value, ...  */
  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks) override;

  /** Sets up the MetaParameterManager to use. This function should be called once shortly after
  this MetaControlledParameter object has been created and the passed manager object should remain
  valid for the whole lifetime of this object. */
  virtual void setMetaParameterManager(MetaParameterManager* newManager);

  /** Attaches this parameter to the MetaParameter with the given index (in the
  MetaParameterManager). The flatMap parameter decides, whether the mapping function should 
  initially be a flat line at the current value (otherwise it will be left as is - which initially 
  is the identity map). */
  virtual void attachToMetaParameter(int metaParameterIndex, bool flatMap = false);

  /** Detaches this parameter from any MetaParameter, it may be attched to. */
  virtual void detachFromMetaParameter();

  /** Overriden to apply the additional meta-mapper. */
  virtual void setNormalizedValue(double newValue, bool sendNotification, 
    bool callCallbacks) override;

  /** This is called from rsMetaMapEditor whenever the mapping curve changes in oder to update the
  underlying dsp parameter and the numerical slider readout. It's not very elegant to require
  rsMetaMapEditor to call this, but it's good enough for the moment. */
  void metaMapChanged();

  /** Initializes the meta map to a flat line at the current value. */
  void initMetaMapToFlat();

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Return a pointer to the MetaParameterManager object that is used here. */
  inline MetaParameterManager* getMetaParameterManager() { return metaParaManager; }

  /** Returns the index of the MetaParameter that this parameter is attached to. If it's not
  attached to any MetaParameter, it returns -1. */
  inline int getMetaParameterIndex() const { return metaIndex; }

  /** Returns true, if this Parameter has an attached MetaParameter, false otherwise. */
  inline bool hasAttachedMeta() const { return metaIndex != -1; }

  /** Returns a pointer to the attached meta-parameter (nullptr, if none). */
  MetaParameter* getAttachedMeta();

  /** Returns the name of the MetaParameter which this Parameter is attached to. */
  String getMetaParameterName(); // may be const?

  /** Returns a pointer to our meta-mapper object. */
  rsMetaParameterMapper* getMetaMapper() { return &metaMapper; }

  //virtual double getNormalizedDefaultValue() const override
  //{ ScopedPointerLock spl(mutex); return 0.5; } // preliminary

  //-----------------------------------------------------------------------------------------------
  // \name State recall

  /** Overriden to possibly store the mapping function, if necessarry. */
  virtual void saveToXml(XmlElement* xml) const override;

  /** Overriden to possibly recall the mapping function, if necessarry. */
  virtual void recallFromXml(const XmlElement& xml) override;

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Applies first the user defined node-based map to the normalized input value (mapping the 
  range 0..1 onto itself) and then the inherited map that maps the normalized 0..1 range to the
  actual parameter range (for example, exponentially from 20..20000 for a cutoff frequency 
  parameter) and returns the resulting value. */
  inline double applyBothMaps(double x) { return mapper->map(metaMapper.map(x)); }

protected:

  virtual void setNormalizedTargetValue(double newTargetValue, bool sendNotification, 
    bool callCallbacks) override;

  int metaIndex = -1;
  MetaParameterManager* metaParaManager = nullptr; // use a Null Object by default
  rsMetaParameterMapper metaMapper;

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

  /** Attaches the given MetaControlledParameter to this MetaParameter. The flatMap parameter decides, whether the mapping function should 
  initially be a flat line at the current value. */
  void attachParameter(MetaControlledParameter* p, bool flatMap = false);

  /** Detaches the given MetaControlledParameter from this MetaParameter, if it is attached. It 
  returns true, if it was attached (i.e. an actualy removal took place). If it wasn't attached, 
  nothing happens and false will be returned. */
  bool detachParameter(MetaControlledParameter* p);

  /** Resets this MetaParameter to its default value of 0.5 (causing callbacks and
  notifications). */
  inline void resetToDefaultValue() { setMetaValue(0.5); }

  /** Gives this MetaParameter a name. */
  void setName(const String& newName) { name = newName; }

  /** Returns the name of this MetaParameter. */
  String getName() const { return name; }

  // callbacks:
  virtual void parameterChanged(Parameter* p) override;
  virtual void parameterWillBeDeleted(Parameter* p) override;


protected:

  /** The value of the meta parameter. It is initialized to 0.5 which is used as the default value.
  The center value seems a reasonable choice for most continuous parameters (which we assume to be
  the majority). Note that changing this may break presets and states because meta values are only
  stored when they differ from the default value - so don't change that. */
  double metaValue = 0.5; // NEVER change the 0.5 initialization, state save/recall relies on that

  /** Name of this MetaParameter. */
  juce::String name;

  std::vector<MetaControlledParameter*> params; // list of pointers to the dependent parameters

  friend class MetaParameterManager;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameter)
};

//=================================================================================================

class JUCE_API MetaParameterManager;

/** Baseclass for objects that need to observe a MetaParameterManager. */

class JUCE_API MetaParameterManagerObserver
{

public:

  MetaParameterManagerObserver() = default;
  virtual ~MetaParameterManagerObserver() = default;

  /** Function that gets called when one of the meta-names hase changed. */
  virtual void metaNameChanged(MetaParameterManager* manager, int index) = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameterManagerObserver)
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
  bool attachParameter(MetaControlledParameter* param, int metaIndex, bool flatMap = false);

  /** If the passed parameter is attached to any of our managed MetaParameters, this function
  will detach it (otherwise it will have no effect). */
  void detachParameter(MetaControlledParameter* param);

  /** Adds a new observer to our list of observers. */
  void registerObserver(MetaParameterManagerObserver* observer)
  {
    appendIfNotAlreadyThere(metaObservers, observer);
  }

  /** Removes the given observer from our list of observers. */
  void deRegisterObserver(MetaParameterManagerObserver* observer)
  {
    removeFirstOccurrence(metaObservers, observer);
  }

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

  /** Updates the name of the meta-parameter with givne index, taking into account the depedent
  parameters and notifies our observers about that change. */
  void updateMetaName(int index);


protected:

  std::vector<MetaParameter*> metaParams;

  std::vector<MetaParameterManagerObserver*> metaObservers;

  bool autoUpdateMetaNames = true;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MetaParameterManager)
};

#endif
