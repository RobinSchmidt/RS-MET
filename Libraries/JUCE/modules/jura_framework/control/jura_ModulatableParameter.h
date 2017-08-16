#ifndef jura_ModulatableParameter_h
#define jura_ModulatableParameter_h

/*

Idea to extend the same concept of MetaParameters to a modulation system for synthesizers:

ModulationSource:
-provide interface to pull out a current value
-subclasses can be fixed values (like gui parameter sliders), envelope generators, low frequency 
 oscillators, step sequencers, midi-note values, etc.
-At each sample, the first thing to do is to compute all current values of all existing modulation 
 sources. This should be done in a central place in some higher level object (maybe some subclass 
 of AudioModule)

ModulationTarget:
-is typically some (lower level) dsp algorithm parameter such as a cutoff frequency of a filter
-keeps a list of assigned sources
-keeps a pointer to ModulationManager where it can sign up to receive inputs from sources
-when a sample is produced, pulls out values of all its connected sources and combines them and 
 then sets up the algorithm parameter accordingly
-perhaps combination can additively for some sources and multiplicatively for others: 1st add up 
 all additive sources, then multiply in all multiplicative sources - keep 2 lists
-can regulate the amount of the influence of each source
-maybe can apply different response curves to each source
-hmm...maybe 4 the last 2 points, we need an additional class: ModulationConnection

ModulationManager:
-allows ModulationsTargets to de/register themselves to ModulationSources

similarities to MetaParameters:
-needs similar sign-up system on the gui (right-click - assign, etc.)

differences to MetaParameters:
-each modulation target can be connected to any number of modulation sources (many-to-many instead 
 of one-to-many)
-uses pull rather than push mechanism to update the dependent parameter

feedback:
-it is desirable to be able for (parameters of) ModulationsSources themselves be modulated by 
 other ModulationSources or even by their own output
-maybe the parameters of a ModulationSource (like an LFO freq) should be subclasses of 
 ModulationTarget ...or just BE a ModulationTarget object - or maybe we'll need a class
 ModulatableParameter


maybe make the class hierarchy like this:
Parameter <- MetaControlledParameter <- ModulatableParameter <- PolyphonicParameter

*/


// forward declarations:
class ModulationManager; 
class ModulationSource;
class ModulationTarget;


//=================================================================================================

/** Baseclass for participants in the modulation system. Subclasses are ModulationSource and 
ModulationTarget and this class factors out what they have in common, which is mainly the access 
functions to the ModulationManager. */

class JUCE_API ModulationParticipant
{

public:

  /** Constructor */
  ModulationParticipant(ModulationManager* managerToUse = nullptr) 
  {
    modManager = managerToUse;
  }

  /** Destructor */
  virtual ~ModulationParticipant() {}
    // maybe in the destructor, we should try to dynamically cast ourselves to ModulationTarget
    // and ModulationSource and if any of these two casts is successfull, de-register ourselves
    // from our modManager to avoid dangling pointers there. ...or maybe do it in the two  subclass
    // destructors to avoid the cast

  /** Sets up the ModulationManager that should be used for registering ModulationSources and 
  ModulationTargets. Should be called sometime soon after construction. */
  void setModulationManager(ModulationManager* managerToUse)
  {
    modManager = managerToUse;
  }
  // maybe get rid of this and demand it to be passed to the constructor

  /** Returns a pointer to the ModulationManager object that is used for registering 
  ModulationSources and ModulationTargets. */
  ModulationManager* getModulationManager()
  {
    return modManager;
  }

  /** Registers the given ModulationSource to make it available to ModulationTargets. */
  void registerModulationSource(ModulationSource* source);

  /** De-registers a ModulationSource. */
  void deRegisterModulationSource(ModulationSource* source);

  /** Registers the given ModulationTarget. */
  void registerModulationTarget(ModulationTarget* target);

  /** De-registers a ModulationTarget. */
  void deRegisterModulationTarget(ModulationTarget* target);

  /** Returns a reference to the array of available ModulationSource objects. */
  const std::vector<ModulationSource*>& getAvailableModulationSources();

  /** Returns a reference to the array of available ModulationTarget objects. */
  const std::vector<ModulationTarget*>& getAvailableModulationTargets();


protected:

  ModulationManager* modManager = nullptr;


  // these are empty dummy arrays to which references will be returned by 
  // getAvailableModulationSources/Targets in case our modManager is a nullptr (this is somehow
  // ugly design, but however):
  static std::vector<ModulationSource*> dummySources;
  static std::vector<ModulationTarget*> dummyTargets;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationParticipant)
};


//=================================================================================================

/** Baseclass for modulation sources. */

class JUCE_API ModulationSource : public ModulationParticipant
{

public:

  /** Constructor */
  ModulationSource(ModulationManager* managerToUse = nullptr) 
    : ModulationParticipant(managerToUse) {}

  /** Destructor */
  virtual ~ModulationSource();

  /** Returns a pointer to the modulation value. Modulation targets should retrieve this pointer 
  when they are connected to */
  double* getModulationValuePointer() { return &modValue; }

  /** Should be overriden by subclasses to update the "modValue" member variable per sample. Once 
  it is updated, connected modulation targets can access it using the pointer-to-double that they 
  have previously retrieved via getValuePointer. */
  virtual void updateModulationValue() = 0;

  /** Adds a modulation target to our list of attached targets. We keep this list here mainly to 
  detach the target in the case, the source gets deleted. */
  void addModulationTarget(ModulationTarget* target)
  {
    appendIfNotAlreadyThere(targets, target);
  }

  /** Removes a modulation target from our list of attached targets. */
  void removeModulationTarget(ModulationTarget* target)
  {
    removeFirstOccurrence(targets, target);
  }


  //juce::String getModulationSourceName() = 0;

protected:

  double modValue = 0;

  std::vector<ModulationTarget*> targets; // do we need this?

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationSource)
};

//=================================================================================================

/** Baseclass for modulation targets.  */

class JUCE_API ModulationTarget : public ModulationParticipant
{

public:

  /** Constructor */
  ModulationTarget(ModulationManager* managerToUse = nullptr) 
    : ModulationParticipant(managerToUse) {}

  /** Destructor */
  virtual ~ModulationTarget() {}

  /** Sets the nominal, unmodulated value. This will be used as reference, when a modulated value 
  will be computed. */
  void setUnmodulatedValue(double newValue)
  {
    unmodulatedValue = newValue;
  }

  /** Adds a ModulationSource to this ModulationTarget with an optional modulation amount. */
  void addModulationSource(ModulationSource* source, double amount = 0)
  {
    source->addModulationTarget(this);
    appendIfNotAlreadyThere(sources,      source);
    appendIfNotAlreadyThere(sourceValues, source->getModulationValuePointer());
  }

  /** Removes the ModulationSource with given index from this ModulationTarget. */
  void removeModulationSource(int index)
  {
    sources[index]->removeModulationTarget(this);
    remove(sources,      index);
    remove(sourceValues, index);
  }

  /** Removes the given ModulationSource from this ModulationTarget. */
  void removeModulationSource(ModulationSource* source)
  {
    int index = find(sources, source);
    removeModulationSource(index);
  }

  /** Sets the amount by which the ModulationSource with given index modulates this 
  ModulationTarget. */
  void setModulationAmount(int sourceIndex, double newAmount)
  {
    amounts[sourceIndex] = newAmount;
  }

  /** This function should be called from some central place in outside code to compute/update a 
  modulated value per sample before the value is used. It starts with the unmodulated value and 
  applies all attached ModulationSources to it (with their appropriate amounts) and stores the 
  result in modulatedValue. */
  void computeModulatedValue()
  {
    double tmp = unmodulatedValue;


    // maybe this block can be optimized out of this function, scaler can be made a member and 
    // assigned elsewhere (in a function that is not called per sample)
    bool relative = false;  // make user adjustable member variable (switch between absolute and relative modulation)
    double scaler = 1;      // maybe this too
    if(relative)
      scaler = unmodulatedValue;


    for(int i = 0; i < size(sourceValues); i++)
      tmp += amounts[i] * (*sourceValues[i]) * scaler;
    modulatedValue = tmp;
  }
  // hmm...maybe it would be more efficient, if we do the gathering of the values in 
  // ModulationManager which would the keep an array of some kind of ModulationConnection. This 
  // way, we could perhaps avoid iterating through a lot of zero-sized sourceValues arrays



protected:

  double unmodulatedValue = 0;
  double modulatedValue = 0;

  std::vector<ModulationSource*> sources;
  std::vector<double*> sourceValues;
  std::vector<double>  amounts;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationTarget)
};


//=================================================================================================

/** This class is responsible for keeping a list of all available ModulationSources and making it 
possible for ModulationTargets to register themselves to any number of these sources. */

class JUCE_API ModulationManager
{

public:

  /** Constructor */
  ModulationManager() {}

  /** Destructor */
  virtual ~ModulationManager() {}

  /** Registers the given ModulationSource to make it available to ModulationTargets. */
  void registerModulationSource(ModulationSource* source)
  {
    appendIfNotAlreadyThere(modulationSources, source);
    source->setModulationManager(this);
  }

  /** De-registers a ModulationSource. */
  void deRegisterModulationSource(ModulationSource* source)
  {
    removeFirstOccurrence(modulationSources, source);
    source->setModulationManager(nullptr); // maybe we should do this conditionally when the passed
                                           // source is actually in the array
  }

  /** Registers the given ModulationTarget. */
  void registerModulationTarget(ModulationTarget* target)
  {
    appendIfNotAlreadyThere(modulationTargets, target);
    target->setModulationManager(this);
  }

  /** De-registers a ModulationTarget. */
  void deRegisterModulationTarget(ModulationTarget* target)
  {
    removeFirstOccurrence(modulationTargets, target);
    target->setModulationManager(nullptr); // maybe we should do this conditionally when the passed
                                           // target is actually in the array
  }

  /** Returns a reference to our list of available ModulationSources. */
  const std::vector<ModulationSource*>& getAvailableModulationSources() { return modulationSources; }

  /** Returns a reference to our list of available ModulationTargets. */
  const std::vector<ModulationTarget*>& getAvailableModulationTargets() { return modulationTargets; }


  void applyModulations()
  {
    // todo: this should be the method that should be called once per sample to compute all the
    // outputs of the ModulationSources and apply them to their attached ModulationTargets. We can
    // then have a subclass of AudioModule which is also a subclass of ModulationManager and 
    // override the per-sample callback in order to first call applyModulations and after that, 
    // proceed with the regular per-sample computations as usual.
    //...sooo we can perhaps get rid of ModulationTarget::computeModulatedValue

    int i;

    for(i = 0; i < size(modulationSources); i++)
      modulationSources[i]->updateModulationValue();

    for(i = 0; i < size(modulationTargets); i++)
      modulationTargets[i]->computeModulatedValue();
      // this loop should be replaced by a loop over ModulationConnections because typically it
      // will be wasteful to loop over all possible targets (i.e. parameters, many of which may
      // not be modulated at all)
  }

protected:

  std::vector<ModulationSource*> modulationSources;  // array of the available sources
  std::vector<ModulationTarget*> modulationTargets;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationManager)
};

//=================================================================================================

/** A subclass of Parameter that is suitable as target for modulations by also being a subclass of
ModulationTarget. */

class JUCE_API ModulatableParameter : public MetaControlledParameter, public ModulationTarget
{

public:

  /** Constructor */
  ModulatableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0)
    : MetaControlledParameter(name, min, max, defaultValue, scaling, interval) {}


  virtual void setValue(double newValue, bool sendNotification, bool callCallbacks)
  {
    MetaControlledParameter::setValue(newValue, sendNotification, callCallbacks);
    ModulationTarget::setUnmodulatedValue(newValue);
  }


  /** We suppose that modulatedValue of the ModulationTarget baseclass has already been computed, 
  i.e. ModulationTarget::computeModulatedValue has been called, such that modulatedValue has a 
  legitimate value with all modulations applied. Here we pull out this modulated value and call our 
  valueChangeCallback with it. The function is supposed to be called per sample for each modulated 
  Parameter. */
  void callCallbackWithModulatedValue()
  {
    if( valueChangeCallbackDouble != nullptr )
      valueChangeCallbackDouble->call(modulatedValue);
  }

protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableParameter)
};





// stuff below may not be needed - we'll see:

//=================================================================================================

/**   */

class JUCE_API ParameterModulator : public ModulationSource
{

public:

  /** Constructor */
  ParameterModulator() {}

  /** Destructor */
  virtual ~ParameterModulator() {}



protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterModulator)
};

//=================================================================================================

/**  */

class JUCE_API ParameterModulationManager
{

public:

  ParameterModulationManager() {};


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterModulationManager)
};

#endif
