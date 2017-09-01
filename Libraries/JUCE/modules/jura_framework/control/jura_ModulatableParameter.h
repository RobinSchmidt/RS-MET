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

ToDo:
-check ModulationManager::removeConnection - if it also correctly removes the affected-target from 
 the array
-let the user set up the min/max amount values
-let the user select absolute or relative modulation
-let the Amount/Depth parameters be meta-controlled
-optimize, where possible
-test in practice

*/

// forward declarations:
class ModulationManager; 
class ModulationSource;
class ModulationTarget;
class ModulationConnection;

//=================================================================================================

/** Baseclass for participants in the modulation system. Subclasses are ModulationSource and 
ModulationTarget and this class factors out what they have in common, which is mainly the access 
functions to the ModulationManager. There can also be other subclasses that need a convenient 
access point to a ModulationManager, such as an AudioModule subclass that contains modulatable
parameters. */

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

  /** Returns a reference to the array of ModulationConnectionTarget objects. */
  const std::vector<ModulationConnection*>& getModulationConnections();


protected:

  ModulationManager* modManager = nullptr;


  static std::vector<ModulationSource*> dummySources;
  static std::vector<ModulationTarget*> dummyTargets;
  static std::vector<ModulationConnection*> dummyConnections;
  // These are empty dummy arrays to which references will be returned by 
  // getAvailableModulationSources/Targets in case our modManager is a nullptr.  This is somehow
  // ugly design, maybe use the "Null Object" pattern instead:
  // https://sourcemaking.com/design_patterns/null_object)
  // somewhere, we should have a default "null" ModulationManager object lying around to which
  // out pointer is initialized



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

  /** Should be overriden by subclasses to update the "modValue" member variable per sample. It 
  should assign modValue to the output signal value of the modulator. */
  virtual void updateModulationValue() = 0;

  /** Sets up a name for this ModulationSource. This should be unique among all the available 
  ModulationSource names, so it can be used to identify the source in state recall. */
  void setModulationSourceName(const juce::String& newName) { modSourceName = newName; }

  /** Returns the name of this ModulationSource. */
  juce::String getModulationSourceName() { return modSourceName; }

protected:

  double modValue = 0;
  juce::String modSourceName = "ModulationSource";

  friend class ModulationConnection;
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
  virtual ~ModulationTarget();

  /** Must be overriden by subclasses to do whatever they need to do after our modulatedValue has 
  been computed (for example, ModulatableParameter invokes the setter-callback which in turn 
  updates the corresponding value in the core dsp algorithm). */
  virtual void doModulationUpdate() = 0;

  /** Sets the nominal, unmodulated value. This will be used as reference, when a modulated value 
  will be computed. */
  void setUnmodulatedValue(double newValue)
  {
    unmodulatedValue = newValue;
  }

  /** Adds a ModulationSource to this ModulationTarget. The amount of modulation is initially 0. */
  void addModulationSource(ModulationSource* source);

  /** Removes a ModulationSource from this ModulationTarget. */
  void removeModulationSource(ModulationSource* source);

  /** Returns true, if there's a connection between this ModulationTarget and the given 
  ModulationSource. */
  bool isConnectedTo(ModulationSource* source);

  /** Returns a vector of pointers to the ModulationSources which are connected to this 
  ModulationTarget. */
  std::vector<ModulationSource*> getConnectedSources();

  /** Returns a vector of pointers to the available ModulationSources which are not connected to 
  this ModulationTarget. */
  std::vector<ModulationSource*> getDisconnectedSources();

  /** Returns a vector of pointers to ModulationConnections that are incoming into this 
  ModulationTarget. */
  std::vector<ModulationConnection*> getConnections();

  /** Initialized the modulated value by setting it to the unmodulated value. */
  inline void initModulatedValue()
  {
    modulatedValue = unmodulatedValue;
  }

  /** This function must be overriden by subclasses to return a unique name that can be used to 
  identify the target in state recall. */
  virtual juce::String getModulationTargetName() = 0;

protected:

  double unmodulatedValue = 0;
  double modulatedValue = 0;

  friend class ModulationConnection;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationTarget)
};

//=================================================================================================

/** This class represents a connection from a ModulationSource to a ModulationTarget with 
adjustable modulation amount. */

class JUCE_API ModulationConnection
{

public:

  /** Constructor. You should pass a ModulationSource, a ModulationTarget, an amount and a flag to
  to indicate whether this amount is absolute or relative (i.e. scaled by the unmodulated 
  value). */
  ModulationConnection(ModulationSource* source, ModulationTarget* target);

  /** Destructor. */
  virtual ~ModulationConnection(); // maybe make non-virtual

  /** Sets the modulation amount for this connection. */
  void setAmount(double newAmount)
  {
    amount = newAmount;
  }

  /** Sets the parameter to relative mode in which case the output signal of the ModulationSource
  will be multiplied by the ModulationTarget's unmodulated nominal value before being added to the 
  target value. */
  void setRelative(bool shouldBeRelative)
  {
    relative = shouldBeRelative;
  }

  /** Returns the Parameter object that controls the amount of modulation */
  MetaControlledParameter* getAmountParameter() { return amountParam; }

  /** Applies the source-value to the target-value with given amount. */
  inline void apply()
  {
    // optimize away later:
    double scaler = 1; // make this a member
    if(relative)
      scaler = target->unmodulatedValue;

    *targetValue += *sourceValue * amount * scaler; // only this line shall remain after optimization
  }


  /** Returns a xml element containing the information about this connection (needed for state 
  save/recall in ModulationManager) */
  XmlElement* getAsXml();


protected:

  ModulationSource* source;
  ModulationTarget* target;
  double* sourceValue;
  double* targetValue;
  double amount;
  bool relative;

  MetaControlledParameter* amountParam; // maybe it should be a ModulatableParameter? but that may
                                        // raise some issues - maybe later...

  friend class ModulationTarget;
  friend class ModulationManager;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationConnection) 
};

//=================================================================================================

/** This class is responsible for keeping a list of all available ModulationSources and making it 
possible for ModulationTargets to register themselves to any number of these sources. */

class JUCE_API ModulationManager
{

public:

  /** Constructor. You must pass a CriticalSection object that will be used to access our arrays of
  sources, targets and connections (they are accessed from the gui and audio-thread, so the 
  accesses have to mutexed). */
  ModulationManager(CriticalSection *lockToUse);

  /** Destructor */
  virtual ~ModulationManager();

  /** Function to do the per-sample updates of all modulation-sources and targets. Should be called
  from outside code once per sample before the per-sample functions of the actual dsp-algorithms 
  (oscs, filters, whatever) are called. */
  void applyModulations();

  /** Same as applyModulations() but it doesn'T acquire the mutex-lock. This should be used in cases 
  where the caller already has acquired the lock. */
  void applyModulationsNoLock();


  /** \name Connection setup */

  /** Adds a connection between the given source and target. */
  void addConnection(ModulationSource* source, ModulationTarget* target);

  /** Adds the passed ModulationConnection to our array. */
  void addConnection(ModulationConnection* connection);

  /** Removes a connection between the given source and target. */
  void removeConnection(ModulationSource* source, ModulationTarget* target);

  /** Removes all ModulationConnections. */
  void removeAllConnections();

  /** Removes all modulation connections that involve the given source. */
  void removeConnectionsWith(ModulationSource* source);

  /** Removes all modulation connections that involve the given target. */
  void removeConnectionsWith(ModulationTarget* target);


  /** \name Registration of sources and targets */

  /** Registers the given ModulationSource to make it available to ModulationTargets. */
  void registerModulationSource(ModulationSource* source);

  /** De-registers a ModulationSource. */
  void deRegisterModulationSource(ModulationSource* source);

  /** Registers the given ModulationTarget. */
  void registerModulationTarget(ModulationTarget* target); // may not be needed

  /** De-registers a ModulationTarget. */
  void deRegisterModulationTarget(ModulationTarget* target); // may not be needed


  /** \name Inquiry */

  /** Returns true if there's a connection between the given source and target. */
  bool isConnected(ModulationSource* source, ModulationTarget* target);

  /** Returns a reference to our vector of available ModulationSources. */
  const std::vector<ModulationSource*>& getAvailableModulationSources() { return availableSources; }

  /** Returns a reference to our vector of available ModulationTargets. */
  const std::vector<ModulationTarget*>& getAvailableModulationTargets() { return availableTargets; }

  /** Returns a reference to our vector of ModulationConnections. */
  const std::vector<ModulationConnection*>& getModulationConnections()
  { return modulationConnections; }

  /** Given a pointer to a ModulationSource of some type, this function returns the number of 
  ModulationSources of the same type that are registered here. This is used to figure out, for 
  example, how many LFOs, envelopes, etc. already exist in order to assign an appropriate name
  to the next one to be added, for example in Chainer (the name should include an index/counter for 
  unique identification of the source). */
  int numRegisteredSourcesOfType(ModulationSource* source);

  /** Returns a pointer to the ModulationSource with given name, if a source with that name exists 
  in our array of registered sources. Otherwise, it will return a nullptr. */
  ModulationSource* getSourceByName(const juce::String& sourceName);

  /** Returns a pointer to the ModulationTarget with given name, if a source with that name exists 
  in our array of registered targets. Otherwise, it will return a nullptr. */
  ModulationTarget* getTargetByName(const juce::String& targetName);


  /** \name Misc */

  /** Recalls a state (i.e. all the connections and their settings) from an XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState);

  /** Returns the state (i.e. all the connections and their settings) in form of an XmlElement. */
  virtual XmlElement* getStateAsXml();

protected:

  std::vector<ModulationSource*> availableSources;
  std::vector<ModulationTarget*> availableTargets; // do we need this? maybe for state recall? if not, remove
  std::vector<ModulationTarget*> affectedTargets;
  std::vector<ModulationConnection*> modulationConnections;

  CriticalSection *modLock = nullptr; 

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationManager)
};

//=================================================================================================

class AudioModule; // ModulatableParameter needs a pointer to its owning AudioModule

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


  /** Sets up the pointer to our owner, i.e. the AudioModule that contains this parameter. */
  void setOwnerAudioModule(AudioModule *newOwner) { ownerModule = newOwner; }

  /** We suppose that modulatedValue of the ModulationTarget baseclass has already been computed, 
  i.e. ModulationTarget::computeModulatedValue has been called, such that modulatedValue has a 
  legitimate value with all modulations applied. Here we pull out this modulated value and call our 
  valueChangeCallback with it. The function is supposed to be called per sample for each modulated 
  Parameter. */
  inline void callCallbackWithModulatedValue()
  {
    if( valueChangeCallbackDouble != nullptr )
      valueChangeCallbackDouble->call(modulatedValue);
  }

  /** Overriden to call our callback function with the modulated value. */
  virtual void doModulationUpdate() override
  {
    callCallbackWithModulatedValue();
  }

  /** Returns the unique name of this modulation target that is used to identify it in state 
  recall. */
  virtual juce::String getModulationTargetName() override;

protected:


  AudioModule* ownerModule = nullptr; // needed to create the unique name for state recall
  // maybe we could use the ParameterManager baseclass of AudioModule here

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableParameter)
};

#endif