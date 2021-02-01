#ifndef jura_ModulatableParameter_h
#define jura_ModulatableParameter_h

// todo: 
// -rename files to ModulationSystem.h/cpp - they contain several classes, not only 
//  ModulatableParameter
// -use rs-prefixes for class names

/*

The Modulation System


How to use it:

Client code that wants to use the modulation system for plugin parameters must do the following 
things:

-derive its AudioModule from ModulatableAudioModule (or from a subclass thereof)
-somewhere have a ModulationManager object lying around, a pointer to that object should be passed
 to the constructor call of the ModulatableAudioModule baseclass
-use objects of type ModulatableParameter for its parameters
-use objects of type rsModulatableSlider for the sliders of the to-be-modulated parameters
-every ModulationSource that should be available must be registered with the ModulationManager
 object
-before producing a sample, your AudioModule must call applyModulations on the ModulationManager
 object


How it works:

There are 4 main classes that are relevant for the system: ModulationSource, ModulationTarget, 
ModulationConnection and ModulationManager. ModulationSource and ModulationTarget have a common 
baseclass ModulationParticipant for factoring out some common stuff but the client is usually not 
concerned with this. Sources are things like envelope generators, LFOs, etc. For example, to make 
an envelope generator available for the modulation system, it would have to be a subclass of 
ModulationSource. The subclass is then required to implement the purely virtual 
getModulatorOutputSample() method. There, it is supposed to compute a new output sample of the 
modulator (that value will later be gathered by the ModulationManager and used for modulating 
things). ModulationTargets are typically parameters of some audio-processing algorithm such as the
cutoff frequency of a filter. Typically, clients will only use ModulatableParameter which is a 
subclass of ModulationTarget and Parameter. A ModulationTarget subclass is required to implement 
doModulationUpdate() which is also called per sample. There, they must use the "modulatedValue" 
member (which was set up by the ModulationManager and now - inside the doModulationUpdate call - 
it contains the parameter value with all modulations applied. A ModulatableParameter will call its 
valueChangeCallback (inherited from Parameter), which, for example, may cause a call to 
setCutoff(modulatedValue) on some core filter object. ModulationConnection represents a connection 
between a given ModualtionSource and ModulationTarget. It also contains the strength of the 
connection, i.e. the depth by which some source modulates some target. The ModulationManager 
maintains arrays of all available sources and targets and, importantly, all the connections. It can
also create new and remove existing connections. It is also responsible for calling 
updateModulationValue on all registered ModulationSources and doModulationUpdate on all affected 
ModulationTargets - between these sequences of calls, it gathers up all the source-outputs, and 
applies them to their targets - as defined via the connections.



old:

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
-It doesn't work in Elan's SpiralGenerator - why? maybe it has to do with creation order? There,
 the modulators whould modulate parameters of their parent module, in ToolChain it's the sibling
 modules - and no sources or targets are created in the constructor
 Solution:
 -call setModulationManager in constructor b4 creating sources and targets
 -..
-figure out what happens if the user changes the range of depthParam - how will this affect the 
 meta-value, how can we make sure that the depth parameter is always consistent with its attached 
 metaparameter? how is patch recall affected?
-optimize, where possible - maybe the scaler can be multiplied into the depth (maybe use a 
 scaledDepth member)
 -maybe the ModulationManager should maintain an array usedSources and iterate only over that in
  applyModulations (similar to the treatment of affectedTargets - maybe call the arrays 
  connectedTargets and connectedSources to make the parallelism obvious)
-let the depth be controlled by key and/or velocity - or better, make key and vel available as
 ModulationSource and let Depth-parameters be ModulationTargets (i.e. ModulatableParameters)
 ...or maybe not the midi key/vel values but something more physical, for example Note-Frequency
 and Amplitude, i.e. freq = pitchToFreq(key+pitchbend), amp = vel/127.

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
  virtual void setModulationManager(ModulationManager* managerToUse) { modManager = managerToUse; }
  // maybe get rid of this and demand it to be passed to the constructor

  /** Returns a pointer to the ModulationManager object that is used for registering 
  ModulationSources and ModulationTargets. */
  ModulationManager* getModulationManager() { return modManager; }

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

  /** Returns a reference to the array of ModulationConnection objects. */
  const std::vector<ModulationConnection*>& getModulationConnections();


protected:

  ModulationManager* modManager = nullptr;

  static std::vector<ModulationSource*> dummySources;
  static std::vector<ModulationTarget*> dummyTargets;
  static std::vector<ModulationConnection*> dummyConnections;
  // These are empty dummy arrays to which references will be returned by 
  // getAvailableModulationSources/Targets/getModulationConnections in case our modManager is a 
  // nullptr. This is somehow ugly design, maybe use the "Null Object" pattern instead:
  // https://sourcemaking.com/design_patterns/null_object
  // somewhere, we should have a default "null" ModulationManager object lying around to which
  // our pointer is initialized

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationParticipant)
};

//=================================================================================================

/** Baseclass for modulation sources. */

class JUCE_API ModulationSource : public ModulationParticipant
{

public:

  /** Constructor */
  ModulationSource(ModulationManager* managerToUse = nullptr) 
    : ModulationParticipant(managerToUse) { }

  /** Destructor */
  virtual ~ModulationSource();


  /** Override this function in your subclass to produce one modulator output sample at a time. */
  virtual double renderModulation() = 0;
  // was formerly getModulatorOutputSample - but the old name was ambiguous with 
  // getModulationValue, the getter for the modulation signal after it has been rendered, 
  // especially in the polyphonic subclass

  /** This is called per sample from the ModulationManager and it updates our modValue member 
  variable by calling the virtual renderModulation() function (that subclasses must override) and 
  using its output for the modValue. The rsVoiceManager parameter is irrelevant in this class here 
  but it becomes relevant in the ModulationSourcePoly subclass, where we override this method in 
  order to update the modulator outputs for all active voices, so the voice manager must be already
  included in the function signature here. */
  virtual void updateModulationValue(rsVoiceManager* voiceManager) 
  { modValue = renderModulation(); }
  // maybe rename - the singular "Value" part does not apply anymore for the poly subclass which
  // updates all values for all voices - maybe updateModulationOutput
  // why is this virtual? this is not supposed to be overriden

  /** Sets up a name for this ModulationSource. This should be unique among all the available 
  ModulationSource names, so it can be used to identify the source in state recall. */
  void setModulationSourceName(const juce::String& newName) { modSourceName = newName; }
  // todo: maybe in debug mode, assert that the name is unique (requires iterating through all 
  // already registered sources and comaring the name)

  /** Sets a name that should be used in dropdown list when connecting a mod-source. If you don't
  set this up, the name that was set by setModulationSourceName will be used by default. */
  void setModulationSourceDisplayName(const juce::String& newName) { displayName = newName; }

  /** Returns the name of this ModulationSource. This is the (supposed to be unique) name that will
  used for identifying the source state recall. */
  juce::String getModulationSourceName() const { return modSourceName; }

  /** Returns the name that should be used to identify the source in the dropdown menu on the gui.
  This may potentially be different from the name used for state recall. */
  juce::String getModulationSourceDisplayName() const;

  /** Returns true, if this source is connected to at least one ModulationTarget. Currently, this 
  is quite expensive to call (requires to loop through the connections). I think, that's acceptable
  at the momen because it's not supposed to be called often (like, per sample). If that situation 
  changes, it may be optimized by using a numConnectedTargets member just like ModulationTarget has
  a numConenctedSources member. */
  bool hasConnectedTargets() const;

  /** Returns the current output value of the modulation source. This is supposed to be called 
  after updateModulationValue has been called. */
  inline double getModulationValue() const { return modValue; }

  /** This is a provision for implementing polyphony in subclass ModulationSourcePoly. There, it 
  returns the modulator output for a given voice. Here, we just return the same value for all 
  voices, because this baseclass is only monophonic and when a monophonic modulator is asked for a
  voice output, it makes most sense to just use the same value for all voices. */
  virtual double getVoiceModulationValue(int voiceIndex) const { return modValue; }

protected:

  double modValue = 0;

  juce::String modSourceName = "ModulationSource";
  juce::String displayName   = "";

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationSource)
};

//=================================================================================================

/** Baseclass for modulation targets.  */

class JUCE_API ModulationTarget : public ModulationParticipant
{

public:

  /** Constructor */
  ModulationTarget(ModulationManager* managerToUse = nullptr) 
    : ModulationParticipant(managerToUse) { }

  /** Destructor */
  virtual ~ModulationTarget();

  /** Must be overriden by subclasses to do whatever they need to do after a modulatedValue has 
  been computed (for example, ModulatableParameter invokes the setter-callback which in turn 
  updates the corresponding value in the core dsp algorithm). */
  virtual void doModulationUpdate(double modulatedValue)
  {
    // We need an empty baseclass implementation because in the destructor of a plugin, 
    // doModulationUpdate would otherwise (in case of a purely virtual function) get called with a 
    // null-reference (or something), in ModulationManager::removeConnection when the modulateble 
    // parameter deletes itself. ToDo: figure out and document why this is the case
  }

  /** Like doModulationUpdate, but for polyphonic modulations. Overriden in 
  ModuldatableParameterPoly. Monophonic modulation targets can ignore this. */
  virtual void doVoiceModulationUpdate(double modulatedValue, int voiceIndex) { }


  /** \name Setup */

  /** Sets the nominal, unmodulated value. This will be used as reference, when a modulated value 
  will be computed. */
  void setUnmodulatedValue(double newValue)
  {
    unmodulatedValue = newValue;

    modulatedValue   = unmodulatedValue; 
    // Elan's smoother needs this - it makes sure that the modulated value is always the same as
    // the unmodulated value, even in cases in case where no modulations are applied (in which case 
    // initModulatedValue is not called (per sample)). It works only, if the smoother calls 
    // setUnmodulatedValue (per sample) before applyModulations is called for the same sample, 
    // otherwise the modulations don't work anymore, because the modulatedValue is reset after 
    // modulations have been applied. It's a bit hacky but probably needs to be taken into account
    // when designing a parameter smoother for jura.
    // do we also need to set more values in the array in the polyphonic case?
  }

  /** Adds a ModulationSource to this ModulationTarget. The amount of modulation is initially 0. */
  void addModulationSource(ModulationSource* source);

  /** Removes a ModulationSource from this ModulationTarget. */
  void removeModulationSource(ModulationSource* source);

  /** Sets the minimum value of the allowed modulation range. These min/max values for the range 
  are used to clip the final modulated value (after all modulations have been applied) to a sane 
  range. For example, a cutoff frequency of a filter should perhaps be clipped at 0 and 
  sampleRate/2, otherwise the filter may go crazy. */
  inline void setModulationRangeMin(double newMin) { rangeMin = newMin; }

  /** Sets the maximum value of the allowed modulation range.  */
  inline void setModulationRangeMax(double newMax) { rangeMax = newMax; }

  /** Convenience function to set up all the range/depth/mode parameters at once. */
  inline void setDefaultModParameters(double newRangeMin, double newRangeMax, 
    double newDefaultDepthMin, double newDefaultDepthMax, int defaultMode, 
    double initialModDepth = 0.0)
  {
    rangeMin = newRangeMin;
    rangeMax = newRangeMax;
    defaultDepthMin = newDefaultDepthMin;
    defaultDepthMax = newDefaultDepthMax;
    defaultModMode  = defaultMode;
    initialDepth    = initialModDepth;
  }


  /** \name Inquiry */

  /** Returns the minimum value of the allowed modulation range. */
  inline double getModulationRangeMin() const { return rangeMin; }

  /** Returns the maximum value of the allowed modulation range. */
  inline double getModulationRangeMax() const { return rangeMax; }

  /** Returns the default value for the modulation depth minimum for a new connection going into 
  this target. The user can later change that. */
  inline double getDefaultModulationDepthMin() const { return defaultDepthMin; }

  /** Returns the default value for the modulation depth maximum for a new connection. */
  inline double getDefaultModulationDepthMax() const { return defaultDepthMax; }

  /** Returns the initial modulation depth that will be used for new incoming connections. */
  inline double getInitialModulationDepth() const { return initialDepth;    }

  /** Returns the default modulation mode which is the mode that a newly connected connection to 
  this target will have. This can be any of the values defined in 
  ModulationConnection::modModes. */
  inline int getDefaultModulationMode() const { return defaultModMode; } 

  /** Returns true, if there's a connection between this ModulationTarget and the given 
  ModulationSource. */
  bool isConnectedTo(ModulationSource* source) const;

  /** Returns a pointer to the ModulationConnection object that connects this target to the given 
  source (if any, nullptr otherwise); */
  ModulationConnection* getConnectionTo(ModulationSource* source);

  /** Returns a vector of pointers to the ModulationSources which are connected to this 
  ModulationTarget. They will appear in the order of the connections in the modManager. */
  std::vector<ModulationSource*> getConnectedSources() const;

  /** Returns a vector of pointers to the available ModulationSources which are not connected to 
  this ModulationTarget. They will appear in the order in which they were registered with the
  modManager. */
  std::vector<ModulationSource*> getDisconnectedSources() const;

  /** Returns a vector of pointers to ModulationConnections that are incoming into this 
  ModulationTarget. */
  std::vector<ModulationConnection*> getConnections() const;

  /** This function must be overriden by subclasses to return a unique name that can be used to 
  identify the target in state recall. */
  virtual juce::String getModulationTargetName() = 0;

  /** Returns true, if this target is connected to at least one ModulationSource. */
  bool hasConnectedSources() const;


  /** \name Misc */

  /** Initializes the modulated value by setting it to the unmodulated value. */
  inline void initModulatedValue()
  {
    modulatedValue = unmodulatedValue;
  }

  /** Function to retrieve the modulated value after all modulations have been applied. This may 
  also include a clipping function, such that the returned value is restricted to some allowable
  range. */
  inline double getModulatedValue() const
  {
    return RAPT::rsClip(modulatedValue, rangeMin, rangeMax);
  }

  /** Returns the unmodulated value, which is the base value when no modulation is applied. */
  inline double getUnmodulatedValue() const
  {
    return unmodulatedValue;
  }

  /** Returns a pointer to the modulated value to facilitate accumulation of the contributions 
  from all the connections */
  inline double* getPointerToModulatedValue() { return &modulatedValue; }


protected:

  double unmodulatedValue = 0;
  double modulatedValue   = 0; 
  double rangeMin = -INF, rangeMax = INF; 
  double defaultDepthMin = -1, defaultDepthMax = 1;
  double initialDepth = 0.0;
  int    defaultModMode = 0; // absolute
  int    numConnectedSources = 0;

  friend class ModulationManager;
  friend class ModulationManagerPoly; 
  // maybe get rid by providing accessors

  // ToDo:
  // -Maybe we don't need to store the modulatedValue here. This oculd perhaps be handled by the 
  //  ModulationManager as a local variable, quite similar to what we do with the modValues array
  //  ModulationManagerPoly.
  // -Maybe it makes sense to keep track of the connected source via a numConnectedSources 
  //  variable....

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationTarget);
};

//=================================================================================================

/** This class represents a connection from a ModulationSource to a ModulationTarget with 
adjustable modulation amount. */

class JUCE_API ModulationConnection
{

public:

  /** Selects one of the modes in which the modulation can be applied. In the formulas below, u is
  the unmodulated base value, m the modulator output and d the modulation depth. */
  enum modModes
  {
    ABSOLUTE = 0,    // d * m
    RELATIVE,        // d * m * u
    EXPONENTIAL,     // u * 2^(d*m) - u
    MULTIPLICATIVE   // u * m^d     - u
  };
  // ToDo: use an enum class

  /** Constructor. You should pass a ModulationSource and a ModulationTarget. If you want to enable
  meta-control for the modulation-depth associated with the connection, you may pass the pointer to
  the MetaParameterManager object that should be used for this. */
  ModulationConnection(ModulationSource* source, ModulationTarget* target, 
    MetaParameterManager* metaManager);
    // maybe disallow metaManager to be a nullptr - or handle the case with a null object in 
    // MetaControlledParameter

  /** Destructor. */
  virtual ~ModulationConnection(); // maybe make non-virtual

  /** Sets the modulation depth for this connection. This will actually update the Parameter object 
  that is associated with the depth. */
  void setDepth(double newDepth)
  {
    depthParam->setValue(newDepth, true, true); // the 2nd true triggers the callback call which 
                                                // typically updates the dsp object (i think)
    // modManager->modulationDepthChanged(this);
  }

  /** Sets the parameter range and value of the modulation depths. */
  void setDepthRangeAndValue(double newMin, double newMax, double newValue)
  { depthParam->setRangeAndValue(newMin, newMax, newValue, true, true); }

  /** Sets the mode in which this connection should apply the modulation source output value to the
  modualtion target. @see modModes. */
  inline void setMode(int newMode) { mode = newMode; }

  /** Returns the modulation mode of this connection. */
  inline int getMode() { return mode; }
    
  /** Returns the Parameter object that controls the amount of modulation */
  MetaControlledParameter* getDepthParameter() const { return depthParam; }


  /** Given a modulator output value m, a modulatio depth d and an unmodulated parameter value u,
  this function applies the modulation to the target value t. The mode controls, which formula
  is used for the application. Factored out common code for applying modulations mono- and 
  polyphonically. */
  static inline void apply(double m, double d, double u, double* t, int mode)
  {
    double z;
    switch(mode)  // maybe use function pointer instead of switch
    {
    case ABSOLUTE:       z = d * m;             break;
    case RELATIVE:       z = d * m * u;         break;
    case EXPONENTIAL:    z = u * exp2(d*m) - u; break; // maybe rename function to pow2
    case MULTIPLICATIVE: 
    {
      z = RAPT::rsClip(u * pow(m, d) - u, -1.e100, +1.e100);
      // Multiplicative mode may produce infinity when it gets an input like m=0, d=-1: 0^-1 = 1/0.
      // That in itself might not be a problem due to later being clipped to the max target value 
      // but if one modulator produces inf and another one -inf because u is negative in the 2nd, 
      // we would get inf-inf = NaN so we need to clip here already.
    } break;
    default:             z = 0;
    }
    *t += z;
  }
  // ToDo:
  // Maybe have a true multiplicative mode that actually multiplies t by m^d because sometimes 
  // that's really a desirable behavior. But then the order of the connections becomes relevant, so
  // we need a means for the user to reorder them. Then, the current faux-multiplicative mode 
  // should be renamed into something else and the real multiplicative mode should get the name 
  // "Multiplicative" - that needs to invoke the patch conversion when old patches are loaded...
  // the real and faux mode should give the same results only if the multiplicative connection
  // is the first of all (and in particular if its the only connection). Maybe we could have 
  // another faux mode that uses the depth not as exponent but as scaler...

  /** Applies the source-value to the target-value that is stored in the ModulationTarget object of
  this connection, taking into account the modulation depth. This "modulatedValue" member of 
  ModulationTarget is used for monophonic modulations. */
  inline void apply() const
  {
    apply(source->getModulationValue(), depth, target->getUnmodulatedValue(),
      target->getPointerToModulatedValue(), mode);
  }

  /** Applies the modulation value of the sources for the given voice index to the target value 
  where the location of the target value can be controlled by the caller. This is used for 
  polyphonic modulations, where the target value is not stored in the ModulationTarget object but
  buffered elesewhere in order to save memory. */
  inline void applyVoice(double* targetValue, int voiceIndex) const
  {
    apply(source->getVoiceModulationValue(voiceIndex), depth, target->getUnmodulatedValue(), 
      targetValue, mode);
  }
  // maybe find a better name



  /** Returns a xml element containing the information about this connection (needed for state 
  save/recall in ModulationManager - todo: maybe move it to there) */
  XmlElement* getAsXml(); // should be const

  /** Returns a (const) pointer to the source, so you can inquire something about it. */
  const ModulationSource* getSource() const { return source; }

  /** Returns a (const) pointer to the target, so you can inquire something about it. */
  const ModulationTarget* getTarget() const { return target; }

  /** Returns a (non-const) pointer to the target, so you can do stuff with it, like invoking the
  callbacks */
  ModulationTarget* getTarget() { return target; }

protected:

  /** Sets the modulation depth for this connection. Used as target callback in the depthParam. To 
  keep our "depth" member consistent with the value of the "depthParam" Parameter object. This is 
  used only internally, client code should use setDepth. */
  void setDepthMember(double newDepth) { depth = newDepth; }

  ModulationSource* source;
  ModulationTarget* target;
  double depth;              // modulation depth
  int mode = ABSOLUTE;       // application mode for modulation signal

  MetaControlledParameter* depthParam; // maybe it should be a ModulatableParameter? but that may
                                       // raise some issues - maybe later...

  friend class ModulationSource;
  friend class ModulationTarget;
  friend class ModulationManager;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationConnection);
};

//=================================================================================================

/** This class is responsible for keeping a list of all available ModulationSources and 
ModulationTargets and making it possible for ModulationTargets to connect themselves to any number 
of these sources. It's also responsible for computing all modulated values and updating the targets
accordingly. */

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

  /** Same as applyModulations() but it doesn't acquire the mutex-lock. This should be used in cases 
  where the caller already has acquired the lock (i.e. it already holds the lock for the same 
  CriticalSection object that was initially passed to our constructor). */
  virtual void applyModulationsNoLock();


  /** \name Connection setup */

  /** Adds a connection between the given source and target. */
  void addConnection(ModulationSource* source, ModulationTarget* target);

  /** Adds the passed ModulationConnection to our array. */
  void addConnection(ModulationConnection* connection);

  /** Removes a connection between the given source and target. */
  void removeConnection(ModulationSource* source, ModulationTarget* target);

  /** Removes all modulation connections that involve the given source. */
  void removeConnectionsWith(ModulationSource* source);

  /** Removes all modulation connections that involve the given target. */
  void removeConnectionsWith(ModulationTarget* target);

  /** Removes all ModulationConnections. */
  void removeAllConnections();

  /** Removes the connection with given index. */
  void removeConnection(int index, bool updateAffectTargets);
  // maybe set updateAffectTargets to true by default

  /** Resets all the range limits for all registered modulation targets to +-inf. */
  //void resetAllTargetRangeLimits();


  /** \name Registration of sources and targets */

  /** Registers the given ModulationSource to make it available to ModulationTargets. */
  void registerModulationSource(ModulationSource* source);

  /** De-registers a ModulationSource. */
  void deRegisterModulationSource(ModulationSource* source);

  /** De-registers all ModulationSources. */
  void deRegisterAllSources();

  /** Registers the given ModulationTarget. */
  void registerModulationTarget(ModulationTarget* target);

  /** De-registers a ModulationTarget. */
  void deRegisterModulationTarget(ModulationTarget* target);

  /** De-registers all ModulationTargets. */
  void deRegisterAllTargets();


  /** \name Inquiry */

  /** Returns true if there's a connection between the given source and target. */
  bool isConnected(const ModulationSource* source, const ModulationTarget* target) const;

  /** Returns a pointer to a ModulationConnection object between the given source and target, if 
  such a conncetion exists. Otherwise, it will return a nullptr. */
  ModulationConnection* getConnectionBetween(const ModulationSource* source, 
    const ModulationTarget* target) const;

  /** Returns the number of ModulationConnections. */
  inline int getNumConnections() const
  { 
    ScopedLock scopedLock(*modLock); 
    return size(modulationConnections); 
  }

  /** Returns a reference to our vector of available ModulationSources. */
  const std::vector<ModulationSource*>& getAvailableModulationSources() 
  { 
    ScopedLock scopedLock(*modLock); 
    return availableSources; 
  }

  /** Returns a reference to our vector of available ModulationTargets. */
  const std::vector<ModulationTarget*>& getAvailableModulationTargets() 
  { 
    ScopedLock scopedLock(*modLock); 
    return availableTargets; 
  }

  /** Returns a reference to our vector of ModulationConnections. */
  const std::vector<ModulationConnection*>& getModulationConnections()
  { 
    ScopedLock scopedLock(*modLock); 
    return modulationConnections; 
  }

  /** Given a pointer to a ModulationSource of some type, this function returns the number of 
  ModulationSources of the same type that are registered here. This is used to figure out, for 
  example, how many LFOs, envelopes, etc. already exist in order to assign an appropriate name
  to the next one to be added (it's actually not currently used). */
  //int numRegisteredSourcesOfType(ModulationSource* source);

  /** Returns a pointer to the ModulationSource with given name, if a source with that name exists 
  in our array of registered sources. Otherwise, it will return a nullptr. */
  ModulationSource* getSourceByName(const juce::String& sourceName);

  /** Returns a pointer to the ModulationTarget with given name, if a source with that name exists 
  in our array of registered targets. Otherwise, it will return a nullptr. */
  ModulationTarget* getTargetByName(const juce::String& targetName);

  /** Returns true, if any of our affected targets has some range limits set up. If so, we need to
  store them in a dedicated xml-element when producing a state xml. */
  bool needsToStoreRangeLimits();


  /** \name Misc */

  /** Sets the MetaParameterManager that should be used to attach meta-parameters to the 
  depth-parameters of the modulation connections. */
  virtual void setMetaParameterManager(MetaParameterManager* managerToUse);

  /** Sets the voice manager that should be used in the polyphonic case. Can be ignored in for
  monophonic plugins. */
  virtual void setVoiceManager(rsVoiceManager* managerToUse);

  /** Recalls a state (i.e. all the connections and their settings) from an XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState);

  /** Returns the state (i.e. all the connections and their settings) in form of an XmlElement. */
  virtual XmlElement* getStateAsXml();

  /** Returns the voice manager that is used for polyphonic modulation, if applicable. */
  rsVoiceManager* getVoiceManager() const { return voiceManager; }


  /** Updates our affectedTargets array. Called from the connection-removal functions in order to 
  remove the target from the affectedTargets, in case it has no incoming connections anymore after 
  connection removal. */
  void updateAffectedTargetsArray();

protected:

  /** Just appends the given connction to the end of our modulationConnections. Intended to be 
  overriden by subclass which need to enforce a certain ordering of the connections. */
  virtual void addConnectionToArray(ModulationConnection* connection)
  { modulationConnections.push_back(connection); }

  /** Just removes the pointer to the connection at the given index from our modulationConnections
  array. Intended to be overriden by subclasses which need to perform additional operations when 
  this happens. */
  virtual void removeConnectionFromArray(int index)
  { remove(modulationConnections, index); }

  /** Tries to cast the passed ModulationTarget into an ObservableModulationTarget and if this is 
  successful, it sends out the modulation change notifiaction for it. */
  void sendModulationChangeNotificationFor(ModulationTarget* target);

  std::vector<ModulationSource*> availableSources;
  std::vector<ModulationTarget*> availableTargets;
  std::vector<ModulationTarget*> affectedTargets;
  std::vector<ModulationConnection*> modulationConnections;

  CriticalSection *modLock = nullptr; 
  MetaParameterManager*  metaManager  = nullptr;  // maybe use a null object instead - or maybe not

  rsVoiceManager*        voiceManager = nullptr; // move to poly subclass

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationManager)
};

//=================================================================================================

/** A baseclass for objects that need to observe ModulationTarget objects. */

class JUCE_API ModulationTargetObserver
{

public:

  ModulationTargetObserver() = default;
  virtual ~ModulationTargetObserver() = default;

  /** Subclasses need to override this in order to respond to changes of the modulation settings of 
  their observed modulation target. */
  virtual void modulationsChanged() = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationTargetObserver)
};

/** A subclass of ModulationTarget that allows to be monitored by observer objects. For example,
a slider on a gui could keep track of whether its underlying parameter has modulations applied to 
it and if so, change its appearance. */

class JUCE_API ObservableModulationTarget : public ModulationTarget
{

public:

  /** Constructor */
  ObservableModulationTarget(ModulationManager* managerToUse = nullptr) 
    : ModulationTarget(managerToUse) {}


  void registerModulationTargetObserver(ModulationTargetObserver* observer)
  {
    appendIfNotAlreadyThere(modTargetObservers, observer);
  }

  void deRegisterModulationTargetObserver(ModulationTargetObserver* observer)
  {
    removeFirstOccurrence(modTargetObservers, observer);
  }

  void sendModulationsChangedNotification()
  {
    for(size_t i = 0; i < modTargetObservers.size(); i++)
      modTargetObservers[i]->modulationsChanged();
  }

protected:

  std::vector<ModulationTargetObserver*> modTargetObservers;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ObservableModulationTarget)
};

//=================================================================================================

class AudioModule; // ModulatableParameter needs a pointer to its owning AudioModule

/** A subclass of Parameter that is suitable as target for modulations by also being a subclass of
ModulationTarget. */

class JUCE_API ModulatableParameter : public MetaControlledParameter, 
  public ObservableModulationTarget
{

public:

  /** Constructor */
  ModulatableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0)
    : MetaControlledParameter(name, min, max, defaultValue, scaling, interval) 
  {
    modulatedValue = unmodulatedValue = value;
  }


  /** \name Setup */


  void setNormalizedValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // Note: setValue is not overriden because the baseclass implementation actually calls the 
  // overriden setNormalizedValue function from here, so it's enough to override this

  void setSmoothedValue(double newValue) override;



  /** Sets up the pointer to our owner, i.e. the AudioModule that contains this parameter (needed 
  for unique identification of this parameter in the tree of AudioModules when a state is 
  recalled). */
  void setOwnerAudioModule(AudioModule *newOwner) { ownerModule = newOwner; }


  /** \name Misc */

  /** comment out of date - todo: update
  We suppose that modulatedValue of the ModulationTarget baseclass has already been computed, 
  such that modulatedValue has a legitimate value with all modulations applied. Here we pull out 
  this modulated value and call our valueChangeCallback with it. The function is supposed to be 
  called per sample for each modulated Parameter. ModulationManager will take care of this. */
  inline void callCallback(double modulatedValue)
  {
    //jassert(RAPT::rsIsFiniteNumber(getModulatedValue()));
    if( valueChangeCallbackDouble != nullptr )
      valueChangeCallbackDouble->call(modulatedValue);
    if( valueChangeCallbackInt != nullptr )
      valueChangeCallbackInt->call((int)modulatedValue);
  }

  /** Overriden to call our callback function with the modulated value. */
  void doModulationUpdate(double modulatedValue) override
  {
    callCallback(modulatedValue);
  }

  /** Returns the unique name of this modulation target that is used to identify it in state 
  recall. */
  virtual juce::String getModulationTargetName() override;

protected:

  AudioModule* ownerModule = nullptr; // needed to create the unique name for state recall
  // maybe we could use the ParameterManager baseclass of AudioModule here

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableParameter)
};

/** This class is like ModulatebleParameter except that it uses the other std::function based
callback, so if you use this callback mechanism, use this class for your parameters */ 
class JUCE_API ModulatableParameter2 : public ModulatableParameter
{
  using ModulatableParameter::ModulatableParameter; // inherit baseclass constructors
  virtual void doModulationUpdate(double modulatedValue) override
  {
    valueChangeCallbackFunction(modulatedValue);
  }
};
// Elan uses this - eventually, i should probably switch to this too and get rid of all the other
// callback types in class Parameter, but the current mechanism may actually be more performant 
// than std::function. Benchmark this first by setting up a realistically complex configuration of
// modulators and parameters with ModulatableParameter and ModulatableParameter2 and compare their
// performances.


//#################################################################################################

// The polyphonic version of the modulation system


class JUCE_API ModulationManagerPoly : public ModulationManager
{

public:

  ModulationManagerPoly(CriticalSection* lockToUse) : ModulationManager(lockToUse) {}

  /** Overriden to apply the polyphonic modulations as well. */
  void applyModulationsNoLock() override;

protected:

  /** Applies the modulations for the given voice. Called from applyModulationsNoLock in a loop
  over the active voices. */
  void applyVoiceModulations(int voiceIndex);

  /** Overriden to make sure that the connections are in an order that ensures that connections 
  with the same target are adjacent in the modulationConnections array. This facilitates looping 
  over connections in applyVoiceModulations for one target at a time. This is needed because here 
  because we render the modulated values into our modulatedValues array and not in the member 
  stored in ModulationTarget. 
  ToDo: maybe do it like this in the monophonic case, too. This would allow to throw away the 
  modulatedValue member of ModulationTarget. The connections would always be ordered by target, 
  even in the baseclass, so we would also get rid of the need to override. */
  void addConnectionToArray(ModulationConnection* connection) override;

  /** Overriden to potentially resize the modulatedValues array and call the per voice callbacks
  for the active voices, in cases when the last connection from a target is removed. */
  void removeConnectionFromArray(int index) override;

  /** Storage buffer for the modulated values. Each ModulationTarget that has incoming modulation 
  connections gets a correpsonding slot in this array here. It serves as a substitute for the 
  modulatedValue member of the ModulationTarget class (that is used for monophonic modualtions) 
  when computing and applying the modulations, because we really don't want to have to store an 
  array in each target because that would starkly increase the memory footprint of a 
  ModulationTarget object. 
  ToDo: maybe use the same meachanism for monophonic modulations, too. That would save one byte
  of storage space (a member of type double) for ModulationTarget objects and the code would be 
  more consistent.  */
  std::vector<double> modulatedValues;


  int numDistinctActiveTargets = 0;
  // Isn't that redundant with numAffectedTargets.size() and also modulatedValues.size()? It 
  // represents the number of targets that are connected...active seems a misnaming since it 
  // suggest having something to do with active voice which it doesn't. It's the number of 
  // parameters that are being modulated. figure out and if so, delete it

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationManagerPoly)
};

//=================================================================================================

/** Baseclass for polyphonic modulation sources.

...

The rules for applying polyphonic sources to monophonic targets are as follows:

  -If at least one voice is active (i.e. at least one note is playing), the monophonic target 
   receives the modulation signal from the newest active voice, i.e. the latest triggered note that
   is still being held.
  -If no voice is playing, the monophonic target receives whatever your overriden renderModulation
   function produces. But you don't have to override it in which case it receives zero, because 
   that's what the inherited baseclass function will then produce. 

So, your subclass is only *required* to override renderVoiceModulation and it *may optionally* 
also override renderModulation in cases where it's not appropriate to produce a zero modulation 
when no voice is playing (which, i guess, is uncommon - so mostly it should suffice to override
renderVoiceModulation -> less boilerplate).

...i'm not sure, if that's the best way to do it - experimentation is needed
maybe it could also make sense to continue producing the last signal value from the latest voice
after it was turned of - i.e. always use value of the newest voice, regardless whether this voice 
is still active or not. but his implies that the sound of a patch will depend on the latest note
that was ever produced - even after the note has long been died out - that may be bad for sound 
design: the synth is in a certain state (due to the last played note), the user creates a sound,
saves it and when they reload it at another day, it will sound different, because the synth is in
a different state. that would be bad! ..it seems like having a state which depends on the last 
played note and which affects the sound output is not a good idea.

*/

class JUCE_API ModulationSourcePoly : public ModulationSource
{

public:


  ModulationSourcePoly(ModulationManager* managerToUse = nullptr)
    : ModulationSource(managerToUse) { }


  virtual void setMonophonic(bool shouldBeMonophonic) { monophonic = shouldBeMonophonic; }


  /** Must be overriden by subclasses to produce a modulator output sample for the given voice 
  index. */
  virtual double renderVoiceModulation(int voiceIndex) = 0;

  /** We override this purely virtual method inherited from ModulationSource here to produce a zero
  value. This rendering function is invoked when there is no active voice, i.e. no note is playing.
  If your modulator subclass wants to produce something other than zero in such a case, you may
  override it in your subclass of ModulationSourcePoly again. */
  double renderModulation() override { return 0.0; }



  double getVoiceModulationValue(int voiceIndex) const override
  { 
    jassert(voiceIndex >= 0 && voiceIndex < (int)modValues.size());
    return modValues[voiceIndex]; 
  }


  void updateModulationValue(rsVoiceManager* voiceManager) override
  { 
    if(voiceManager == nullptr) {
      modValue = renderModulation();
      return; }
      // This is supposed to happen when either the monophonic ModulationManager baseclass is used
      // or when ModulationManagerPoly is used but there are currently no active voices (it then
      // falls back to calling the basclass method)

    jassert(modValues.size() >= voiceManager->getMaxNumVoices());
    if(!monophonic) {
      for(int i = 0; i < voiceManager->getNumActiveVoices(); i++)  {
        int k = voiceManager->getActiveVoiceIndex(i);
        modValues[k] = renderVoiceModulation(k); }}
    else {
      int newestVoice = voiceManager->getNewestActiveVoice();
      double val = renderVoiceModulation(newestVoice);
      for(int i = 0; i < voiceManager->getNumActiveVoices(); i++)  {
        int k = voiceManager->getActiveVoiceIndex(i);
        modValues[k] = val; }}
  }


  virtual void allocateVoiceOutputPins(rsVoiceManager* voiceManager)
  {
    if(voiceManager == nullptr)
      modValues.resize(0);       // or maybe use size 1? ...nah - array should not be used in this case
    else
      modValues.resize(voiceManager->getMaxNumVoices());
  }

protected:

  std::vector<double> modValues; // replaces "modValue" member of monophonic baseclass
  bool monophonic = false; 
  // if true, the modulation source should behave as if it was a monophonic one


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulationSourcePoly)
};

//=================================================================================================

class JUCE_API ModulatableParameterPoly : public ModulatableParameter
{

public:


  /** Constructor */
  ModulatableParameterPoly(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0)
    : ModulatableParameter(name, min, max, defaultValue, scaling, interval) {}


  //void setValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // called from gui?

  void setNormalizedValue(double newValue, bool sendNotification, bool callCallbacks) override;
  // called from automation system?

  void setSmoothedValue(double newValue) override;
  // called from smoothing manager?


  void setMonophonic(bool shouldBeMonophonic) { monophonic = shouldBeMonophonic; }


  void setValueChangeCallbackPoly(std::function<void(double, int)> cb)
  {
    ScopedPointerLock spl(mutex);
    valueChangeCallbackPoly = cb;
    //callValueChangeCallbackPoly(value);
  }
  // maybe make the parameter a const reference


  bool isMonophonic() const { return monophonic; }

  int getNewestActiveVoice()
  {
    ModulationManager* modMan = getModulationManager();    jassert(modMan);
    rsVoiceManager*  voiceMan = modMan->getVoiceManager(); jassert(voiceMan);
    return voiceMan->getNewestActiveVoice();
    // todo: try to optimize this and make it const
  }


  /** Overriden to call our callback function with the modulated value. */
  void doVoiceModulationUpdate(double modulatedValue, int voiceIndex) override
  {
    jassert(voiceIndex >= 0);
    if(!monophonic)  // experimental
      valueChangeCallbackPoly(modulatedValue, voiceIndex);
    else if(voiceIndex == getNewestActiveVoice())  
      valueChangeCallbackPoly(modulatedValue, voiceIndex);
  }


  // maybe we also need to override setNormalizedValue, setSmoothedValue? ..but maybe not - these
  // should modifiers should be baked into the unmodulatedValue before modulation is applied..but:
  // the per voice callbacks will never get called on parameter changes, when the parameter has no
  // connected modulator. maybe we should call them for all active voices there ourselves


  virtual juce::String getModulationTargetName() override
  {
    return ModulatableParameter::getModulationTargetName();
  }
  // idk, why we don't inherit that method -> figure out

  /** Calls the value change callback for the given voice with the unmodulated value. This is 
  needed when a new note was triggered and the parameter has no connected modulation sources. 
  Called from AudioModulePoly::noteOn. */
  void callCallbackForVoice(int voiceIndex)
  { valueChangeCallbackPoly(modulatedValue, voiceIndex); }

  /** Calls the value change callback for all curently active voices with the unmodulated value. 
  This is needed when the parameter was changed on the GUI when it has no connected modulation 
  sources. Called from our override setValue etc. methods. */
  void callCallbacksForActiveVoices();


protected:



  //typedef GenericMemberFunctionCallback1<void, double> SetValueCallback;
  //std::vector<SetValueCallback*> valueChangeCallbacks;
  // for testing, we only use a "double" callback - it's very ugly design in Parameter, to have 
  // pointers to all 3 kinds of callbacks (double, int, bool) - maybe refactor and/or templatize 
  // the design...or maybe just use std::function, as Elan does - but before doing such a switch, 
  // figure out if the old way is more performant - if that's the case, it may be better to keep 
  // it - and it may actually even better to not use std::function here, too. -> Do performance 
  // tests that compare ModulatableParameter with ModulatableParameter2. The higher memory cost
  // for std::function is probably irrelevant because Parameter is already quite heavyweight 
  // anyway...da macht das biﬂchen mehr den kohl auch nicht mehr fett

  //typedef std::function<void(double)> SetValueCallback;
  //std::vector<SetValueCallback*> valueChangeCallbacks;


  //typedef std::function<void(double, int)> SetValueCallbackPoly;
  //std::vector<SetValueCallbackPoly> valueChangeCallbackPoly;


  std::function<void(double, int)> valueChangeCallbackPoly;

  bool monophonic = false; 
  // if true, the modulation target should behave as if it was a monophonic one


private:

  virtual void callValueChangeCallbacks(double argument) override {}
  // this inherited function should not be used anymore - instead, the function with same name but 
  // different signature: callValueChangeCallbacks(int, double*, int*) should be used here


  //GenericMemberFunctionCallback1<void, double> *valueChangeCallbackDouble

};
// -maybe it should have its own callCallbacks function that takes the voice index as second 
//  parameter
// -or we have an array of callback objects
// -or, use the same simple callback object but pointer arithemetic (based on the voice index) to 
//  address the appropriate object - this may be a bit dirty and hacky internally but would be the 
//  leanest solution - maybe we would have to change the callback class like this:
// ReturnType call(const ArgumentType argument, int index)
// {
//    CalleeObjectType* indexedCallee = calleeObject + sizeof(CalleeObjectType);
//    return (*indexedCallee.*memberToCall)(argument);
// }
//  maybe we would have to make a subclass of the callback class
//  BUT: this strategy assumes that the modulated dsp objects are all in a contiguous array - which 
//  may not be the case (think of recursive composition of objects) ...unless we allow for a sort of 
//  object-stride that is not necessarily equal to sizeof(CalleeObjectType)
//  -but maybe that restriction to dsp objects in contiguous arrays is not a big issue? we probably 
//   always want them to be in an array anyway?
// -we probably need an array of callee-objects instead of a single object
//  -maybe have an calleeArray - actually, we just need to augment SpecificMemberFunctionCallback1 
//   by a numCallees field

// see:
// https://github.com/RobinSchmidt/RS-MET/wiki/The-Modulation-System
// https://github.com/RobinSchmidt/RS-MET/issues/65
// https://github.com/RobinSchmidt/RS-MET/issues/269
// https://github.com/RobinSchmidt/RS-MET/discussions/312

#endif