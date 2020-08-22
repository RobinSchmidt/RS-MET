#ifndef romos_ContainerModule_h
#define romos_ContainerModule_h

namespace romos
{

// declarations to satisfy gcc:
void processContainerBlockFrameWiseMixed(Module *moduleAsVoid, int blockSize);
void processContainerBlockFrameWiseMono(Module *moduleAsVoid, int blockSize);
void processContainerBlockFrameWisePoly(Module *moduleAsVoid, int blockSize);
void processContainerMixedMonoPoly(Module *module, int voiceIndex);
void processContainerAllMono(Module *module, int voiceIndex);
void processContainerAllPoly(Module *module, int voiceIndex);
void processContainerMixedMonoPolyBlock(Module *module, int voiceIndex, int blockSize);
void processContainerAllMonoBlock(Module *module, int voiceIndex, int blockSize);
void processContainerAllPolyBlock(Module *module, int voiceIndex, int blockSize);

/** This is the baseclass for all modules that have embedded child-modules (and may also have a 
variable number of inputs and/or outputs). This class is where the meat of the modular capabilities
is implemented.

\todo:
-use the command pattern to add/remove child-modules and maybe various other potentially damgeful 
 actions to enable undo/redo */

class AudioInputModule;
class AudioOutputModule;

class ContainerModule : public Module
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:



  //-----------------------------------------------------------------------------------------------
  // Setup:

  /** Switches polyphony for this module on or off. For containers, this means that the poly-flag 
  will be set for the container itself as well as its input- and output modules. */
  virtual void setPolyphonic(bool shouldBePolyphonic);

  /** Switches polyphony for this module and recursively for all its child modules on or off. */
  virtual void setPolyphonicRecursively(bool shouldBePolyphonic);

  /** Overriden from Module to additionaly set up the pin-data in the input-module that corresponds 
  to the given pin-index. */
  virtual void connectInputPinTo(int inputPinIndex, Module *sourceModule, 
    int sourceOutputPinIndex) override;

  /** Overriden from Module to additionaly reset the pin-data in the input-module that corresponds 
  to the given pin-index. */
  virtual void disconnectInputPin(int inputPinIndex) override;

  /** Adds an audio input to this module. A name for the pin (and module) can optionally be passed. 
  If empty, the function will assign a default name. The return value is a pointer to the added 
  module. */
  virtual Module* addAudioInputModule(std::string name = std::string(), int x = 1, int y = 1, 
    bool sortModuleArrayAfterInsertion = true);

  /** Adds an audio input to this module. A name for the pin can optionally be passed. */
  //void addAudioInput(const rosic::rsString &pinName = rosic::rsString());

  /** Adds an audio output to this module. A name for the pin (and module) can optionally be 
  passed. If empty, the function will assign a default name. The return value is a pointer to the 
  added module. */
  virtual Module* addAudioOutputModule(std::string name = std::string(), int x = 1, int y = 1, 
    bool sortModuleArrayAfterInsertion = true);


  //virtual Module* addChildModule(NULL, ModuleTypeRegistry::ADD, rosic::rsString(""),  10,  2, false);


  /** Adds a child-module to this one. This object takes over the responsibility to delete the 
  object such that sub modules are automatically deleted when their parent gets deleted. The 
  optional parameter determines whether or not the arrays of embedded modules shall be sorted after
  insertion. Normally, you will want to sort in order to always have a well defined evaluation 
  order of the child modules - however, when adding a bunch of child modules at once, you might 
  want to pass false (inside the adding-loop) and afterwards call sortChildModuleArray manually in 
  order to avoid repeated sorting. The return value is a pointer to the added module.  */
  virtual Module* addChildModule(Module *moduleToAdd, bool sortModuleArraysAfterInsertion = true);

  /** Adds a child- or I/O module of the kind given by the identifier 
  (@see romos::moduleIdentifiers) at the given coordinates. 
  @see addChildModule(Module*, bool, bool). The return value is a pointer to the added module. */
  //virtual Module* addChildModule(int moduleIdentifier, rosic::rsString name = rosic::rsString(), 
  //  int x = 0, int y = 0, bool polyphonic = false, bool sortChildModulesAfterInsertion = true);
  // old version - deprecated

  // new version:
  virtual Module* addChildModule(const std::string& fullTypeName, const std::string& name = "", 
    int x = 0, int y = 0, bool polyphonic = false, bool sortChildModulesAfterInsertion = true);


  /** Deletes a child-module from this one. */
  virtual void deleteChildModule(Module *moduleToDelete, 
    bool updateHasDelayedConnectionFlag = true);

  /** Deletes all child modules of this module. */
  //virtual void deleteAllChildModules();

  /** Deletes all child modules of this module. */
  virtual void deleteAllChildModules();
  //virtual void deleteAllEmbeddedModules();

  /** Deletes a bunch of embedded module from this one (these can be child modules or in-/out 
  modules). */
  virtual void deleteModules(std::vector<Module*> modulesToDelete);

  /** Sets the polyphony for the passed array of modules. */
  virtual void setPolyphonyForModules(std::vector<Module*> modules, bool shouldBePolyphonic);

  /** Puts the passed array of modules into a container. It will also take care of keeping the 
  connections intact by equipping the to-be-created container with an appropriate number of in/out
  modules and connecting them according to the connectivity of the to-be-containerized modules. It
  returns a pointer to the just created container. */
  virtual ContainerModule* containerizeModules(std::vector<Module*> modulesToContainerize);

  /** Extracts all container modules in the passed array and puts their content modules as direct 
  child-modules of "this" one. */
  virtual void unContainerizeModules(std::vector<Module*> modulesToUnContainerize);

  /** Unpacks the modules that are contained inside the passed container and makes them direct 
  child modules of "this" one. */
  virtual void unContainerize(ContainerModule *container);

  /** Establishes an audio connection ('wire') between some output slot of the source module and an
  input slot of the target module. The target module will establish its total input as the weighted 
  sum over all incoming signals. */
  virtual void addAudioConnection(Module *sourceModule, int outputIndex, Module *targetModule,
    int inputIndex);

  /** Adds an audio connection that was created outside this object. We assume that the source and 
  target modules of the passed connection are actually child modules of "this" module. This module 
  takes over ownership of the connection, so don't delete it, once added. */
  virtual void addAudioConnection(AudioConnection *connectionToAdd);

  /** Deletes a connection between the passed source and target module, if such a connection is 
  present. The return value informs about whether such a connection was present. */
  virtual bool deleteAudioConnection(Module *sourceModule, int outputIndex, Module *targetModule, 
    int inputIndex);

  /** Deletes the passed audio connection. */
  virtual void deleteAudioConnection(AudioConnection connectionToDelete);

  /** Deletes a bunch of audio connections at once. */
  virtual void deleteAudioConnections(std::vector<AudioConnection> connectionsToDelete);

  /** Minimizes the number of input pins by investigating which pins are superfluos, reconfiguring 
  the connections accordingly and deleting the now obsolete pins. */
  virtual void minimizeNumberOfAudioInputs();

  /** Sorts the array of the child modules according to their desired evaluation order which is
  determined form the coordinates of the modules. */
  virtual void sortChildModuleArray();


  /** Returns true when this container uses the proxy versions of the input/output modules.
  \todo implement this proxy-stuff
  These are more efficient but can't be used in certain circumstances ...move to class ProxyInputModule/OutputModule. */
  virtual bool usesInOutProxyModules()
  {
    return false; // preliminary
  }


  virtual void mapApparentSourceToProcessingSource(Module * &sourceModule, int &sourceOutputPinIndex);
  virtual void mapProcessingSourceToSourceApparent(Module * &sourceModule, int &sourceOutputPinIndex);


  /** Sorts the passed array of modules such that modules with lower lower x-coordinate come before 
  those with higher x-coordinates. When x-coordinates are equal, the y-coordinate determines the 
  ordering in the same way. If both coordinates are equal for tow modules, their order is undefined
  (->avoid this situation in the first place). */
  static void sortModuleArrayByCoordinates(std::vector<romos::Module*> &modulesToSort);

  /** Removes modules from the passed vector reference that match the given type-code. */
  //static void removeModulesOfType(std::vector<romos::Module*> &modules, int typeCodeToRemove);
  static void removeModulesOfType(std::vector<romos::Module*> &modules, const std::string& typeName);


  // replaceChildModule, duplicateChildModule, copySelectedCildModulesIntoClipboard,
  // cutSelectedCildModulesIntoClipboard, pasteClipboardContent(float x, float y),
  // clearClipboardContent, getState/setState,
  // removeSelectedChildModules, disconnectAllInputs(Module*), disconnetctAllOupts(Module*),
  // diconnectAllInputsAndOutputs(Module*)

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Overloads the inherited function from Module in order to return the address of the memory of 
  one of the embedded input-modules The inherited getAudioInputAddress(void) doesn't make sense 
  anymore for containers because the memory for different pins/channels isn't contiguous in 
  containers. */
  //double* getAudioInputAddress(int pinIndex) const { return getAudioInputModule(pinIndex)->getAudioInputAddress(); }
  //double* getAudioOutputAddress(int pinIndex) const { return getAudioOutputModule(pinIndex)->getAudioOutputAddress(); }


  /** Returns a pointer to one of the input modules. */
  virtual AudioInputModule* getAudioInputModule(int index) const;

  /** Returns a pointer to one of the output modules. */
  virtual AudioOutputModule* getAudioOutputModule(int index) const;

  /** Assuming that the passed module is one of our output modules, this function returns its 
  pin-index as visible from outside the container. */
  int getOutputPinIndexOf(AudioOutputModule *outputModule) const;

  /** Assuming that the passed module is one of our input modules, this function returns its 
  pin-index as visible from outside the container. */
  int getInputPinIndexOf(AudioInputModule *inputModule) const;


  /** Returns the number of child modules. */
  virtual unsigned int getNumChildModules() const { return (unsigned int)childModules.size(); }

  /** Returns a pointer to one of the child modules. */
  virtual Module* getChildModule(int index) const;

  /** Returns true when all child-modules (including I/O modules) are monophonic, false 
  otherwise. */
  virtual bool areAllChildModulesMonophonic() const;

  /** Returns true when all child-modules (including I/O modules) are polyphonic, false 
  otherwise. */
  virtual bool areAllChildModulesPolyphonic() const;

  /** Returns the depth of container nesting of the module - atomic modules have a depth of 0, 
  containers one level above the have a depth of 1 and so on. */
  virtual int getContainerNestingDepth() const;

  /** Returns an array with pointers to our child modules. */
  virtual std::vector<Module*> getChildModules() const { return childModules; }

  /** Returns the index of the passed module inside our array of child modules. The passed module
  is thus assumed to point to one of our children - if it doesn't, then you are probably doing 
  something wrong and the function will trigger a debug-break and return -1. */
  virtual int getIndexOfChildModule(romos::Module *moduleToFindIndexFor);

  /** Returns an array with pointers to our child modules (excluding I/O modules) .*/
  virtual std::vector<romos::Module*> getNonInOutChildModules() const;

  /** Returns an array with pointers to our child modules that match the given type-identifier 
  which should be one of the values enumerated in ModuleTypeRegistry::moduleIdentifiers. */
  //virtual std::vector<romos::Module*> getChildModulesWithTypeOld(int typeIdentifier) const;
    // deprecated

  virtual std::vector<romos::Module*> getChildModulesWithTypeId(int typeIdentifier) const;

  virtual std::vector<romos::Module*> getChildModulesWithType(const std::string& type) const;



  /** Returns a vector of modules that are connected as target-modules to the given 
  source-module. */
  std::vector<romos::Module*> getConnectedTargetModulesOf(const romos::Module* sourceModule) const;

  /** Returns a vector of modules that are connected as target-modules to the given output pin of 
  the given source-module. */
  std::vector<romos::Module*> getConnectedTargetModulesOf(const romos::Module* sourceModule, 
    int outputPinIndex) const;

  /** Overriden here because it does not equal the outFrameStride member variable for 
  containers. */
  virtual unsigned int getNumOutputPins() const
  {
    //return (unsigned int)getChildModulesWithTypeOld(ModuleTypeRegistry::AUDIO_OUTPUT).size();
    // old version

    return (unsigned int)getChildModulesWithType("AudioOutput").size(); // new version doesn't work because of wrong id?
  }

  /** Overriden her because for containers, this does not necessarily equal the number of output 
  pins (only in the special case of one output pin) */
  virtual double* getOutputPointer(int pinIndex) const
  {
    return audioOutputs + pinIndex * getRequiredOutputBufferSizePerPin();
  }

  /** Returns the name of one of our pins. Overriden from Module in order to retrun the name of one 
  of our I/O modules. */
  virtual rosic::rsString getPinName(int kind, int direction, int pinIndex) const;

  /** Assuming that the passed module is one of our input- or output modules, the function assigns 
  the passed reference variables to the kind, direction and index of the pin that corresponds to 
  that module. When you pass a module that is none of our I/O modules, the variables will all be 
  assigned to -1. */
  //virtual void getPinDataForModule(Module *module, int &kind, int &direction, int &pinIndex) const;

  /** Returns the name of one of our pins. */


  /** Returns true if any of the connections inside this container is a connection with implicit 
  delay, false otherwise. Note: in order to determine whether we can do block processing, this info
  is not sufficient - we must also ensure that the parent does not contain a feedback loop that 
  contains "this" module. */
  //virtual bool containsConnectionsWithImplicitDelay() const;

  /** Returns true if the position given by the coordinates x, y is occupied by some module, false
  otherwise. */
  virtual bool isPositionOccupied(int &x, int &y) const;

  /** Given coordinates x, y, the function checks if the position is already occupied by some 
  module and if so, re-assigns them to a nearby a non-occupied position. If the position is not 
  occupied, it does nothing. The functin is useful for ensuring that no position is occupied by 
  more than one module (in which case evaluation order would be undefined) */
  virtual void getNonOccupiedPositionNear(int &x, int &y) const;

  //-----------------------------------------------------------------------------------------------
  // misc:

  /** Callback that is supposed to be called from child-modules when they change theri polyphony 
  setting. This information is used here in order to apply optimized processing functions when all 
  children have the same polyphony setting. */
  virtual void childPolyphonyChanged(Module *childThatHasChangedPolyphony);

  /** Overriden to call resetVoiceState on all our child-modules. */
  virtual void resetVoiceState(int voiceIndex);

  /** Allocates the memory area to be used for input/output and as work area by this module. */
  //virtual void allocateMemory();

  /** Frees the memory area to be used for input/output and as work area by this module. */
  virtual void freeMemory();

  /** (Re) allocates the memory for the audio input pins and does the necessary callbacks to keep 
  everything informed about these
  pointers. */
  //virtual void allocateAudioInputs();

  /** Overriden from Module to update the corresponding variables in the input modules, too. */
  virtual void updateInputPointersAndInFrameStrides();

  /** The input/output modules store pointers to memory locations which are allocated and 
  deallocated here - this function is used to pass the respective pointers to our input modules. */
  virtual void updatePointersInInputModules();


  /** Returns the required siize for the output buffer in number-of-values. */
  virtual int getRequiredOutputBufferSize() const;

  /** Returns the required size of the memory area that is required for each output pin in 
  number-of-values (each value is a double) */
  virtual int getRequiredOutputBufferSizePerPin() const; // maybe move to Module baseclass

  /** (Re) allocates the memory for the audio output pins and does the necessary callbacks to 
  keep everything informed about these pointers.  */
  virtual void allocateAudioOutputs();

  /** The input/output modules store pointers to memory locations which are allocated and 
  deallocated here - this function is used to pass the respective pointers to our output 
  modules. */
  virtual void updatePointersInOutputModules();

  /** Assigns the processing functions. */
  virtual void assignProcessingFunctions();

  /** A callback function that is called from one of our child-modules to notify us that the 
  respective child-module has re-allocated the memory for its output signals. This will invalidate
  all input-pointers inside our other child-modules that are connected to these outputs, so we 
  respond to this callback with letting all child-modules update their input pointers. */
  virtual void outputsWereReAllocated(Module *moduleThatHasReAllocated);

  //-----------------------------------------------------------------------------------------------
  // processing:

  /** After you have established the input signals to the module by retrieving the input-pointer 
  via getAudioInputAddress() and assigning values to the returned array/pointer, you should call 
  this function to trigger the proceesing. When the function returns, the output signals will be 
  available in the array that you can retrieve via getAudioOutputAddress(). You should make sure to
  assign the right number of inputs and retrieve the right number of outputs - what these numbers 
  are can be inquired via getNumAudioInputs/getNumAudioInputs. */
  /*
  INLINE void processSampleFrame(int voiceIndex)
  {
    // invoke our function-pointer member:
    //process(getAudioInputAddress(), getAudioOutputAddress(), getNumAudioInputs(), getNumAudioOutputs(), data);
    process(this, NULL, NULL, voiceIndex);
  }
  */

  /** Similar to processSampleFrame but processes a whole block. */
  /*
  INLINE void processBlockOfSamples(int voiceIndex, int blockSize)
  {
    processBlock(this, blockSize);
  }
  */


protected:

  //-----------------------------------------------------------------------------------------------
  // internal functions:

  /** Checks whether or not the passed module is among our direct child modules (including in/out 
  modules). */
  virtual bool hasAsDirectlyEmbeddedModule(Module *moduleToSearchFor);

  /** Checks whether or not the passed module is among our descendants - that is, our 
  child-modules, their children and so on. */
  virtual bool hasAsDescendant(Module *moduleToSearchFor);

  /** Updates our hasDelayedConnections member variable - should be called from all operations that 
  could cause a change in this flag, for example adding connections, insertion/removal/re-ordering
  (sorting) of child-modules, etc. */
  virtual void updateHasDelayedConnectionsFlag();

  /** Constructor. */
  ContainerModule(const std::string& name = std::string(), int x = 0, int y = 0, 
    bool polyphonic = false);

  /** Destructor. Deletes all sub modules as well. */
  virtual ~ContainerModule();

  /** To be overriden by concrete subclasses to initialize the names of the pins, allocate 
  ressources, etc.. */
  virtual void initialize();

  /** To be overriden by concrete subclasses to clean up allocated ressources. */
  virtual void cleanUp();

  //-----------------------------------------------------------------------------------------------
  // data-members:

  std::vector<Module*> childModules;

  double *tmpOutFrame; // needed, when processing a block in frames

  //bool   canDoBlockProcessing;
  bool hasDelayedConnections;

  // friend functions:

  friend void processContainerMixedMonoPoly(Module *module, int voiceIndex);
  friend void processContainerAllMono(Module *module, int voiceIndex);
  friend void processContainerAllPoly(Module *module, int voiceIndex);

  friend void processContainerMixedMonoPolyBlock(Module *module, int voiceIndex, int blockSize);
  friend void processContainerAllMonoBlock(Module *module, int voiceIndex, int blockSize);
  friend void processContainerAllPolyBlock(Module *module, int voiceIndex, int blockSize);


  friend void processContainerBlockFrameWiseMono(Module *moduleAsVoid, int blockSize);
  friend void processContainerBlockFrameWisePoly(Module *moduleAsVoid, int blockSize);
  friend void processContainerBlockFrameWiseMixed(Module *moduleAsVoid, int blockSize);

  friend class ContainerModuleTypeInfo;
  friend class ModuleFactory; // old

private:

  // copy and assignment - not possible:
  ContainerModule(const ContainerModule& /*other*/) {}
  ContainerModule& operator=(const ContainerModule& /*other*/) { return *this; }

};
class ContainerModuleTypeInfo : public ModuleTypeInfo
{
public:
  ContainerModuleTypeInfo() {
    shortName    = "Container";
    fullName     = "Container";
    description  = "Wraps a bunch of interconnected modules into a single module";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new ContainerModule; };
  }
};


} // end namespace romos

#endif
