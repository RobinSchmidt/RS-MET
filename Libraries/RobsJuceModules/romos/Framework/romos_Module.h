#ifndef romos_Module_h
#define romos_Module_h

namespace romos
{

enum connectionKinds
{
  AUDIO,
  EVENT
  // maybe have also:
  // DATA-connections (for example to connect a wavetable to an oscillator)
  // MULTI_CHANNEL-connections
  // ARRAY-connections (for buffered stuf, like FFT, etc.)
  // MATRIX-connections
  // STRING-connections
};

enum connectionDirections
{
  INCOMING,
  OUTGOING
};

class AudioConnection;
class Module;

class ContainerModule;                          // test to satisfy gcc ...seems to work
void  retrieveModuleState(void *moduleAsVoid);  // test to satisfy gcc


//=================================================================================================
// class AudioInputPinData:

/** In its current state, it's only a data class with rudimentary behavior - maybe turn this into 
a proper class later */

class AudioInputPinData  // maybe rename to IncomingConnection
{

public:

  AudioInputPinData();
  AudioInputPinData(const AudioInputPinData &other);
  AudioInputPinData& operator=(const AudioInputPinData &other);


  void setDefaultValue(double newDefaultValue);
  void reset();

  romos::Module *sourceModule;
  double        *outputPointer;     // rename to sourcePointer
  unsigned int  outputIndex;        // rename to sourceOutIndex
  unsigned int  outputFrameSize;    /*<<< distance between successive output frames */
  unsigned int  outputVoiceStride;  /*<<< distance between samples of 2 voices */

  double defaultValue; // value to let outputPointer refer to, when input pin is disconnected

  // \todo: make sourceModule and outputIndex private
  // -> provide accessors getApparent/ActualSourceModule/Index


private:

  void copyDataFrom(const AudioInputPinData &other);

};

//=================================================================================================
// class Module:

/** This is the baseclass for all modules. A module is seen from the outside as a black box that 
has a number of inputs and a number of outputs for both, audio signals and events. The module 
baseclass defines an interface for all concrete modules to adhere to and provides methods for 
keeping track of the wiring between modules.

\todo:
-look into pd, csound, max/msp, reaktor for inspiration
-Module* getClone() operation (virtual)
-array of slave-modules for polyphonic operation
-use unsigned int consistently...or maybe size_t */

//class ModuleTypeInfo;

class Module
{

  friend class ContainerModule;
  friend class ModuleFactory;        // for old creation code
  friend class ModuleFactory;  // for new creation code

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  /** Sets the name of the module. */
  virtual void setModuleName(const std::string& newName);


  virtual void setModuleName(std::string&& newName) noexcept { name = std::move(newName); }
  // optimization for rvalue references, see: https://www.youtube.com/watch?v=xnqTKD8uD64&t=1h9m50s
  // todo: test, benchmark...

  /** Changes the name of one of the input- or output pins. */
  //virtual void setPinName(int kind, int direction, int pinIndex, const rosic::rsString &newName);

  /** Sets the x- and y-coordinates where the module appears on the screen/GUI. The optional 
  boolean parameter sortSiblingsAfterMove decides whether or not the array of siblings of "this" 
  module (that is, the array of child modules of "this" module's parent module) should be sorted 
  after the move. If there is no parent, this has no effect, but if there is a parent, this ensures 
  that after the move, the parent has the child-array sorted in the desired evaluation order as 
  determined from the coordinates of its children. */
  virtual void setPositionXY(int newX, int newY, bool sortSiblingsAfterMove = true);

  /** Switches polyphony for this module on or off. */
  virtual void setPolyphonic(bool shouldBePolyphonic);

  /** Adds a number of pins where the 1st argument specifies the kind, the second the direction, 
  the 3rd the number of pins to add and then follows a variable argument list of type 
  rosic::rsString that define the names of the pins to be created. */
  //void addPins(int kind, int direction, int number, ...);

  /** Connects the input pin of this Module with given inputPinIndex to the output pin of some 
  source-module with the given output pin-index. */
  virtual void connectInputPinTo(int inputPinIndex, Module *sourceModule, 
    int sourceOutputPinIndex);

  /** Disconnects the input pin of this Module with given inputPinIndex. */
  virtual void disconnectInputPin(int inputPinIndex);

  /** Disconnects all our input pins. */
  virtual void disconnectAllInputPins();

  /** Disconnects all our input pins that receive their input from the given 
  sourceMdouleToDisconnect. */
  virtual void disconnectInputPinsWithInputFrom(romos::Module *sourceModuleToDisconnect);

  /** Disconnects all our input pins that receive their input from the given 
  sourceModuleToDisconnect and given outputPinIndex. */
  virtual void disconnectInputPinsWithInputFrom(romos::Module *sourceModuleToDisconnect, 
    int outputPinIndex);

  /** Disconnects the output pin of this Module with given outputPinIndex. */
  virtual void disconnectOutputPin(int outputPinIndex);

  /** Disconnects all our the output pins. */
  virtual void disconnectAllOutputPins();

  /** Updates some shorthand variables in our inputPins array according to some return values of
  the respective pin's source module such as memory location pointer and strides (memory distances) 
  between frames and voices. outputFrameSize is the distance between successive output frames in 
  the source module (this is equal to the source-modules' number of output pins. 
  outputVoiceStride is the memory distance between outputs of 2 voices. */
  virtual void updateInputPointersAndInFrameStrides();
    // rename to updateInPointersAndStrides

  /** A hook function that can be overriden by Module subclasses to set themselves up from a map 
  of key/value pairs. This map is supposed to have been once created by the (also overriden) 
  getState() function. The boolean return value should inform the caller, if setting up the state 
  was successful. */
  virtual bool setState(const std::map<std::string, std::string>& state); // { return true; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry:

  /** Returns the x-coordinate. */
  virtual int getPositionX() const { return x; }

  /** Returns the y-coordinate. */
  virtual int getPositionY() const { return y; }

  /** Informs, whether this module is polyphonic or not. */
  INLINE bool isPolyphonic() const { return polyphonic; }  // we can get rid of the inlinee

  /** Returns voiceAllocator.getNumVoices() when "this" module is and 1 when it is monophonic. */
  virtual int getNumVoices() const;

  /** Informs whether or not the visual rendering of the module needs a header. */
  virtual bool hasHeader() const { return hasHeaderFlag; }

  virtual bool hasEditor() const { return typeInfo->hasEditor; }

  //virtual bool hasTreeNode() const { return typeInfo->hasTreeNode; }



  /** Returns the name. */
  virtual std::string getName() const { return name; }

  /** Returns the identifier of the module class. This is one of the values defined in the 
  moduleIdentifiers enumeration and can be used to infer the kind of the module at runtime (aka 
  runtime type information (RTTI)). Knowing the concrete module (sub)class is necessary for 
  typecasts and total recall. */
  //INLINE int getTypeIdentifierOld() const { return moduleTypeIdentifier; }

  /** A function that can be overriden by Module subclasses to return a state in the form
  of list of key/value pairs. */
  virtual std::map<std::string, std::string> getState() const;
  //{
  //  return std::map<std::string, std::string>(); // return empty map by default
  //}



  INLINE int getTypeId() const { return typeInfo->id; }

  // new:
  virtual std::string getTypeName() const { return typeInfo->fullName; }

#ifdef RS_BUILD_OLD_MODULE_FACTORY

  /** Returns the type name of this module. */
  virtual rosic::rsString getTypeNameOld() const
  {
    return ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(getTypeIdentifierOld());
  }
  // old, deprecated

  /** Returns true when this module is an input module, false otherwise. */
  INLINE bool isInputModule() const { return getTypeIdentifierOld() == ModuleTypeRegistry::AUDIO_INPUT; }

  /** Returns true when this module is an output module, false otherwise. */
  INLINE bool isOutputModule() const { return getTypeIdentifierOld() == ModuleTypeRegistry::AUDIO_OUTPUT; }

  /** Returns true when this module is a container module, false otherwise. */
  INLINE bool isContainerModule() const { return getTypeIdentifierOld() == ModuleTypeRegistry::CONTAINER; }

  /** Returns true when this module is the topl-level module, false otherwise. */
  INLINE bool isTopLevelModule() const { return getTypeIdentifierOld() == ModuleTypeRegistry::TOP_LEVEL_MODULE; }

#else

  // todo:
  // first make it simple:
  INLINE bool isParameterModule()   const { return typeInfo->fullName == "Parameter";   }
  INLINE bool isContainerModule()   const { return typeInfo->fullName == "Container";   }

  INLINE bool isInputModule()       const { return typeInfo->fullName == "AudioInput";  }
  INLINE bool isOutputModule()      const { return typeInfo->fullName == "AudioOutput"; }
  INLINE bool isInputOrOutput()     const { return isInputModule() || isOutputModule(); }

  virtual bool isTopLevelModule()   const { return false; } // overriden in TopLevelModule
  // maybe use the overriding technique for the other module-type inquiries, too



  // later make it fast:
  // INLINE bool isParameter()   const { return typeInfo.id == moduleFactory.getTypeIdForParameter();   }
  // INLINE bool isAudioInput()  const { return typeInfo.id == moduleFactory.getTypeIdForAudioInput();  }
  // INLINE bool isAudioOutput() const { return typeInfo.id == moduleFactory.getTypeIdForAudioOutput(); }
  // INLINE bool isContainer()   const { return typeInfo.id == moduleFactory.getTypeIdForContainer();   }
  // INLINE bool isTopLevel()    const { return typeInfo.id == moduleFactory.getTypeIdForTopLevel();    }

#endif




  /** Returns the parent module of "this" module. */
  virtual romos::ContainerModule* getParentModule() const { return parentModule; }
    // may be a nullptr in cas of the top-level module - maybe also in a transitional state when 
    // the module was not yet given a parent?

  /** Returns the top-level module in the whole module hierarchy. It will resolve "this" for 
  modules that don't have a parent and to a recursive call to getTopLevelModule on the parent for 
  modules that do have a parent. */
  virtual romos::Module* getTopLevelModule();

  /** If "this" module is a child module of some parent module, the function returns the index of 
  "this" module with the parent's array of children (which determines evaluation order), -1 is 
  returned, when this module has no parent. */
  virtual int getIndexWithinParentModule();

  /** Returns the depth of container nesting of the module - atomic modules have a depth of 0, 
  containers one level above them have a depth of 1 and so on. */
  virtual int getContainerNestingDepth() const { return 0; }  
    // get rid of that in the baseclass - move to ContainerModule

  /** Returns the name of one of our pins. */
  virtual rosic::rsString getPinName(int kind, int direction, int pinIndex) const = 0;

  /** Returns the name of one of our audio input pins. Convenience function. */
  inline rosic::rsString getAudioInputPinName(int pinIndex) const
  { return getPinName(AUDIO, INCOMING, pinIndex); }




  /** Returns the number of audio input pins. */
  virtual unsigned int getNumInputPins() const { return (unsigned int)inputPins.size(); }
    // rename to getNumAudioInputs

  /** Returns the data for the input pin with given index. */
  virtual AudioInputPinData getAudioInputPinData(int pinIndex) const;

  /** Returns the number of audio outputs. */
  virtual unsigned int getNumOutputPins() const { return outFrameStride; }
    // redundant with getOutputFrameStride

  /** This is an inlined version of getNumInputs and it is for internal use only - it will return 
  wrong results in case of I/O modules. The internal code knows this and handles it appropriately - 
  client code should only use the virtual getNumInputs function instead. We have this inlined 
  version to make it fast to call from DSP code. But inlining also implies that it can't be virtual
  which is why we can't directly inline getNumInputs. The reason we need to override 
  getNumInputs/Outputs in the I/O modules is because for I/O modules, the numInputs member does 
  actually represent the number of inputs of the parent container-module. This is certainly a major
  quirk, but it is required to make the DSP code fast.  ...maybe someday we find a cleaner solution
  to this problem... */
  INLINE unsigned int getNumInputsInlined() const { return numInputs; }  // not needed anymore?

  /** For internal use only @see getNumInputsInlined */
  //INLINE unsigned int getNumOutputsInlined() const { return outFrameStride; }

  /** Returns a pointer to the audio inputs. */
  //INLINE double** getAudioInputAddress() const { return audioInputs; }

  /** Returns a pointer to the audio outputs. */
  //INLINE double* getAudioOutputAddress() const { return audioOutputs; }


  /** Returns the module that should be used as actual source for the output-signals. Normal 
  modules will return themselves, but some modules act as proxy and return the module that is 
  connected to their inputs instead. When this is called for a container, for example, the 
  container will return itself but instead invoke getConnectableSource */
  //virtual Module* getConnectableSourceModule(Module *attemptedSourceModule);
  // can be used for sender/receiver pairs 


  /** Similar to getConnectableSourceModule but returns the index that is actually be used as 
  source output index, when this is called */
  //virtual int getConnectableSourceOutputPinIndex(int attemtedOutputPinIndex);

  /** Maps a pointer to an apparent source-module/pin-index for a connection to the actual 
  source-module/pin-index that is used in the processing code, resolving all proxy-modules that may
  be in between. It should be called like 
  sourceModule->mapApparentSourceToProcessingSource(sourceModule, sourceOutputPinIndex)
  such that the member-function of "sourceModule" is called with a reference to the 
  "sourceModule"-pointer itself. Normal (non-proxy) modules will do nothing and leave the 
  "sourceModule" pointer-reference at "this" and also leave the pin-index as is. But a container, 
  for example, when acting as source for a connection, will set the sourceModule-pointer reference 
  to the output module that corresponds to the attempted "sourceOutputPinIndex" and set the 
  "sourceOutputPinIndex" to 0 (output modules have only one (invisible) output pin). But that's
  not yet the end of the chain: after setting the modulePointer/pinIndex, the container calls 
  sourceModule->mapApparentSourceToProcessingSource(sourceModule, attemtedOutputPinIndex) again - 
  and because the module-pointer points now to an AudioOutputModule, the overriden implementation 
  there will be called. This in turn will again redirect the pointer, this time to the 
  modulePointer/pinIndex that acts as the source for the output module. This chain of 
  modulePointer/pinIndex manipulations will end when it arrives at an elementary non-proxy module - 
  this is then the actual modulePointer/pinIndex pair that will be used by the signal processing 
  code. */
  virtual void mapApparentSourceToProcessingSource(Module* /*&sourceModule*/,
    int /*&sourceOutputPinIndex*/) { }

  /** Given a modulePointer/pinIndex pair, this function will possibly modify the pair such that it
  represents the module/index pair as apparent on the GUI. This realizes the inverse mapping of 
  mapApparentSourceToProcessingSource - see explanation there for more details. */
  virtual void mapProcessingSourceToSourceApparent(Module* /*&sourceModule*/, 
    int /*&sourceOutputPinIndex*/) { }

  // oevrride both of these in ModuleProxy and ContainerModule



  /** Returns a pointer to the double value that represents the value of the given output-pin in 
  voice 0 and frame 0. */
  virtual double* getOutputPointer(int pinIndex) const; // to be overriden by PointerRedirectModule
    // rename to getAudioOutputAddress

  /** Returns the distance between succsessive sample frames in the output signal buffers. 
  Normally, this is the same as the number of output-pins, but PointerRedirectModules (and 
  subclasses thereof) require some special treatment which is why the function is overriden there. 
  A PointerRedirectModule returns the frame-stride for its connected source-module, by delegating 
  the call to it. */
  virtual int getOutputFrameStride() const; // to be overriden by PointerRedirectModule

  /** Returns the distance between succsesive voices in the output signal buffers. Normally, this 
  is the same as the number of output-pins times the output buffersize (each voice writes into the 
  next voice-buffer), but PointerRedirectModules (and subclasses thereof) require some special 
  treatment @see getOutputFrameStride */
  virtual int getOutputVoiceStride() const;  // to be overriden by PointerRedirectModule

  // \todo: replace these functions with getChannelStride


  /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured
  from our first input memory location which is returned by getAudioInputAddress(). */
  /*
  virtual int getInputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
  {
    rassert( !isInputModule() && !isOutputModule() );  // returns wrong results for I/O modules - override it there, if needed
    return (frameIndex * processingStatus.getNumAllocatedVoices() + voiceIndex) * numInputs + pinIndex;
  }
  */

  /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured
  from our first output memory location which is returned by getAudioOutputAddress().
  ...or getOutputPointer(0)...check this - maybe rename
  */
  virtual int getOutputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
  {
    rassert(!isInputModule() && !isOutputModule());  // returns wrong results for I/O modules - override it there, if needed

    //return (frameIndex * processingStatus.getNumAllocatedVoices() + voiceIndex) * numOutputs + pinIndex; // voice-index last
    return (voiceIndex * processingStatus.getBufferSize() + frameIndex) * outFrameStride + pinIndex;  // frame index last
  }
  // rename to getAudioOutputOffset, write getAudioOutputAddress = base-address + offset

  /** Returns true if the input pin with given index is connected to some output pin, false 
  otherwise. */
  bool isInputPinConnected(int pinIndex) const;

  /** Returns the value to which the signal at the given input pin is set when the pin is not 
  connected to any output. */
  double getInputPinDefaultValue(int pinIndex) const;

  // const char* getInputPinShortName(int pinIndex) { return typeInfo.inputShortNames[pinIndex]; }
  // or:
  // virtual const std::string& getInputPinShortName(int pinIndex) const 
  // { return typeInfo->inputShortNames[pinIndex]; }

  /** Returns true when this module has an AudioOutputModule among its target modules. */
  bool isConnectedToAudioOutput() const;

  /** Returns true when this module has an incoming connection from the given source module. */
  bool hasIncomingConnectionFrom(const romos::Module *sourceModule) const;

  /** Returns true when this module has an incoming connection from the given output pin of the 
  given source module. */
  bool hasIncomingConnectionFrom(const romos::Module *sourceModule, int sourceOutputPinIndex) const;

  /** Returns true when this module has one or more delayed incoming connections. Neede to detect
  feedback loops to switch to sample-wise processing */
  bool hasDelayedIncomingConnection() const;

  /** Returns a vector with all modules that are connected to this module where this module is the 
  source of the connections and the elements of the vector are the targets of the connections. */
  std::vector<romos::Module*> getConnectedTargetModules() const;

  /** Similar to getConnectedTargetModules() but includes only those modules that are connected to 
  the given output pin. */
  std::vector<romos::Module*> getConnectedTargetModulesOfPin(int outputPinIndex) const;


  /** Returns the number of audio connections that are coming in into this module. */
  virtual unsigned int getNumIncomingAudioConnections() const; // { return incomingAudioConnections.size(); }

  /** Returns a pointer to the incoming audio connection with the given index. */
  //virtual AudioConnection* getIncomingAudioConnection(int index) const { return incomingAudioConnections.at(index); }

  /** Returns an array with all incoming audio connections into this module. */
  virtual std::vector<AudioConnection> getIncomingAudioConnections();
    // wouldn't it be more efficient to return an array of pointers? or is that dangerous? or
    // does it not matter bcs this function is not called per sample? ...figure out - add comments
    // actually it seems that an array of AudioConnections is nowhere maintained




  /** Returns an array of pointers to all incoming audio connections to the passed input pin. */
  //virtual std::vector<AudioConnection*> getIncomingAudioConnectionsToPin(int pinIndex);

  /** Returns an array of pointers to all outgoing audio connections. Because modules do not keep 
  track of their outgoing connections, the function must iterate to "this" module's parent's 
  children and for each child, check it's incoming connections if they have "this" module as 
  source-module. */
  virtual std::vector<AudioConnection> getOutgoingAudioConnections();

  /** Returns an array of pointers to all outgoing audio connections that emanate from the passed 
  output pin. */
  virtual std::vector<AudioConnection> getOutgoingAudioConnectionsFromPin(int pinIndex);

  /** Returns true when this module has any incoming audio connections, false otherwise. */
  //virtual bool hasIncomingAudioConnections() const { return !getIncomingAudioConnections().empty(); }

  /** Returns true when this module has any outgoing audio connections, false otherwise. */
  virtual bool hasOutgoingAudioConnections() { return !getOutgoingAudioConnections().empty(); }

  /** Given the passed array of modules, it assigns the reference variables to the minimum and 
  maximum x- and y coordinates of all modules in the array. */
  static void getExtremeCoordinates(std::vector<Module*> &modules, int &xMin, int &yMin, 
    int &xMax, int &yMax);

  /** Given the passed array of modules, it assigns the reference variables to the midpoint of the
  rectangle that is determined by the extreme coordinates (as given by getExtremeCoordinates()). */
  static void getMidpointCoordinates(std::vector<Module*> &modules, int &xMid, int &yMid);

  //-----------------------------------------------------------------------------------------------
  // processing:

  /** Resets the internal state for all voices. */
  virtual void resetStateForAllVoices()
  {
    for(int voiceIndex = 0; voiceIndex < getNumVoices(); voiceIndex++)
      resetVoiceState(voiceIndex);
  }
  // rename to resetAllVoices

  /** Function that is called whenever the module should resets the internal state of the module 
  (for example filter history, oscillator phase, etc) for a particular voice. */
  virtual void resetVoiceState(int voiceIndex)
  {
    rassert(voiceIndex < getNumVoices());
    clearVoiceBuffer(voiceIndex);
  }
  // rename to resetVoice

  /** Zeros the outputs for all voices in all frames. */
  virtual void clearBufferForAllVoices()
  {
    for(int voiceIndex = 0; voiceIndex < getNumVoices(); voiceIndex++)
      clearVoiceBuffer(voiceIndex);
  }

  /** Zeros the outputs for in all frames for a particular voice. */
  virtual void clearVoiceBuffer(int voiceIndex)
  {
    if(audioOutputs == nullptr) // occurs on costruction of AudioOutputModules
      return;

    int voiceStride;
    if(isPolyphonic())
      voiceStride = outFrameStride * processingStatus.getBufferSize();
    else
      voiceStride = 0;
    double *startPointer = audioOutputs + voiceIndex * voiceStride;
    memset(startPointer, 0, outFrameStride * processingStatus.getBufferSize() * sizeof(double));
  }

  // comment seems outdated
  ///** After you have established the input signals to the module by retrieving the input-pointer 
  //via getAudioInputAddress() and assigning values to the returned array/pointer, you should call 
  //this function to trigger the proceesing. When the function returns, the output signals will be 
  //available in the array that you can retrieve via getAudioOutputAddress(). You should make sure to 
  //assign the right number of inputs and retrieve the right number of outputs - what these numbers 
  //are can be inquired via getNumAudioInputs/getNumAudioInputs. */
  INLINE void processSampleFrame()
  {
    processFrame(this, 0);
  }

  /** Similar to processSampleFrame but processes a whole block. */
  INLINE void processBlockOfSamples(int blockSize)
  {
    processBlock(this, 0, blockSize);
  }

  std::vector<AudioInputPinData> inputPins;  // temporarily moved to public for debug
  double *audioOutputs;
  unsigned int outFrameStride, numInputs; // isn't numInputs redundant with inputPins.size()?
  // rename outFrameStride to numOutputs, or numAudioOutputs
  // try to move them into the protected section - currently, doing so gives compiler errors
  // figure out why and try to fix



  //===============================================================================================

protected:


  /** Assigns the two processing function-pointers (for frame-based and block-based processing). 
  Subclasses need to override this and in their overriden function check for the "polyphonic" flag
  in order to assign the monophonic or polyphonic version of each of the functions. */
  virtual void assignProcessingFunctions() = 0;
    // is this still true? should all modules implement all 4 combinations of mono/poly, 
    // block/sample

  /** Allocates the memory area to be used for input/output and as work area by this module. */
  virtual void allocateMemory();

  /** Frees the memory area to be used for input/output and as work area by this module. */
  virtual void freeMemory();

  /** (Re) allocates the memory for the audio output pins and does the necessary callbacks to keep 
  everything informed about these pointers. */
  virtual void allocateAudioOutputs();

  /** Constructor. Protected because instances should be created only through the ModuleFactory. */
  Module(const std::string& name = std::string(), int x = 0, int y = 0, 
    bool polyphonic = false);

  /** Destructor. Protected because instances should be deleted only through the ModuleFactory. */
  virtual ~Module();

  /** To be overriden by concrete subclasses to initialize the names of the pins, allocate 
  resources, etc.. */
  virtual void initialize() = 0;

  /** May be overriden by concrete subclasses to clean up allocated ressources. The baseclass 
  method just calls freeMemory(). */
  virtual void cleanUp();

  //-----------------------------------------------------------------------------------------------
  // data members:


  romos::ContainerModule *parentModule;

  //int  moduleTypeIdentifier; // one of the values defined in moduleIdentifiers - for RTTI and total recall
  // to be deprecated (replaced by pointer to ModuleTypeInfo structure)

  bool polyphonic;           // flag to indicate that this module is polyphonic
  int  x, y;                 // position on the GUI block-diagram


  std::string name;
  //rosic::rsString name;      // use std::string
  // have shortName, longName, description, write a function getLongName etc. - inside this
  // function, check, if the variables here are empty and if so, show the long/short/etc names
  // from the typeInfo pointer

  // hmm...maybe it was not such a good idea to switch to std::string
  // https://stackoverflow.com/questions/5058676/stdstring-implementation-in-gcc-and-its-memory-overhead-for-short-strings
  // http://jovislab.com/blog/?p=76
  // https://www.oreilly.com/library/view/optimized-c/9781491922057/ch04.html
  // maybe i really should do a measurement (sizeof)

  ModuleTypeInfo *typeInfo = nullptr; 
  // not yet used - should later be set on construction by the factory and replace the
  // moduleTypeIdentifier variable

  bool  hasHeaderFlag;
  // determines, if the visual rendering needs a header - actually a GUI-thing that does not 
  // really belong here - maybe write a function instead that switches on the typeCode or similar
  // remove - use info from typeInfo object

  // function pointers for audio processing:
  void (*processFrame) (Module *module, int voiceIndex);
  void (*processBlock) (Module *module, int voiceIndex, int blockSize);

  // friends:
  friend void processContainerMixedMonoPoly(Module *module, int voiceIndex);
  friend void processContainerAllMono(Module *module, int voiceIndex);
  friend void processContainerAllPoly(Module *module, int voiceIndex);
  friend void processContainerMixedMonoPolyBlock(Module *module, int voiceIndex, int blockSize);
  friend void processContainerAllMonoBlock(Module *module, int voiceIndex, int blockSize);
  friend void processContainerAllPolyBlock(Module *module, int voiceIndex, int blockSize);
  friend void processContainerBlockFrameWiseMono(Module *module, int blockSize);
  friend void processContainerBlockFrameWisePoly(Module *module, int blockSize);
  friend void processContainerBlockFrameWiseMixed(Module *module, int blockSize);

  friend void writeModuleStateToConsole(void *module, bool waitForKeyAfterOutput);
  friend void retrieveModuleState(void *module); // for debugging



private:

  // copy and assignment - not possible:
  Module(const Module& /*other*/) {}
  Module& operator=(const Module& /*other*/) { return *this; }

};

/** Decides for two modules (given by their pointers), if the first is considered to become before 
the second or not. This is needed to sort module-pointers according to their desired evaluation 
order which shall be governed by the x- and y coordinates of the respective modules. Modules that 
are to the left of other modules (that is, with lower x-coordinate) are evaluated before those 
other modules. Inside a column with equal x-values, the upper modules (those with lower y-values) 
are evaluated before others. When two modules have equal x- and y coordinates, the evaluation order
is undefined - this should be prevented by presenting the user a warning in such cases. Or the GUI
code should prevent the user from placing a module on top of another. */
bool modulePointerLessByXY(const romos::Module *left, const romos::Module *right);

/** Returns true when the the passed vector of modules contains one (or more) modules of the given 
type (the code is interpreted according to the ModuleTypeRegistry::moduleTypeCodes enumerator. */
//bool containsModuleOfType(const std::vector<romos::Module*> &modules, int typeCode);
// old

bool containsModuleOfType(const std::vector<romos::Module*> &modules, 
  const std::string& fullTypeName);
// new


/** Similar to modulePointerLessByXY but considers y-coordinate as more important than 
x-coordinate. */
//bool modulePointerLessByYX(romos::Module *modulePointer1, romos::Module *modulePointer2);

/** A helper function for degubbing purposes - writes the state (name, numChildren, parent, 
numIns/Outs, values of ins/outs, etc.) of the module to the standard output. The module should be 
passed as void pointer so we can call it from the process functions. The second parameter specifies
whether the program should pause and wait for the user to hit a key until it continues. */
void writeModuleStateToConsole(void *module, bool waitForKeyAfterOutput);

/** A function that can be called when things occur that shouldn't. It will open a message box with 
the error text and give the user the option to write the message into a logfile. ....later. */
void triggerRuntimeError(const char *errorMessage);


} // end namespace romos

#endif
