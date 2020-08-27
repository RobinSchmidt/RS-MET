#ifndef romos_AudioConnection_h
#define romos_AudioConnection_h

/** This class represents a connection from one Module to another over which audio signals are
transmitted. 

seems like connectiosn objects are just used temporarily to conveniently pass connection data 
around

*/

class Module;

class AudioConnection
{

  friend class ContainerModule;

public:

  /** Enumeration of some criterions under which two connections can be considered as "alike" - 
  used by areConnectionsAlike. */
  enum alikenessCriterions
  {
    INPUT_PINS_EQUAL,
    OUTPUT_PINS_EQUAL,
    ALL_PINS_EQUAL
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. Creates an audio connection ('wire') between some output slot of the source 
  module and an input slot of the target module. 
  
  obsolete?:
  The target module will establish its total input  as the weighted sum over all incoming 
  signals ...(i think that is not true anymore). */
  AudioConnection(Module *sourceModule = nullptr, int outputIndex = 0, 
    Module *targetModule = nullptr, int inputIndex = 0);

  /** Destructor. 
  
  obsolete?:
  It will take care to remove itself from the source- and target module's array of 
  outgoing/incoming connections. */
  ~AudioConnection();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Updates the pointer to the output address of the source module. This pointer may change when 
  the source module creates a new output slot in which case this function should be called in order
  to have the connection the correct pointer again. */
  //void updateSourcePointer()
  //{ sourcePointer = sourceModule->getAudioOutputAddress() + outIndex; }

  /** Updates the pointer to the input address of the target module. @see setSourcePointer */
  //void updateTargetPointer() 
  // { targetPointer = targetModule->getAudioInputAddress() + inIndex; }

  /** Sets source- and target modules to NULL and the pin-indices to 0 */
  void resetToNull();

  /** Sets the module which is the source of this connection. */
  void setSourceModule(Module *newSourceModule);

  /** Sets the module which is the target of this connection. The function will NOT take care to 
  remove "this" connection from the array of incoming connections of the old target module and add 
  it to the incoming connections of the new target module. */
  void setTargetModule(Module *newTargetModule);

  /** Sets the index of the output pin of the source module that is connected by this 
  connection. */
  void setSourceOutputIndex(int newIndex);

  /** Sets the index of the input pin of the target module that is connected by this connection. */
  void setTargetInputIndex(int newIndex);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the source module of this connection. */
  INLINE Module* getSourceModule() const { return sourceModule; }

  /** Returns the targetModule of this connection. */
  INLINE Module* getTargetModule() const { return targetModule; }

  /** Returns the index of the output pin of the source module that is connected by this 
  connection. */
  INLINE int getSourceOutputIndex() const { return outIndex; }

  /** Returns the index of the input pin of the target module that is connected by this 
  connection. */
  INLINE int getTargetInputIndex() const { return inIndex; }

  /** Returns a pointer to the memory location that holds the source signal of this connection. */
  //INLINE double *getSourcePointer() const { return sourcePointer; }

  /** Returns a pointer to the memory location that holds the target signal of this connection. */
  //INLINE double *getTargetPointer() const { return targetPointer; }

  /** Returns true when both sourceModule and targetModules are NULL. */
  virtual bool isNull()
  {
    return sourceModule == NULL && targetModule == NULL;
  }

  /** Returns true if this connection involves the passed Module (either as source or as target), 
  false otherwise. */
  virtual bool involvesModule(Module *moduleToCheckFor) const
  {
    return moduleToCheckFor == sourceModule || moduleToCheckFor == targetModule;
  }

  /** Returns true when "this" connection has the same source pin (i.e. the same source-module and 
  source output-pin index) as the passed connection. */
  virtual bool hasSameSourcePinAs(AudioConnection *otherConnection);

  /** Returns true when "this" connection has the same target pin (i.e. the same target-module and 
  target input-pin index) as the passed connection. */
  virtual bool hasSameTargetPinAs(AudioConnection *otherConnection);

  /** Returns true when "this" connection connects the same pins as the passed connection. */
  virtual bool connectsSamePinsAs(AudioConnection *otherConnection);

  /** Returns true when this connection contains an implicit delay. This occurs whenever the target
  module of the connection is evaluated before the source module of the connection, or the source 
  and target modules are the same (the connection connects an output of some module with an input 
  of the same module.) Such an implicit delay is sometimes inavoidable in graphs that contain 
  feedback. */
  virtual bool hasImplicitDelay() const;

  /** Returns whether or not the two passed connections are considered "alike" as defined by some 
  criterion which should be one of the values defined in the alikenessCriterions enumeration. */
  static bool areConnectionsAlike(AudioConnection *connection1, AudioConnection *connection2, 
    int criterion);

  /** Returns true when (1) both arrays are of the same length and (2) for each connection in the 
  first array, there is a connection in the second array that is "alike" and vice versa (according
  to some criterion @see alikenessCriterions). For example, this function is useful to figure out 
  if some input pin receives the same signal as some other input pin, in which case the arrays of 
  incoming connections would be "alike" for both pins according to the INPUT_PINS_EQUAL criterion.
  The function does not assume any ordering of the arrays - if there's a connection at position 7 
  in array 1 that is alike a connection at position 4 in array 2, the function will figure it out. 
  In the example above, this is to say, both pins receive the same signal regardless in which 
  order they sum over their incoming wires.  */
  static bool areConnectionArraysAlike(std::vector<AudioConnection*> connections1, 
    std::vector<AudioConnection*> connections2, int criterion);

  /** Compares two AudioConnection of equality. Two connections are euqal if their source- and 
  target modules and indices match. */
  bool operator==(const AudioConnection& other) const
  {
    if(sourceModule == other.sourceModule && targetModule == other.targetModule
      && outIndex     == other.outIndex     && inIndex      == other.inIndex)
      return true;
    else
      return false;
  }


protected:

  Module *sourceModule, *targetModule;
  int    outIndex, inIndex;

  //double *sourcePointer, *targetPointer; 
  // maybe get rid of these and compute them on the fly - provide (inlined) functions 
  // getSource/TargetPointer()


  /** Structure to be used for the breakpoints. */
  struct ConnectionNode
  {
    int x, y;
  };

  //rosic::Array<ConnectionBreakpoint> breakpoints; 
  // when we indeed use such an array here, we must take care to override our assignment operator 
  // and copy constructor to create a deep copy of such an array
  // better - just store a pointer to such an array here - when creating temporary connection 
  // objects, we don't need to copy everything

  friend void processContainerInSampleFrames(romos::Module *module, double *ins, double *outs, 
    int voiceIndex);
  friend void processContainerInBlocks(romos::Module *module, double *ins, double *outs, 
    int voiceIndex);
  //friend void processContainer(romos::Module *module, double *ins, double *outs, int voiceIndex);

};


/*
INLINE void AudioConnection::accumulateMonophonicToMonophonic()
{
  *(getTargetPointer()) += *(getSourcePointer());
}
*/

/*
INLINE void AudioConnection::accumulateMonophonicToPolyphonic()
{
  int numPlayingVoices           = ProcessingStatus::getSoleInstance()->getNumPlayingVoices();
  const int *playingVoiceIndices = ProcessingStatus::getSoleInstance()->getPlayingVoiceIndices();
  int inStride = targetModule->getInputVoiceStride();
  for(int v=0; v<numPlayingVoices; v++)
    *(getTargetPointer() + inStride*playingVoiceIndices[v]) += *(getSourcePointer());
}
*/

/*
INLINE void AudioConnection::accumulatePolyphonicToMonophonic()
{
  int numPlayingVoices           = ProcessingStatus::getSoleInstance()->getNumPlayingVoices();
  const int *playingVoiceIndices = ProcessingStatus::getSoleInstance()->getPlayingVoiceIndices();
  int outStride = sourceModule->getOutputVoiceStride();
  for(int v=0; v<numPlayingVoices; v++)
    *(getTargetPointer()) += *(getSourcePointer() + outStride*playingVoiceIndices[v]);  // shouldnt it be numOuts*allocatedBlockSize?
}
*/

/*
INLINE void AudioConnection::accumulatePolyphonicToPolyphonic()
{
  int numPlayingVoices           = ProcessingStatus::getSoleInstance()->getNumPlayingVoices();
  const int *playingVoiceIndices = ProcessingStatus::getSoleInstance()->getPlayingVoiceIndices();
  int inStride  = targetModule->getInputVoiceStride();
  int outStride = sourceModule->getOutputVoiceStride();
  for(int v=0; v<numPlayingVoices; v++)
    *(getTargetPointer() + inStride*playingVoiceIndices[v]) += *(getSourcePointer() + outStride*playingVoiceIndices[v]);
}
*/

#endif
