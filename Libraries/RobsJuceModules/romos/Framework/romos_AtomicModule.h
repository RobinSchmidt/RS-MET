#ifndef romos_AtomicModule_h
#define romos_AtomicModule_h

namespace romos
{

// maybe rename file to ModuleBaseClasses - we have more specific baseclasses here too
// or better: move those subclasses into some file in the "Modules" folder

//=================================================================================================
// class AtomicModule:

/** This is the baseclass for atomic modules (except for some very special ones). Atomic modules 
differ from composite modules in the way they deal with connections. Atomic modules deal with their 
connections directly while composite modules delegate their connections to their I/O 
child-modules. */

class AtomicModule : public Module
{

  friend class ContainerModule;

public:

  //-----------------------------------------------------------------------------------------------
  // setup of pins:



  //-----------------------------------------------------------------------------------------------
  // inquiry about pins:

  /** Returns the name of one of our pins. */
  virtual rosic::rsString getPinName(int kind, int direction, int pinIndex) const;


  // these two are already defined in Module:
  /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured 
  from our first input memory location which is returned by getAudioInputAddress(). */
  /*
  INLINE int getInputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
  {
    return (frameIndex * ProcessingStatus::getSoleInstance()->getNumAllocatedVoices() 
    + voiceIndex) * getNumAudioInputs() + pinIndex;
  }
  */
  /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured 
  from our first output memory location which is returned by getAudioOutputAddress(). */
  /*
  INLINE int getOutputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
  {
    return (frameIndex * ProcessingStatus::getSoleInstance()->getNumAllocatedVoices() 
    + voiceIndex) * getNumAudioOutputs() + pinIndex;
  }
  */


protected:

  /** Creates an initial number of input pins. To be called from initialize only. When more pins 
  should be added later, use addAudioInput. */
  //virtual void initInputPins(int numberOfPins, const char*, ...);
  // deprecated

  virtual void initInputPins(const std::vector<std::string>& pinNames);
  // new version - avoids the error-prone vararg mechanism

  /** Creates an initial number of output pins. To be called from initialize only. When more pins 
  should be added later, use addAudioOutput. */
  //virtual void initOutputPins(int numberOfPins, const char*, ...); 
  // deprecated

  virtual void initOutputPins(const std::vector<std::string>& pinNames);
  // new version

  /** Adds an audio input. A name for the pin can optionally be passed. */
  virtual void addAudioInput(const char* shortName = "", const char* longName = "", 
    const char* description = "");

  /** Adds an audio output. A name for the pin can optionally be passed. */
  virtual void addAudioOutput(const char* pinName = "");
   // these functions are in AtomicModule and not in Module because in containers the inputs and
   // outputs are actually modules

  /** Deletes the audio input with the given index. */
  virtual void deleteAudioInput(int index);

  /** Deletes the audio output with the given index. */
  virtual void deleteAudioOutput(int index);

  /** Constructor. Protected because instances should be created only through the ModuleFactory. */
  AtomicModule(const std::string& name = std::string(), int x = 0, int y = 0, 
    bool polyphonic = false);

  /** Destructor. Protected because instances should be deleted only through the ModuleFactory. */
  virtual ~AtomicModule();

  //-----------------------------------------------------------------------------------------------
  // data members:

  std::vector<rosic::rsString> audioInputNames;   // rename to audioInputShortNames
  std::vector<rosic::rsString> audioInputLongNames;
  std::vector<rosic::rsString> audioInputDescriptions;

  std::vector<rosic::rsString> audioOutputNames;
  //std::vector<rosic::rsString> audioOutputLongNames;
  //std::vector<rosic::rsString> audioOutputDescriptions;

  // it is actually a bad idea to store the input names in an array of strings because it implies
  // that each module of the same types has to store the same strings - what a waste!
  // Better: let Module have virtual functions getInputName(int index), getInputLongName(),..
  // and override these
  // maybe have a class ModuleTypeInfo - done - the new ModuleRegistry2 also serves as factory
  // ..when all code is changed to use this new class, get rid of the arrays above
  // ...hmm...but this wouldn't allow renaming of the pins at runtime - but such a feature is
  // only needed for very few modules ...hmmm


private:

  // copy and assignment - not possible:
  AtomicModule(const AtomicModule& /*other*/) {}
  AtomicModule& operator=(const AtomicModule& /*other*/) { return *this; }

};

//=================================================================================================
// class ParameterMixIn:

/** A mix-in class that can be used for modules that need parameters that should be presented on 
the GUI editor. We use a mix-in approach to make it possible to derive one specific module type 
with GUI-parameters from another specific module type that doesn't have GUI-parameters without 
resorting to multiple virtual inheritance. */

class ParameterMixIn
{

public:

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the parameter with given name to the given new value. If no parameter with the given 
  name is found, it returns false, otherwise true. */
  bool setParameter(const rosic::rsString &parameterName, const rosic::rsString &newValue, 
    bool callInternalCallback = true);

  /** Sets the parameter with given index to the given new value. When the "callInternalCallback" 
  parameter is true, it will call the virtual parameterChanged() callback function after setting 
  the value. You may override this function in order to do any internal variable updates that may 
  be necessary when a parameter changes. Normally, you want this but in certain circumstances it 
  may be necessary to supress the callback by passing false to avoid endless (indirect) 
  recursions. */
  void setParameter(int index, const rosic::rsString &newValue, bool callInternalCallback = true);
  // why the hell are parameters set via a string? that totally inefficient

  /** Appends a parameter with given name and default value to this module - the current value will
  also be initialized with the default value. */
  void addParameter(const rosic::rsString &parameterName, const rosic::rsString &defaultValue);
  // maybe we should distinguish between continuous, discrete/choice, boolean parameters and
  // have separate addContinuousParameter, .. etc. functions for them


  /** Internal callback that is triggered from setParameter - you may override it when you need to 
  re-compute some internal variables when a parameter was changed. */
  virtual void parameterChanged(int /*index*/) { }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of parameters that this module has. */
  int getNumParameters() const { return (int)parameters.size(); }

  /** Returns the name of the parameter with given index. */
  rosic::rsString getParameterName(int index) const;

  /** Returns the value of the parameter with given index. */
  rosic::rsString getParameterValue(int index) const;

  /** Returns the value of the parameter with the given name. If no parameter with the given name 
  is found, it return 0.0. */
  rosic::rsString getParameterValue(const rosic::rsString &parameterName) const;

  /** Returns the default value of the parameter with given index. */
  rosic::rsString getParameterDefaultValue(int index) const;

  /** Returns the index of the parameter with given name. Returns -1, if no matching name is 
  found. */
  int findIndexOfParameterWithName(const rosic::rsString &nameToFind) const;


protected:

  // data:
  struct Parameter
  {
    rosic::rsString name;
    rosic::rsString value;
    rosic::rsString defaultValue;
  };

  std::vector<Parameter> parameters;

};

//=================================================================================================
// class ModuleWithParameters:

/** Baseclass for modules that have some parameters that are adjusted from the GUI. This class is 
suitable as baseclass for modules with parameters that do not need to derive from other modules 
(without parameters). */

class ModuleWithParameters : public AtomicModule, public ParameterMixIn
{

protected:

  ModuleWithParameters(const std::string& name = std::string(), int x = 0, int y = 0, 
    bool polyphonic = false) : AtomicModule(name, x, y, polyphonic)
  {

  }

  virtual ~ModuleWithParameters()
  {

  }

private:

  // copy and assignment - not possible:
  ModuleWithParameters(const ModuleWithParameters& /*other*/) {}
  ModuleWithParameters& operator=(const ModuleWithParameters& /*other*/) { return *this; }

};

//=================================================================================================
// class ModuleProxy:

/** Baseclass for modules that must redirect/delegate the pointers (and related variables) for the
input signals to their connected source-modules and notify target modules when their source-module 
has changed. Input- and and output modules fall into this class. */

class ModuleProxy : public AtomicModule
{

public:

  virtual void mapApparentSourceToProcessingSource(Module * &sourceModule, int &sourceOutputPinIndex);
  virtual void mapProcessingSourceToSourceApparent(Module * &sourceModule, int &sourceOutputPinIndex);
  virtual void initialize();


  /** Assigns dummy functions that do nothing. */
  virtual void assignProcessingFunctions();

  static void processFrameDummy(Module* /*module*/, int /*voiceIndex*/)
  {
    //DEBUG_BREAK; // dummy function - should not be called
  }
  static void processBlockDummy(Module* /*module*/, int /*voiceIndex*/, int /*blockSize*/)
  {
    //DEBUG_BREAK; // dummy function - should not be called
  }


protected:


  ModuleProxy(const std::string& name = std::string(), int x = 0, int y = 0, 
    bool polyphonic = false);
  virtual ~ModuleProxy();


private:

  // copy and assignment - not possible:
  ModuleProxy(const ModuleProxy& /*other*/) {}
  ModuleProxy& operator=(const ModuleProxy& /*other*/) { return *this; }

};

} // end namespace romos

#endif 
