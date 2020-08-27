#ifndef romos_InfrastructuralModules_h
#define romos_InfrastructuralModules_h

//-------------------------------------------------------------------------------------------------

/** We use identity modules to realize input and output - but we make some special subclasses in 
order to override the default module names. maybe factor out a baseclass PinModule */

class AudioInputModule : public IdentityModule
{
  ENFORCE_FACTORY_USAGE(AudioInputModule);
public:
  //virtual void initialize();
  //virtual void resetState() { } // overriden to avoid heap-corruption - baseclass function would write into invalid memory areas
  virtual unsigned int getNumInputPins()  const { return 0; }
  virtual unsigned int getNumOutputPins() const { return 1; }
protected:
  //virtual void allocateMemory();
  //virtual void freeMemory();
  //virtual void setAudioInputAddress(double *newAddress);
  //virtual void setAudioOutputAddress(double *newAddress);
  friend class ContainerModule;
};
class AudioInputTypeInfo : public ModuleTypeInfo
{
public:
  AudioInputTypeInfo() {
    shortName    = "In";
    fullName     = "AudioInput";
    description  = "Input module for feeding signals into containers";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new AudioInputModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class AudioOutputModule : public IdentityModule
//class AudioOutputModule : public ModuleProxy
{
  ENFORCE_FACTORY_USAGE(AudioOutputModule);
public:
  //virtual void initialize();
  //virtual void resetState() { }
  virtual unsigned int getNumInputPins()  const { return 1; }
  virtual unsigned int getNumOutputPins() const { return 0; }
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  //virtual void resetState();
  friend class ContainerModule;
};
class AudioOutputTypeInfo : public ModuleTypeInfo
{
public:
  AudioOutputTypeInfo() {
    shortName    = "Out";
    fullName     = "AudioOutput";
    description  = "Output module for feeding signals out of containers";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new AudioOutputModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

class SystemSampleRateModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(SystemSampleRateModule);
};
class SystemSampleRateTypeInfo : public ModuleTypeInfo
{
public:
  SystemSampleRateTypeInfo() {
    shortName    = "SR";
    fullName     = "SampleRate";
    description  = "Outputs the system's sample rate in Hz";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new SystemSampleRateModule; };
    hasHeader    = false;
    hasEditor    = false;
  }
};

//-------------------------------------------------------------------------------------------------

/** Outputs the current system sampling period, that is: the distance between two samples in 
seconds. This is equal to the reciprocal of the samplerate. */

class SystemSamplePeriodModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(SystemSamplePeriodModule);
};
class SystemSamplePeriodTypeInfo : public ModuleTypeInfo
{
public:
  SystemSamplePeriodTypeInfo() {
    shortName    = "1/SR";
    fullName     = "SampleInterval";
    description  = "Outputs the system's sampling interval in seconds (= 1/SampleRate)";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new SystemSamplePeriodModule; };
    hasHeader = false;
  }
};
// rename to SampleInterval

//-------------------------------------------------------------------------------------------------

class NoteGateModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(NoteGateModule);
};
class NoteGateTypeInfo : public ModuleTypeInfo
{
public:
  NoteGateTypeInfo() {
    shortName    = "NoteGate";
    fullName     = "NoteGate";
    description  = "Outputs 1, if the voice is currently playing a note, 0 otherwise";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new NoteGateModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class NoteOnTriggerModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(NoteOnTriggerModule);
};
class NoteOnTriggerTypeInfo : public ModuleTypeInfo
{
public:
  NoteOnTriggerTypeInfo() {
    shortName    = "NoteOn";
    fullName     = "NoteOnTrigger";
    description  = "Outputs 1, if there was a note-on at this sample for this voice, 0 otherwise";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new NoteOnTriggerModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class NoteOffTriggerModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(NoteOffTriggerModule);
};
class NoteOffTriggerTypeInfo : public ModuleTypeInfo
{
public:
  NoteOffTriggerTypeInfo() {
    shortName    = "NoteOff";
    fullName     = "NoteOffTrigger";
    description  = "Outputs 1 if there was a note-off at this sample for this voice, 0 otherwise";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new NoteOffTriggerModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

/** Kills the voice when the maximum absolute value of the input signal remains below a threshold 
for some specified time. */

class VoiceKillerModule : public ModuleWithParameters
{
  CREATE_COMMON_DECLARATIONS_1(VoiceKillerModule);
public:
  virtual void resetVoiceState(int voiceIndex);
  virtual void parameterChanged(int index);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  unsigned int *sampleCounters;

  double threshold, timeOut;
};
class VoiceKillerTypeInfo : public ModuleTypeInfo
{
public:
  VoiceKillerTypeInfo() {
    shortName    = "VoiceKill";
    fullName     = "VoiceKiller";
    //description  = "Kills the voice when the maximum absolute value of the input signal remains below a threshold for some specified time";
    description  = "Kills voice when max-abs of input remains below threshold for some time";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new VoiceKillerModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class VoiceCombinerModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(VoiceCombinerModule);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void connectInputPinTo(int inputPinIndex, Module *sourceModule, 
    int sourceOutputPinIndex);
  virtual void updateInputPointersAndInFrameStrides();
  virtual void setPolyphonic(bool shouldBePolyphonic);
  virtual void clearVoiceBuffer(int voiceIndex);
};
class VoiceCombinerTypeInfo : public ModuleTypeInfo
{
public:
  VoiceCombinerTypeInfo() {
    shortName    = "}";
    fullName     = "VoiceCombiner";
    description  = "Adds polyphonic signals at its input into a monophonic signal at its output";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new VoiceCombinerModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class NoteFrequencyModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(NoteFrequencyModule);
};
class NoteFrequencyTypeInfo : public ModuleTypeInfo
{
public:
  NoteFrequencyTypeInfo() {
    shortName    = "NoteFreq";
    fullName     = "NoteFrequency";
    description  = "Frequency of the most recent note that is/was played in this voice";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new NoteFrequencyModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

class NoteVelocityModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_0(NoteVelocityModule);
};
class NoteVelocityTypeInfo : public ModuleTypeInfo
{
public:
  NoteVelocityTypeInfo() {
    shortName    = "NoteVelo";
    fullName     = "NoteVelocity";
    description  = "Velocity of the most recent note that is/was played in this voice";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new NoteVelocityModule; };
    hasHeader = false;
  }
};

// \todo AfterTouchModules (or ChannelPressure and NotePressure (or KeyPressure) modules),

//-------------------------------------------------------------------------------------------------

/** Outputs a number that represents some user parameter. It facilitates presentation to the user 
by defining a range and mapping. */

class ParameterModule : public ConstantModule, public ParameterMixIn
{
  ENFORCE_FACTORY_USAGE(ParameterModule);

public:

  enum mappingFunctions
  {
    LINEAR_MAPPING = 0,
    EXPONENTIAL_MAPPING
  };

  virtual void initialize();

  virtual void resetVoiceState(int voiceIndex); 
  // later: set current-value to target value there (smoothing)

  virtual void parameterChanged(int index);

  // re-overriden to fall back to Module::setModuleName(rosic::rsString(value)); - we don't want 
  // to set the value here
  virtual void setModuleName(const std::string& newName);

  /** Sets the minimum- and maximum-values and the mapping function simultaneously. We have this 
  function because min, max and mapping  must satisfy certain contraints (for example min <= max, 
  min > 0 for exponential mapping, etc.) and the ParameterModule class may refuse to take values 
  when setting one at a time that violates the constraint. For example, when the currrent maximum 
  is 1 and you try to set a new minimum > 1, it will actually set the new minimum to 1, etc. */
  virtual void setMinMaxAndMapping(double newMin, double newMax, int newMapping);

  /*
  virtual void setValue(       double newValue);
  virtual void setMinValue(    double newMinValue);
  virtual void setMaxValue(    double newMaxValue);
  virtual void setDefaultValue(double newDefaultValue);
  virtual void setQuantization(double newQuantization);
  // setMappingFunction
  */

  /** sets up the current value of the parameter according to the passed controller-value. */
  virtual void setValueFromController(double controllerValue);

  /** Given the 4 snapshot-indices s1...s4 at 4 corners of a square and an x- and y-coordinate 
  (between 0...1), this function sets up the value of the parameter to an interpolated value 
  between the 4 corner snapshot values. The interpolation is bilinear in the normalized parameter 
  domain. */
  virtual void setValueFromSnapshots(int topLeftIndex, int topRightIndex, int bottomLeftIndex, 
    int bottomRightIndex, double x, double y);

  virtual double getValue()           const { return value; }
  virtual double getMinValue()        const { return minValue; }
  virtual double getMaxValue()        const { return maxValue; }
  virtual double getDefaultValue()    const { return defaultValue; }
  virtual double getQuantization()    const { return quantization; }
  virtual int    getMappingFunction() const { return mappingFunction; }

protected:

  /** Enforces the values of minValue, maxValue, etc. to satisfy certain constrains (for exmple 
  minValue <= maxValue, etc.). */
  virtual void enforceConsistencyOfValues();

  /** Maps a normalized value in the range 0...1 (both inclusive) to the actual parameter range 
  minValue...maxValue (both inclusive) according to the selected mappingFunction. */
  double mapNormalizedValue(double normalizedValue);

  /** Given an actual parameter value in the range minValue...maxValue (both inclusive), this 
  function returns the corresponding normalized value in the range 0...1 (both inclusive). */
  double unmapToNormalizedValue(double mappedValue);

  double minValue;
  double maxValue;
  double defaultValue;
  double quantization;

  // maybe have a normalizedValue member - it would make snapshot-interpolation morde efficient 
  // (because this is done in the normalized domain), but it also introduce redundancy in the 
  // data ...we'll see

  int mappingFunction;

  int assignedController;

  std::vector<double> snapshotValues;

};
class ParameterModuleTypeInfo : public ModuleTypeInfo
{
public:
  ParameterModuleTypeInfo() {
    shortName    = "Param";
    fullName     = "Parameter";
    description  = "User parameter that is set up on the GUI";
    category     = "Infrastructure";
    createModule =  []()->Module* { return new ParameterModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------


/** Outputs a number that represents some user parameter. It facilitates presentation to the
user by defining a range and mapping. */
/*
class ParameterModule : public NumberModule
{
// todo: incorporate range and scaling/mapping functions
};
*/

// \todo: NoteOnTrigger, NoteOffTrigger, NoteOnTriggerIfVoiceWasSilentBefore, NoteKey, 
// VoiceKiller -> Threshold, TimeOut, In


#endif 
