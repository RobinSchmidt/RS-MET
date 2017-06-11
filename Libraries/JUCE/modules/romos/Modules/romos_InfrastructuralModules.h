#ifndef romos_InfrastructuralModules_h
#define romos_InfrastructuralModules_h

#include "romos_ArithmeticModules.h"

namespace romos
{

//-------------------------------------------------------------------------------------------------
// infrastructural modules:

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
  friend class ModuleContainer;
};

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
  friend class ModuleContainer;
};

/** Outputs the current system sample rate. */
class SystemSampleRateModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(SystemSampleRateModule);
};

/** Outputs the current system sampling period, that is: the distance between two samples in 
seconds. This is equal to the reciprocal of the samplerate. */
class SystemSamplePeriodModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(SystemSamplePeriodModule);
};

/** Outputs 1 if the voice is currently playing a note, 0 otherwise. */
class NoteGateModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(NoteGateModule);
};

/** Outputs 1 if the voice has just been triggered (note-on), 0 otherwise. */
class NoteOnTriggerModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(NoteOnTriggerModule);
};

/** Outputs 1 if the voice has just been released (note-off), 0 otherwise. */
class NoteOffTriggerModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(NoteOffTriggerModule);
};

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


/** Adds polyphonic signals at its input into a monophonic signal at its output. */
class VoiceCombinerModule : public ModuleAtomic
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


/** Outputs the frequency of the most recent note that is/was played in this voice. */
class NoteFrequencyModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(NoteFrequencyModule);
};

/** Outputs the velocity the note that is currently played in this voice. */
class NoteVelocityModule : public ModuleAtomic
{
  CREATE_COMMON_DECLARATIONS_0(NoteVelocityModule);
};


// \todo AfterTouchModules (or ChannelPressure and NotePressure (or KeyPressure) modules),


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
  virtual void setModuleName(const rosic::rsString& newName);

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

}

#endif 
