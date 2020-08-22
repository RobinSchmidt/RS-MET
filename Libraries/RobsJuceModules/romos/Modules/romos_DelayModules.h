#ifndef romos_DelayModules_h
#define romos_DelayModules_h

//-------------------------------------------------------------------------------------------------

class UnitDelayModule : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(UnitDelayModule);
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  double *buffer;
};
class UnitDelayTypeInfo : public ModuleTypeInfo
{
public:
  UnitDelayTypeInfo() {
    shortName    = "D";
    fullName     = "UnitDelay";
    description  = "Delays the input signal by one sample";
    category     = "Delays";
    createModule =  []()->Module* { return new UnitDelayModule; };
    hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------

// todo: MultiUnitDelay/UnitDelayChain: provides outputs for x[n-1], x[n-2], ... - useful for 
// building filters

// DelayLine: user can choose interpolation (round, linear, allpass, etc.) on gui
// Inputs: Audio-In, Delay (in seconds)
// Outputs: delayed audio

// MultiTapDelay: x[n-N1], x[n-N2], ...
// Inputs: Audio-In, Delay-times
// Outputs: delayed audio signals

#endif 
