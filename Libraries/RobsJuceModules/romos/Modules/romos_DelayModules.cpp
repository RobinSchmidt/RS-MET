//#include "romos_DelayModules.h"

  void UnitDelayModule::initialize()
  { 
    initInputPins({ "" });
    initOutputPins({ "" });
    hasHeaderFlag = false;
  }
  INLINE void UnitDelayModule::process(Module *module, double *in, double *out, int voiceIndex)
  {
    UnitDelayModule *unitDelay = static_cast<UnitDelayModule*> (module);
    double *buf = unitDelay->buffer + voiceIndex;  // actually  "+ voiceIndex * numPins" but numPins is unity
    *out = *buf;  
    *buf = *in;
  }
  void UnitDelayModule::resetVoiceState(int voiceIndex)
  {
    AtomicModule::resetVoiceState(voiceIndex);
    buffer[voiceIndex] = 0.0;
  } 
  void UnitDelayModule::allocateMemory()
  {
    AtomicModule::allocateMemory();
    buffer = new double[getNumVoices()]; // x[n-1]
  }
  void UnitDelayModule::freeMemory()
  {
    AtomicModule::freeMemory();
    delete[] buffer;
    buffer = NULL;
  }
  CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(UnitDelayModule); 

