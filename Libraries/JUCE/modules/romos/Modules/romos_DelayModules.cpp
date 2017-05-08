#include "romos_DelayModules.h"

namespace romos
{


  void UnitDelayModule::initialize()
  { 
    initInputPins( 1, rosic::String());
    initOutputPins(1, rosic::String());
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
    ModuleAtomic::resetVoiceState(voiceIndex);
    buffer[voiceIndex] = 0.0;
  } 
  void UnitDelayModule::allocateMemory()
  {
    ModuleAtomic::allocateMemory();
    buffer = new double[getNumVoices()]; // x[n-1]
  }
  void UnitDelayModule::freeMemory()
  {
    ModuleAtomic::freeMemory();
    delete[] buffer;
    buffer = NULL;
  }
  CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(UnitDelayModule); 

}
