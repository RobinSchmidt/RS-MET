#ifndef romos_DelayModules_h
#define romos_DelayModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{

  /** Delays the input signal by one sample. */
  class UnitDelayModule : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_1(UnitDelayModule);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:
    virtual void allocateMemory();
    virtual void freeMemory();
    double *buffer;    
  };

} 

#endif 
