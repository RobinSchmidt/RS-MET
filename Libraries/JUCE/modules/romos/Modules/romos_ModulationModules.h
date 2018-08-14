#ifndef romos_ModulationModules_h
#define romos_ModulationModules_h

#include "../Framework/romos_ModuleAtomic.h"
#include "romos_ModuleDefinitionMacros.h"

namespace romos
{




  /** Generates a classic Attack-Decay-Sustain-Release (ADSR) envelope with variable shapes for each phase.
  Parameters:
   none
  Inputs:
   0: Att -> attack time (in seconds)
   1: AtSh -> attack shape (in number of time constants between the points, sing decides between growth/decay)
   2: Dec -> decay time (in seconds)
   3: DcSh -> decay shape
   4: Sus -> sustain level (as value between 0 and 1)
   5: Rel -> relase time in seconds
   6: RlSh -> release shape

   // \todo input for:
   7: Trg -> retrigger envelope if nonzero
   8: TmScl -> time scaling factor for the whole envelope (as raw factor) - nah -> get rid of that, redundant with A,D,R

  Outputs:
   0: the envelope signal

  References: Moore - Elements of Computer Music, page 183...
  todo: change the meaning of the Shape parameters to be the normalized y-value at x=0.5 shifted
  and scaled to the range -1..+1  */
  class EnvelopeADSR : public ModuleAtomic
  {
    CREATE_COMMON_DECLARATIONS_7(EnvelopeADSR);
  public:
    virtual void resetVoiceState(int voiceIndex);
  protected:


    static void processWithoutTriggerFlagCheck(Module *module, const double *in1, const double *in2, const double *in3, const double *in4, 
      const double *in5, const double *in6, const double *in7, double *out, const int voiceIndex);
      // to be invoke from the actual process function


    virtual void allocateMemory();
    virtual void freeMemory();

    unsigned long *counters;
    double        *startValues;
  };


  






} 

#endif 
