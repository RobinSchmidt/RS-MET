#ifndef rosic_ModulatorModules_h
#define rosic_ModulatorModules_h

// rosic-indcludes:
#include "rosic_Module.h"
#include "../modulators/rosic_BreakpointModulator.h"

namespace rosic
{

  /**

  This file defines wrapper classes that wrap some core modulation-source objects into
  ModulationSource objects to facilitate their use in a (semi) modular framework.

  */

  class BreakpointModulatorModule : public ModulationSource, public BreakpointModulator
  {
  public:
    virtual void setSampleRate(float newSampleRate)
    { BreakpointModulator::setSampleRate(newSampleRate); }
    virtual void setTempoInBPM(float newTempo)
    { BreakpointModulator::setBeatsPerMinute(newTempo); }
    virtual void trigger() { BreakpointModulator::noteOn(); } // todo: include key/velocity
    virtual double getSample() { return BreakpointModulator::getSample(); }
  };

} // end namespace rosic

#endif
