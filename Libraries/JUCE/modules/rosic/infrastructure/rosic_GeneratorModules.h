#ifndef rosic_GeneratorModules_h
#define rosic_GeneratorModules_h

//// rosic-indcludes:
//#include "rosic_Module.h"
//#include "../generators/rosic_OscillatorStereo.h"

namespace rosic
{

  /**

  This file defines wrapper classes that wrap some core signal-processing objects into 
  Module objects to facilitate their use in a (semi) modular framework such as Quadrifex.

  */

  class OscillatorStereoModule : public Module, public OscillatorStereo
  {
  public:
    virtual void setSampleRate(double newSampleRate) 
    { OscillatorStereo::setSampleRate(newSampleRate); }
    virtual void reset() { OscillatorStereo::reset(); }
    OscillatorStereoModule()
    {
      modulatableParameters.appendElement( ModulatableParameter("Amplitude", 1.0) );
      modulatableParameters.appendElement( ModulatableParameter("Detune",    0.0) );
    }
    virtual void processSampleFrame(double *inOutL, double *inOutR) 
    {
      OscillatorStereo::setModulatedAmplitude( modulatableParameters[0].getInstantaneousValue() );
      // set Detune, ...
      OscillatorStereo::getSampleFrameStereo(inOutL, inOutR); 
    }
  };




} // end namespace rosic

#endif 
