#ifndef rosic_asSineOscillatorStereo_h
#define rosic_asSineOscillatorStereo_h

// rosic-indcludes:
#include "../../generators/rosic_SineOscillatorStereo.h"

#include "../../_third_party/angelscript/angelscript.h"

namespace rosic
{

  /**

  This class 

  */

  class asSineOscillatorStereo  
  {
  public:
    asSineOscillatorStereo();  
    void addRef();
    void release();
    static void registerObjectType(asIScriptEngine *engine);

    void   setSampleRate(double newSampleRate);
    void   setFrequency(double newFrequency);
    void   setAmplitude(double newAmplitude);
    void   trigger();
    void   getSampleFrameStereo(double *outL, double *outR);
    double getSample();

    SineOscillatorStereo sineOscillatorStereo;

  protected:
    ~asSineOscillatorStereo();  
    int refCount;
  };

} // end namespace rosic

#endif // #ifndef rosic_asSineOscillatorStereo_h
