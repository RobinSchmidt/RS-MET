#ifndef jura_FourOscSection_h
#define jura_FourOscSection_h

/** This class wraps rosic::FourOscSection into a rosof::AudioModule to facilitate its use as 
plugIn. */

class FourOscSectionAudioModule : public AudioModule
{

  friend class FourOscSectionModuleEditor;

  friend class StraightlinerAudioModule;
    // this is a backward compatibility artifact: required for state recall from older versions 
    // states of Straightliner

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FourOscSectionAudioModule(CriticalSection *newPlugInLock, 
    rosic::FourOscSection *fourOscSectionToWrap);

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedFourOscSection->setSampleRate(newSampleRate);
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    //if( wrappedFourOscSection != NULL )
    //  wrappedFourOscSection->getSampleFrameStereo(inOutL, inOutR); 
  }

protected:

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::FourOscSection *wrappedFourOscSection;

  // we maintain wrappped versions (into rosof::AudioModules) of the synth's bzuilding blocks 
  // here in order to make them automatable:
  OscillatorStereoAudioModule *osc1Module, *osc2Module, *osc3Module, *osc4Module;


  juce_UseDebuggingNewOperator;
};

//=============================================================================================



#endif
