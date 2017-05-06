#ifndef jura_SimpleSampler_h
#define jura_SimpleSampler_h

/** This class wraps rosic::SimpleSampler into a rosof::AudioModule to facilitate its use as 
plugIn.  */

class SimpleSamplerAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class SimpleSamplerModuleEditor;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  SimpleSamplerAudioModule(CriticalSection *newPlugInLock, 
    rosic::SimpleSampler *simpleSamplerToWrap);

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Do we really need to override this ?! */
  virtual void setSampleRate(double newSampleRate)
  {
    wrappedSimpleSampler->setSampleRate(newSampleRate);
  }

  virtual void reset()
  {
    wrappedSimpleSampler->resetAllVoices();
  }

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedSimpleSampler->getSampleFrameStereo(inOutL, inOutR);
  }


protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the sampler's building blocks 
  // here in order to make them automatable:
  SamplePlayerAudioModule        *samplePlayerModule;
  MultiModeFilterAudioModule     *filterModule;
  BreakpointModulatorAudioModule *pitchEnvModule, *filterEnvModule, *ampEnvModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::SimpleSampler *wrappedSimpleSampler;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif 
