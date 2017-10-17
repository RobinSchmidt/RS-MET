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

  virtual void setSampleRate(double newSampleRate) override
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

  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override
  {
    jassertfalse; // no code yet dragged over - maybe we should accumulate the outputs of
    // our embedded OscillatorStereoAudioModule objects - i.e. clear the buffer and let
    // OscillatorStereoAudioModule have a function:
    // accumulateIntoBlock(double **inOutBuffer, int numChannels, int numSamples)
    // and call that here...hmm...but how is that supposed to work when rosic::Straightliner
    // has a rosic::FourOscSection object - or maybe it hasn't? check that out....
  }

protected:

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::FourOscSection *wrappedFourOscSection;

  // we maintain wrappped versions (into rosof::AudioModules) of the synth's building blocks 
  // here in order to make them automatable:
  OscillatorStereoAudioModule *osc1Module, *osc2Module, *osc3Module, *osc4Module;


  juce_UseDebuggingNewOperator;
};

//=============================================================================================

/** An editor for a section of 4 oscillators. */

class FourOscSectionModuleEditor : public AudioModuleEditor 
{

public:
  
  /** Constructor. You must pass 4 valid pointers to OscillatorStereoAudioModule objects. */
  FourOscSectionModuleEditor(CriticalSection *newPlugInLock, 
    FourOscSectionAudioModule *newFourOscSectionToEdit);

  // callbacks:
  virtual void resized();
  virtual void updateWidgetsAccordingToState();

protected:

  // the 4 sub-editors for the 4 oscillators:
  OscillatorStereoEditor *osc1Editor, *osc2Editor, *osc3Editor, *osc4Editor;

  juce_UseDebuggingNewOperator;
};

#endif
