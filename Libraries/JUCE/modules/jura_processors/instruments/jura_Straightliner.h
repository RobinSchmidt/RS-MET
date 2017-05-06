#ifndef jura_Straightliner_h
#define jura_Straightliner_h


/** This class wraps rosic::Straightliner into a rosof::AudioModule to facilitate its use as 
plugIn. */

class StraightlinerAudioModule : public PolyphonicInstrumentAudioModule
{

  friend class StraightlinerModuleEditor;

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StraightlinerAudioModule(CriticalSection *newPlugInLock, rosic::Straightliner *straightlinerToWrap);

  //---------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->setSampleRate(newSampleRate);
  }

  /** Checks, if this is a cracked version and if so, it sets up the appendix for the headline
  accordingly. Return value informs also whether or not a cracked version was detected. */
  virtual bool checkForCrack();

  /** Loads the user's custom preferences such as the sample-content path. */
  //virtual void loadPreferences();

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedStraightliner->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    if(wrappedStraightliner->isSilent())
    {
      rosic::fillWithZeros(left, numSamples);
      rosic::fillWithZeros(right, numSamples);
    }
    else
    {
      for(int n=0; n<numSamples; n++)
      {
        double dL = (double)left[n];
        double dR = (double)right[n];
        wrappedStraightliner->getSampleFrameStereo(&dL, &dR);
        left[n]  = (float)dL;
        right[n] = (float)dR;
      }
    }
  }

  //---------------------------------------------------------------------------------------------
  // others:

  virtual void reset()
  {
    if(wrappedStraightliner != NULL)
      wrappedStraightliner->resetAllVoices();
  }

  /** We must override this here because we have - for historical reasons - not the FourOscScion
  as child-module but each osc separately, so we must manually tak care of updating the preset
  field of the FourOscSectionEditor here. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean);

protected:

  // we maintain wrappped versions (into rosof::AudioModules) of the synth's building blocks 
  // here in order to make them automatable:
  FourOscSectionAudioModule      *oscSectionModule;
  MultiModeFilterAudioModule     *filterModule;
  BreakpointModulatorAudioModule *pitchEnvModule, *filterEnvModule, *ampEnvModule;

  /** Pointer to the underlying rosic object which is wrapped. */
  rosic::Straightliner *wrappedStraightliner;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


#endif
