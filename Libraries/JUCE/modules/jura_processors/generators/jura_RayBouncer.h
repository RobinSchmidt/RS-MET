#ifndef jura_RayBouncer_h
#define jura_RayBouncer_h
  
//=================================================================================================

/**  */

class JUCE_API RayBouncerAudioModule : public jura::AudioModule
{

public:

  RayBouncerAudioModule(CriticalSection *lockToUse);
    
  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createParameters();

  /** Creates the GUI editor (returns an object of an appropriate subclass of AudioModuleEditor) */
  //AudioModuleEditor *createEditor() override;

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override;

  // target functions for callbacks that are called on parameter changes:
  void setFrequency(double RayBouncerAudioModule);

protected:

  // embedded core DSP objects from the RAPT library:
  //rsRayBouncerDriverD rayBouncer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Ladder)
};

//=================================================================================================

/** This is the GUI editor class for the jura::RayBouncerAudioModule. */
/*
class JUCE_API RayBouncerAudioModuleEditor : public AudioModuleEditor
{

public:

  RayBouncerAudioModuleEditor(jura::RayBouncerAudioModule *newRayBouncerToEdit);
  virtual void resized() override;

protected:

  RayBouncerAudioModule *rayBouncerToEdit;

  AutomatableSlider *frequencySlider;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RayBouncerAudioModuleEditor)
};
*/

#endif 