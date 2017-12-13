#ifndef jura_RotationOscillator_h
#define jura_RotationOscillator_h


class JUCE_API RotationOscillatorAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  RotationOscillatorAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);

  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  virtual void noteOn(int noteNumber, int velocity) override;

protected:

  RAPT::rsLissajousOscillator3D<double> oscCore;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RotationOscillatorAudioModule)
};

#endif