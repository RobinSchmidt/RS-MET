#ifndef jura_EllipseOscillator_h
#define jura_EllipseOscillator_h


class JUCE_API EllipseOscillatorAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  EllipseOscillatorAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);
    // maybe make a constructor without the managers

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

  rosic::rsEllipseOscillator oscCore;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(EllipseOscillatorAudioModule)
};

#endif