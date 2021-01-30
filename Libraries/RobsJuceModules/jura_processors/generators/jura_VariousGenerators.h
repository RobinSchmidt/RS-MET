#ifndef jura_EllipseOscillator_h
#define jura_EllipseOscillator_h
// rename to Oscillators
//=================================================================================================

/** A thin wrapper around rosic::SineOscillator to facilitate use as an AudioModule. Maybe get 
rid of that. */

class JUCE_API SineOscCore : private rosic::SineOscillator
{

public:

  typedef rosic::SineOscillator Base;
  SineOscCore() {}
  void setSampleRate(double newRate)   { Base::setSampleRate(newRate); }
  void setAmplitude( double newAmp)    { amp = newAmp; }
  void setDetune(    double newDetune) { detune = newDetune; updateFreq(); }
  void setNoteKey(   double newKey)    { key    = newKey;    updateFreq(); }

  double getSample() { return amp * Base::getSample(); }
  void   reset()     { Base::trigger(); }

protected:

  void updateFreq() { Base::setFrequency(RAPT::rsPitchToFreq(key+detune)); }

  double key = 64;
  double detune = 0.0, amp = 1.0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SineOscCore);
  //JUCE_LEAK_DETECTOR(SineOscCore);
};

//=================================================================================================

/** A simple sine oscillator. Serves mostly as a simple/minimal example oscillator during 
development of the framework - may not be included into any product. */

class JUCE_API SineOscAudioModule : public jura::AudioModuleWithMidiIn
{

public:

  SineOscAudioModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  // overriden from AudioModule baseclass:
  virtual void processBlock(double** inOutBuffer, int numChannels, int numSamples) override
  {
    jassert(numChannels == 2);
    for(int n = 0; n < numSamples; n++)
      inOutBuffer[0][n] = inOutBuffer[1][n] = core.getSample();
  }
  virtual void processStereoFrame(double* L, double* R) override { *L = *R = core.getSample(); }

  virtual void setSampleRate(double newRate) override { core.setSampleRate(newRate); }
  virtual void reset()                       override { core.reset(); }
  virtual void noteOn(int key, int vel)      override { core.setNoteKey(key); }

protected:

  virtual void createParameters();

  SineOscCore core;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SineOscAudioModule)
};

class JUCE_API SineOscAudioModulePoly : public jura::AudioModulePoly
{

public:

  SineOscAudioModulePoly(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  virtual void noteOn(int key, int vel, int voice) override
  {
    jassert(voice >= 0 && voice < voiceManager->getMaxNumVoices());
    // see comment in AttackDecayEnvelopeModulePoly

    // We want to retrigger the phase of the oscillator but only if the voice was allocated freshly
    // and not revived from relase

    voices[voice].reset();
    // hmm...but i think this retriggering should only be done, if the voice was previously silent
    // ...in cases where a releasing voice is retriggered ()

  }

  void processStereoFrameVoice(double* left, double* right, int voice) override
  { 
    *left = *right = voices[voice].getSample(); 
  }


  void setSampleRate(double newSampleRate) override { sampleRate = newSampleRate; }


  // parameter callback targets:
  void setFrequency(double newFrequency, int voice);
  void setAmplitude(double newAmplitude, int voice);
  void setDetune(   double newDetune,    int voice);
  // maybe just have frequency and amplitude parameters and handle the adjustment of frequency via
  // the mod-system -> provide a polyphonic Moudlator for NotePitch (but hwat range?) and/or 
  // NoteFrequency that is 
  // always available in ToolChain - maybe also NoteVelocity (0..1), PitchBend (-1..+1)


protected:

  virtual void createParameters();
  void allocateVoiceResources() override;

  SineOscCore core;  // maybe rename to master
  //std::vector<RAPT::rsSineOscillator<double>> voices;
  std::vector<RAPT::rsSineOscillatorNaive<double>> voices;

  double sampleRate;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SineOscAudioModulePoly)
};


//class JUCE_API AttackDecayEnvelopeModulePoly : public AudioModulePoly

//=================================================================================================

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

//=================================================================================================

/** Oscillator based on morph between saw-up/triangle/saw-down waveform. The raw waveform is then
further shaped to produce sinusoidish and squareish waveforms. */

class JUCE_API TriSawOscModule : public jura::AudioModuleWithMidiIn
{

public:

  TriSawOscModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, ModulationManager* modManagerToUse = nullptr);
  // maybe make a constructor without the managers

  /** Creates the static parameters for this module (i.e. parameters that are not created
  dynamically and are thus always there). */
  virtual void createParameters();

  // overriden from AudioModule baseclass:
  //virtual void processBlock(double **inOutBuffer, int numChannels, int numSamples) override;
  virtual void processStereoFrame(double *left, double *right) override;
  virtual void setSampleRate(double newSampleRate) override;
  virtual void reset() override;
  virtual void noteOn(int noteNumber, int velocity) override;

  // parameter callback targets:
  void setAmplitude(double newAmplitude) { amplitude = newAmplitude; }
  void setBend(double newBend);
  void setBendAsym(double newAsym);
  //void setSigmoid(double newSigmoid);  // or maybe call it sinusoidality/smoothness
  //void setSigmoidAsym(double newAsym);

  // asym: triangle vs saw, or up vs down, pulse-width, hi vs low, "Upness"
  // bend: squareness vs thinness, thick vs thin, "Thickness"
  // sigmoid: sinusoidality/smoothness, "Smoothness"

protected:

  void updateBending();

  RAPT::rsTriSawOscillator<double> oscCore;

  // parameters:
  double freq = 0, sampleRate = 44100;
  double bend = 0, bendAsym = 0;
  double sigmoid = 0, sigmoidAsym = 0;
  double amplitude = 1;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TriSawOscModule)
};

//=================================================================================================





#endif