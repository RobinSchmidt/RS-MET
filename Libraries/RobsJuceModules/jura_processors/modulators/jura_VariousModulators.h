#ifndef jura_VariousModulators_h
#define jura_VariousModulators_h

class JUCE_API TriSawModulatorModule : public AudioModuleWithMidiIn, public ModulationSource
{

public:

  TriSawModulatorModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);


  virtual void setSampleRate(double newSampleRate) override { core.setSampleRate(newSampleRate); }
  virtual void reset() override { core.reset();  }
  virtual void noteOn(int noteNumber, int velocity) override { core.reset(); }
  virtual double getModulatorOutputSample() override { return core.getSample(); }

  //virtual void updateModulationValue() override { modValue = core.getSample();  }
  // todo: change interface in order to return the modValue and let the framework take care of 
  // where to store it

protected:

  virtual void createParameters();


  rosic::rsTriSawModulator core; // maybe use a pointer to be able to wrap existing objects

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TriSawModulatorModule)
};

//=================================================================================================

class JUCE_API AttackDecayEnvelopeModule 
  //: public AudioModulePoly, public ModulationSourcePoly  // for later use
  : public AudioModuleWithMidiIn, public ModulationSource 
{

public:

  AttackDecayEnvelopeModule(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);


  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override { core.reset();  }

  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOff(key, 0);  }

  virtual double getModulatorOutputSample() override { return core.getSample(); }

  // parameter callback targets:
  void setAttack(double newAttack);
  void setDecay( double newDecay);
  // maybe for the polyphonic case, these functions should just receive a second argument for the 
  // voice index? this will be a quite flexible design but will require to write a lot of 
  // boilerplate code to make objects polyphonic - for each modulatable parameter, we must write a 
  // boilerplate per-voice callback function

protected:

  virtual void createParameters();

  RAPT::rsAttackDecayEnvelope<double> core;

  // parameters:
  double sampleRate = 44100;
  double attack = 10.0, decay = 100.0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AttackDecayEnvelopeModule)
};

// maybe it's best to have a monophonic and a polyphonic version of this class - just to show the 
// differences - for future modules, we'll only do the polyphonic one

//=================================================================================================

class JUCE_API AttackDecayEnvelopeModulePoly : public AudioModulePoly
{

public:

  AttackDecayEnvelopeModulePoly(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr,
    rosic::rsVoiceManager* voiceManagerToUse = nullptr);

  virtual ~AttackDecayEnvelopeModulePoly() {}

  /*
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override { core.reset();  }

  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOff(key, 0);  }

  virtual double getModulatorOutputSample() override { return core.getSample(); }
  */

  // parameter callback targets:
  void setAttack(double newAttack, int voice);
  void setDecay( double newDecay,  int voice);
  // maybe for the polyphonic case, these functions should just receive a second argument for the 
  // voice index? this will be a quite flexible design but will require to write a lot of 
  // boilerplate code to make objects polyphonic - for each modulatable parameter, we must write a 
  // boilerplate per-voice callback function

protected:

  // called once in the constructor:
  //virtual void createCores();
  virtual void createParameters();


  void allocateVoiceResources() override;


  //std::vector<RAPT::rsAttackDecayEnvelope<double>*> cores;
  // no - that's not good - we should have one full rsAttackDecayEnvelope object with the full set
  // of parameters and the voices should be leaner objects that only store the voice-dependent 
  // state

  // maybe have an array of direct objects, not pointers - easier to handle:
  std::vector<RAPT::rsAttackDecayEnvelope<double>> cores;


  // parameters:
  double sampleRate = 44100;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AttackDecayEnvelopeModulePoly)
};






#endif