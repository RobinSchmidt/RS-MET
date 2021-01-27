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
  virtual double renderModulation() override { return core.getSample(); }

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

  virtual double renderModulation() override { return core.getSample(); }

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

// some very simple modulators

/** A modulator module that just outputs the constant value 1. That may be useful for use with 
exponentially scaled parameters such as a frequency that goes from 20 to 20000. We can just 
subtract 20 by connecting the constant with depth -20 and the use additive envelopes to directly
determine the freq without having the minimum getting in the way. */
class JUCE_API rsConstantOneModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsConstantOneModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    setModulationSourceName("ConstantOne");
  }

  double renderModulation()                    override { return 1.0; } 
  double renderVoiceModulation(int voiceIndex) override { return 1.0; }

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsConstantOneModulatorModulePoly)
};

/** A modulator module that just outputs the current note pitch of a given voice. */
class JUCE_API rsNotePitchModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNotePitchModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    //setModuleTypeName("NotePitch");
    setModulationSourceName("NotePitch");
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    return voiceManager->getVoicePitch(voiceIndex);
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNotePitchModulatorModulePoly)
};

/** A modulator module that just outputs the current normalized velocity of a given voice. */
class JUCE_API rsNoteVelocityModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNoteVelocityModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    //setModuleTypeName("NormalizedVelocity");
    setModulationSourceName("NoteVelocity");
  }
  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    return voiceManager->getVoiceNormalizedVelocity(voiceIndex);
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNoteVelocityModulatorModulePoly)
};

/** A modulator module that just outputs the current frequency of a given voice. */
class JUCE_API rsNoteFreqModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNoteFreqModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    //setModuleTypeName("NotePitch");
    setModulationSourceName("NoteFrequency");
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    return RAPT::rsPitchToFreq(voiceManager->getVoicePitch(voiceIndex));
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNoteFreqModulatorModulePoly)
};


/** A modulator module that outputs the current value of the pitch-wheel in the range -1..+1. */
class JUCE_API rsPitchBendModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsPitchBendModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse, voiceManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    setModulationSourceName("PitchBend");
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    return voiceManager->getPitchBend();
  }
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsPitchBendModulatorModulePoly)
};


// todo: pitch-wheel, aftertouch, etc




//=================================================================================================

class JUCE_API AttackDecayEnvelopeModulePoly : public ModulatorModulePoly
{

public:

  AttackDecayEnvelopeModulePoly(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr,
    rsVoiceManager* voiceManagerToUse = nullptr);

  virtual ~AttackDecayEnvelopeModulePoly() {}

  void handleMidiMessage(MidiMessage msg) override;

  /*
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override { core.reset();  }

  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOff(key, 0);  }

  virtual double getModulatorOutputSample() override { return core.getSample(); }
  */


  //virtual void handleMidiMessageWithVoiceInfo(MidiMessage message, int voiceInfo)

  virtual void noteOnForVoice( int key, int vel, int voice) override 
  { 
    cores[voice].noteOn(key, vel); 
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceIndex <= cores.size());
    return cores[voiceIndex].getSample();
  }

  // parameter callback targets:
  void setAttack(double newAttack, int voice);
  void setDecay( double newDecay,  int voice);
  // maybe for the polyphonic case, these functions should just receive a second argument for the 
  // voice index? this will be a quite flexible design but will require to write a lot of 
  // boilerplate code to make objects polyphonic - for each modulatable parameter, we must write a 
  // boilerplate per-voice callback function

protected:


  virtual void createParameters();
  void allocateVoiceResources(rsVoiceManager* voiceManager) override;


  std::vector<RAPT::rsAttackDecayEnvelope<double>> cores;
  // todo: use a single full rsAttackDecayEnvelope object and an array of 
  // rsAttackDecayEnvelopeVoice objects that store only voice specific state and have a pointer to
  // the single "master" object - gets rid of redundant values


  // parameters:
  double sampleRate = 44100;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AttackDecayEnvelopeModulePoly)
};

// also have a module that just produces a constant value - that may be useful for exponentially
// scaled parameters - we can just subtract off their minimum 





#endif