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

// some very simple modulators - maybe they should be moved into the jura_framework module

/** A modulator module that just outputs the constant value 1. That may be useful for use with 
exponentially scaled parameters such as a frequency that goes from 20 to 20000. We can just 
subtract 20 by connecting the constant with depth -20 and the use additive envelopes to directly
determine the freq without having the minimum getting in the way. */
class JUCE_API rsConstantOneModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsConstantOneModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    setModulationSourceName("ConstantOne");
  }

  double renderModulation()                    override { return 1.0; } 
  double renderVoiceModulation(int voiceIndex) override { return 1.0; }

  void allocateVoiceModResources() override {}

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsConstantOneModulatorModulePoly)
};

/** A modulator module that just outputs the current note pitch of a given voice. */
class JUCE_API rsNotePitchModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNotePitchModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    //setModuleTypeName("NotePitch");
    setModulationSourceName("NotePitch"); // maybe rename to NoteKey, MidiKey
  }

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    return voiceManager->getVoicePitch(voiceIndex);

    // ToDo: maybe use something like:
    // smoother[voiceIndex]->getSample( voiceManager->getVoicePitch(voiceIndex) );
    // i think, this will also provide glide

  }

  void allocateVoiceModResources() override {}

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNotePitchModulatorModulePoly)
};


/** A modulator module that just outputs the current normalized velocity of a given voice. */
class JUCE_API rsNoteVelocityModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNoteVelocityModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
  {
    ScopedLock scopedLock(*lock);
    //setModuleTypeName("NormalizedVelocity");
    setModulationSourceName("NoteVelocity");
  }

  double renderModulation() override
  {
    jassert(voiceManager != nullptr);
    int voiceIndex = voiceManager->getNewestVoice();
    double vel = voiceManager->getVoiceNormalizedVelocity(voiceIndex); // for debug
    return voiceManager->getVoiceNormalizedVelocity(voiceIndex);
  }
  // not called when we use AcidDevil and try to modulate its cutoff via velocity -> check, if 
  // ToolChain calls applyModulations

  double renderVoiceModulation(int voiceIndex) override
  {
    jassert(voiceManager != nullptr);
    double vel = voiceManager->getVoiceNormalizedVelocity(voiceIndex); // for debug
    return voiceManager->getVoiceNormalizedVelocity(voiceIndex);
  }

  void allocateVoiceModResources() override {}


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNoteVelocityModulatorModulePoly)
};

/** A modulator module that just outputs the current frequency of a given voice. */
class JUCE_API rsNoteFreqModulatorModulePoly : public ModulatorModulePoly
{
public:

  rsNoteFreqModulatorModulePoly(CriticalSection* lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr,
    ModulationManager* modManagerToUse = nullptr)
    : ModulatorModulePoly(lockToUse, metaManagerToUse, modManagerToUse) 
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

  void allocateVoiceModResources() override {}

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsNoteFreqModulatorModulePoly)
};


/** A modulator module that outputs the current value of the pitch-wheel in the range -1..+1. */
/*
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

*/

// todo: controller, aftertouch, etc. - maybe these modulators should have a built-in smoothing, 
// maybe independent form meta-smoothing system. in allocateVoiceModResources, just allocate a 
// bunch of smoothing-filters. Maybe have a baseclass SmoothedMidiModulator from which smoothed
// versions of the pitch/vel/etc modulators derive




//=================================================================================================

class JUCE_API AttackDecayEnvelopeModulePoly : public ModulatorModulePoly
{

public:

  AttackDecayEnvelopeModulePoly(CriticalSection *lockToUse,
    MetaParameterManager* metaManagerToUse = nullptr, 
    ModulationManager* modManagerToUse = nullptr);

  virtual ~AttackDecayEnvelopeModulePoly() {}

  void handleMidiMessage(MidiMessage msg) override;

  /*
  virtual void setSampleRate(double newSampleRate) override; 
  virtual void reset() override { core.reset();  }

  virtual void noteOn( int key, int vel) override { core.noteOn(key, vel); }
  virtual void noteOff(int key)          override { core.noteOff(key, 0);  }

  virtual double getModulatorOutputSample() override { return core.getSample(); }
  */

  void noteOn(int key, int vel, int voice) override 
  { 
    jassert(voice >= 0 && voice < voiceManager->getMaxNumVoices());
    // I think getMaxNumVoices makes more sense than getNumVoices because a note-off can also
    // be passed as note-on (with zero velocity) and if the user reduces the number of voices
    // between note-on and note-off, the condition could occur that we receive a note-off for a 
    // non-existent voice. ...but no: the dispatche already makes sure that only *actual* 
    // note-ons are passed to this - and those should really be constrained to the current
    // "NumVoices" setting...hmmm...figure out....

    cores[voice].reset();
    // maybe do this only if the voice is inactive - if it's in release, don't reset
    // i think, the voice manager needs an optional reuseReleasingVoice mode - maybe call it
    // retriggerMode with options: useNewVoice, reuseReleasingVoice - are there other meaningful 
    // options? maybe reuseIfBelow, reuseIfAbove, etc...

    cores[voice].noteOn(key, vel);
  }

  // what about noteOffForVoice? if the envelope has sustain, we need to override that, too

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
  void allocateVoiceModResources() override;


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