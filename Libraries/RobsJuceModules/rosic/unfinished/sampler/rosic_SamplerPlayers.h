#ifndef rosic_SamplerPlayers_h  // rename files to SamplePlayers.h/cpp
#define rosic_SamplerPlayers_h

namespace rosic { namespace Sampler {

//=================================================================================================

class EffectChain
{
public:

  using uchar = unsigned char;

  void processFrame(float* L, float* R);
  void processBlock(float* L, float* R, int N);
  void prepareToPlay(uchar key, uchar vel, double fs) 
  { for(auto & p : processors) p->prepareToPlay(key, vel, fs); }
  //void resetState()    { for(auto & p : processors) p->resetState();    }
  //void resetSettings() { for(auto & p : processors) p->resetSettings(); }
  //void reset() { resetState(); resetSettings(); }

  void reserve(size_t num) { processors.reserve(num); }
  void addEffect(Effect* p) { processors.push_back(p); }
  void clear() { processors.clear(); }

  bool isEmpty() const { return processors.empty(); }

  /** Returns the total number of processors in the chain. */
  size_t getNumEffects() const { return processors.size(); }

  /** Returns the number of processors of given type in the chain. */
  size_t getNumEffects(DspType type) const;


  Effect* getEffect(int i) { return processors[i]; } 
  // is this needed? it's confusing to have this and the function below because the indices mean
  // different things in both cases

  /** Returns the sfzIndex-th processor of the given type within the chain or nullptr, if there are
  not enough (i.e. less than i) processors of the given type in the chain. This is a 1-based index
  as it occurs in the sfz files. To get the 3rd filter, you would pass type = Dsp::Filter, 
  index = 3. For certain opcodes, an index is not applicable. We usually encode this by setting the
  value to -1 in the data-record. Such a -1 will then be interpreted as "first-and-only" and in 
  this case, it doesn't really matter, if the caller passes -1 or +1 into this function. */
  Effect* getEffect(DspType type, int sfzIndex);

protected:

  std::vector<Effect*> processors;  // rename to effects
  // Where is this allocated, i.e. where do we call resize on this? I mean the vector itself, not
  // the pointed-to objects. This is not yet well defined, i think. It should probably also happen
  // on sfz load

};

//=================================================================================================

class ModulationConnection
{

private:
  std::function<void(double)> targetSetter; // target callback that sets some parameter
  double amount = 0.0;  // strength of modulation
  double refVal = 0.0;  // unmodulated reference value
};

//=================================================================================================

/** UNDER CONSTRUCTION

A class used to represent the current playback status of the engine. There exists exactly one 
object of this type as member of SamplerEngine. It holds information about the most recently 
received midi events of particular kinds (such as controllers and pitchbend) to facilitate 
responding to these events in the effects etc.  ..tbc...

It is also used to hold and accumulate some intermediate values that arise in the calculation of
the Player's member variables to control certain playback aspects such as pitch, amplitude etc. 
There are often multiple opcodes that influence such a setting and we must override or accumulate 
them all seperately before we can compute the final resulting member variable (such as the 
per-sample time increment or the final amplitude scaler). To facilitate this, we pass a pointer to
such a struct to setupPlayerSetting. */

class PlayStatus
{
public:





  //-----------------------------------------------------------------------------------------------
  // the members used to facilitate accumulation of intermediate calculation results:

  // Initialize members to evil values to make sure they are set up correctly in reset:
  double transpose = 666;
  double tune = 666;
  // todo: bendUp, bendDown
  // We use double not mainly because we expect a lot of error accumulation in our computations
  // (although that may be the case as well) but rather because we can easily afford it. This 
  // struct is not used anywhere where minimizing space requirement matters anyway. Its just 
  // created on the stack on note-on and then a pointer to is is passed around to all functions 
  // that need it until the event has been fully consumed and the struct goes out of scope again.
  // We don't store arrays of these things or anything like that anywhere. ...Or maybe the engine
  // should have a member variable of the type that is re-used whenever a new note is triggered
  // ...hmmmm - we'll get a lot more settings - maybe switch to float

  //float pitch_keycenter = 666;
  float pitch_keytrack  = 666;

  std::vector<float> ampN_veltrack; // ToDo: ampN_keytrack, pitchN_veltrack, pitchN_keytrack
  // I think, these are only needed to support key/vel tracking in levels_are_busses mode. In 
  // override mode, we can get away without them. I think, they are needed, because we need to 
  // first accumulate all key- and veltrack values into final total key- and veltrack values and
  // then use the respective tracking formula on these totals. But maybe the formulas could be 
  // implemented incrementally such that each contribution to the amplitude factor gets directly
  // accumulated into the final amplitude multiplier? Then we could perhaps replace the arrays by
  // single values or even get rid of them entirely and accumulate directly into the respctive
  // parameter...we'll see...some math and design questions to figure out....this is still under 
  // construction...

private:

  PlayStatus()
  {
    int dummy = 0;
    ampN_veltrack.resize(5);  // preliminary: ToDo: in SamplerEngine::preAllocateDspMemory, allocate
    // as many as needed in the maximum case defined by the region/group/instrument that uses the 
    // largest number of amplifiers
  }
  

  void reset()
  {
    transpose       = 0.0;
    tune            = 0.0;
    //pitch_keycenter = 60.f;
    pitch_keytrack  = 100.f;

    using namespace RAPT;
    rsFill(ampN_veltrack, 0.f);
  }
  // maybe needed when we don't create it on the stack but rather re-use a member

  // under construction - not yet used:
  char  controllers[128];  // most recnetly received values of all controllers
  short pitchWheel;        // most recently received value of pitch wheel
  // todo: aftertouch, etc.
  bool  dirty = false;
  // ToDo: The dirty flag should be set to true by the midi event handlers when control, wheel, 
  // etc. messages are received, it should be inspected and used by effects to update their coeffs 
  // if needed and it should be set back to false at the end of the engine's process function 
  // because when that function exits, is assumed that all effects have updated themselves in the 
  // process. Maybe we should have additional flags controllersDirty, pitchwheelDirty etc. the 
  // global dirty flag will be true when any of theses are true, so we can inspect only the global
  // dirty flag per sample in the effects and only if it's true, inspect the other flags to figure 
  // out what kind of event happened (to avoid recalculations when they are not necessary - some
  // effects may only want to respond to controllers but not to other kinds of events). Maybe we
  // should use std::atomic<bool> to communicate the fact that it is used for communication between
  // the audio- and midi-thread even though we assume that they are actually the same thread (but 
  // maybe we shouldn't assume that anyway - in a plugin, this will generally be true but maybe 
  // within a DAW, the midi events would indeed be received in another thread?)




  friend class rsSamplerEngine; 
  // rsSamplerEngine has an object of this class as member, so it needs to be able to call the 
  // (private) constructor. This class is a bit like a singleton, but it doesn't live in the 
  // global scope but is a member of the sampler engine. It doesn't make sense to have objects 
  // of that class anywhere else. ToDo: make copy/move constructors/assignments unavailable
  // Can we restrict the friendship further to only the constructor of rsSamplerEngine such that
  // we don't accidiently do something wron inside rsSamplerEngine?

};
// -rename to SetupIntermediates and include also intermediates for the DSPs...hmm...but there's
//  an indefinite number of them...we'll see
// -Maybe this class should also contain all the information about the most recently received
//  controller, pitch_wheel etc. and maybe each SignalProcessor should maintain a pointer to the 
//  (single) object. Maybe it should be called SamplerMidiStatus. Maybe it should maintain an 
//  std::atomic<bool> dirty flag which is set to true in noteOn (or handleMidiEvent), inspected in
//  the DSPs, if true, they recalc their coeffs and we set it to to false at the end of 
//  processFrame. Some mechanism like this is needed anyway to support midi-cc.  Maybe this flag
//  could make the preparToPlay call in the DSPs superfluous...but maybe not - the preparation may 
//  actually do more than just recalc coeffs (such as resetting buffers). It doesn't realy need to
//  be atomic because we expect to receive midi on the audio thread anyway but it would kinda 
//  communicate that this flag is used as information channel from the midi to the audio side of 
//  things.

//=================================================================================================

class RegionPlayer;

/** Baseclass for RegionPlayer, GroupPlayer and InstrumPlayer to factor out the common stuff. The 
functionality for sample playback is distributed mostly between this baseclass and the subclass 
RegionPlayer. GroupPlayer and InstrumPlayer are also subclasses of SamplePlayer and their purpose 
is mostly to sum the signals of their embedded lower level players and apply additional DSP
processes to these sums. */

class SamplePlayer
{

public:

  /** Sets up the pool of DSP resource objects (such as filters, modulators, etc.) from which the
  player may grab items in order to prepare itself for playing back a region. This should be set
  up by the engine soon after it has created all its players. */
  void setDspResourcePool(DspResourcePool* newPool) { dspPool = newPool; }

protected:


  /** Returns a pointer to an effect of given type, if available, otherwise a nullptr. Used in
  assembleEffectChain. */
  Effect* getEffect(DspType type)
  {
    RAPT::rsAssert(dspPool, "This pointer should be assigned soon after creation");
    return dspPool->effectPool.grabEffect(type);
  }

  /** Adds the DSPs of the given types to the chain of actual DSP objects, if needed. Adding a 
  particular DSP is needed, if no suitable such DSP is already there in our dspChain where 
  "suitable" means: "with right type and index". The return value informs, whether or not adding
  the desired DSPs was succesful. It may fail due to not having enough DSPs of required types 
  available. In such cases, any partially assembled dspChain will be disassembled again and 
  false is returned. This potential disassembly is what is meant by the "or clean-up" part. */
  bool augmentOrCleanEffectChain(const std::vector<DspType>& dspTypeChain);

  /** This is supposed to be overriden by subclasses to actually assemble the DSP chain they 
  need. The implementation should return true, if assembling the chain was successful and false 
  otherwise (when not enough DSPs are available).  */
  virtual bool assembleEffectChain(bool busMode) = 0;
  // rename to assemberEffectChain
  // -maybe use an int mode parameter later when more flexibility is needed
  // -maybe provide default argument false for busMode

  /** A helper function that is called from GroupPlayer::assembleDspChain(bool) and
  InstrumentPlayer::assembleDspChain(bool). */
  bool assembleEffectChain(const std::vector<DspType>& dspTypes);

  /** Reposits all the DSP objects back into the dspPool and clears our dspChain. */
  void disassembleEffectChain();




  using PlaybackSetting = SfzInstrument::PlaybackSetting; // for convenience
  // mayb rename to OpcodeSetting

  /** Given a playback setting (i.e. opcode, value, possibly index) that is supposed to be 
  applicable to the sample playback source, the overriden version of this function in the 
  subclasses manipulate the corresponding state of the RegionPlayer, i.e. the lowest level 
  player that may be managed by higher level players, accordingly. RegionPlayer itself sets
  up its own member variables, GroupPlayer manipulates one of its embedded RegionPlayers, 
  etc.  */
  virtual void setupPlayerSetting(const PlaybackSetting& s, double sampleRate, 
    RegionPlayer* rp, PlayStatus* iv) = 0;
  // rename to setPlayerOpcode

  /** Given a playback setting (i.e. opcode, value, possibly index) that is supposed to be 
  applicable to the DSP chain, it finds the processor in our dspChain member to which this 
  setting applies and sets the corresponding parameter in the DSP. It assumes that a suitable 
  processor exists in our chain - if not, then something went wrong with assembling the 
  dspChain in a step before and an assert is triggered. */
  virtual void setupProcessorSetting(const PlaybackSetting& s);
  // rename to setDspOpcode or setEffectOpcode...but maybe it applies to modulators, too,
  // maybe the more general "processor" is actually appropriate

  /** ToDo: add documentation */
  virtual void setupDspSettings(const std::vector<PlaybackSetting>& settings, double sampleRate, 
    RegionPlayer* rp, bool busMode, PlayStatus* iv);
  // Maybe return a bool to indicate, if the setting was handled (if false, the subclass may
  // want to do something in its override)...??? comment obsolete?


  EffectChain effectChain;
  /** This is the chain of our effect processor objects. An EffectChain is basically an array
  of pointers to polymorphic effect classes (i.e. subclasses of the Effect baseclass) that can be 
  assembled at runtime, typically on noteOn. */

  DspResourcePool* dspPool = nullptr;
  /**< A pool of DSP processor objects from which we can grab some as needed to assemble our DSP
  chain. The assembly task is mostly done in the subclasses making use of addDspsIfNeeded. The 
  pointer should be set by the engine once and for all when it creates its Players. */

};

//=================================================================================================

/** A class for playing back a given Region object. */

class RegionPlayer : public SamplePlayer
{

public:

  // For convenience:
  using uchar = unsigned char;
  using Region = SfzInstrument::Region; // todo: make a subclass here that adds the stream field
  using Group  = SfzInstrument::Group;
  using PlaybackSetting = SfzInstrument::PlaybackSetting;

  /** Sets up the region object that this player should play. You need to also pass the output
  sample-rate which is the sample rate at which the player should run (not the sample rate of the
  audio file associated with the region). The groupSettingsOverride parameter determines whether
  the group settings should override the instrument settings (true) or be applied on top of them
  (false). Likewise, regionSettingsOverride = true lets the region settings override the group
  settings. */
  rsReturnCode setRegionToPlay(const Region* regionToPlay, 
    const AudioFileStream<float>* sampleStream, uchar key, uchar vel, double outputSampleRate, 
    bool busMode, PlayStatus* iv);
  // todo: later maybe have default values (false) for the busMode 

  const Region* getRegionToPlay() const { return region; }

  /** Sets the midi note number for which this player was started. This needs to be set up when
  receiving a noteOn. This information is used later when receiving a noteOff to identify which
  players need to stop. */
  //void setKey(uchar newKey) { key = newKey; }

  /** Generates one stereo sample frame at a time. */
  void processFrame(float* L, float* R);
  // use float pointers for signals consistently

  /** Writes a block of given length into the outBuffer. */
  void processBlock(float* L, float* R, int length);

  /** Returns true, iff this player has finished its playback job, for example by having reached
  the end of the sample stream and/or amplitude envelope. */
  bool hasFinished(); // should be const?

  /** Retrieves the information about the midi note for which this player was started. Used to
  identify players that need to stop, when a noteOff is received. @see setKey */
  uchar getKey() const { return key; }

  /** Returns the loop mode that this player is in. */
  LoopMode getLoopMode() const { return loopMode; }


  /** Releases all resources that were acquired for playback such as signal processors, 
  modulators, etc. */
  void releaseResources();


  /** Allocates memory for the pointers to the processors, modulators and modulation 
  connections. Should be called soon after creation. */
  void allocateMemory();


protected:

  /** A basic sanity check for the given region. Mostly for catching bugs. */
  bool isPlayable(const Region* region);

  /** Acquires the ressources that are required for the playback of the given region such as the
  required DSPs and modulators, sets up the internal values for the playback settings (including 
  DSP objects) according to the assigned region and resets all DSP objects. When it returns
  rsReturnCode::success, it means the player is now ready to play. If it returns anything else,
  it means that something went wrong  - presumably not enough ressources were available - and the
  engine should discard the player object, i.e. put it back into the pool. */
  rsReturnCode prepareToPlay(uchar key, uchar vel, double sampleRate, bool busMode,
    PlayStatus* iv);


  bool assembleEffectChain(bool busMode) override;


  bool setupModulations();

  //void resetDspState();
  void resetPlayerSettings();

  // move to baseclass, if possible and/or maybe have a virtual detupDspSettings function in 
  // baseclass that we override here and in the GroupPlayer:
  void setupDspSettingsFor(const Region* r, double sampleRate, bool busMode, 
    PlayStatus* iv);

  void setupFromIntemediates(const PlayStatus& iv, double sampleRate);

  void setupPlayerSetting(const PlaybackSetting& s, double sampleRate, 
    RegionPlayer* rp, PlayStatus* iv) override;

  const Region* region;                 //< The Region object that this object should play
  const AudioFileStream<float>* stream; //< Stream object to get the data from

  // Maybe we should use some sort of fixed-point format for this instead?
  double sampleTime = 0.0;  //< Time index in the sample. Negative values used for delay.
  double increment  = 1.0;  //< Increment of sampleTime per sample
  double loopStart  = 0;    //< Start of loop in samples
  double loopEnd    = 0;    //< End of loop in samples



  float  endTime    = 0;
  // Maybe it should be int? Or maybe double to ease comparison? I actually think, we should get
  // rid of it and inquire the info from the stream. We need one double -> int conversion anyway
  // for the interpolator. The result of that can be compared to the length of the stream

  float  offset     = 0;    //< Start offset into the sample
  // maybe rename to startTime or startSample

  float  amplitude  = 1;    
  // This amplitude factor is determined by key/vel crossfades, etc. but not by key/veltrack 
  // because these are is done by the Amplifier DSP



  LoopMode loopMode = LoopMode::no_loop;
  uchar key = 0;                 //< Midi note number used for starting this player

  std::vector<Modulator*> modulators;
  std::vector<ModulationConnection*> modMatrix;  // not a literal matrix but conceptually




  friend class SampleBusPlayer;
  // So it can accumulate the group and instrument settings into our increment, sampleTime,
  // etc. variables. (ToDo: maybe provide functions applyAdditionalDelay, applyAdditionalDetune 
  // later and unfriend the SamplePlayer again)


  // ToDo: 
  // -arrange members to avoid padding to minimize memory footprint of this object
  // -maybe use an int and a float to represent integer and fractional parts of sampleTime 
  //  separately (we need then also 2 increments). This has two advantages:
  //  -we don't need to compute the fractional part at each sample
  //  -the pitch accuracy will not degrade for loops later in the sample
  // -maybe use float instead of uchar for the key
  // -add more DSP objects: envelope generators, LFOs, we also need to somehow take care
  //  of the effect sends (reverb, chorus)
  // -Maybe in addition to elapsed time, we also need to maintain position inside the sample 
  //  stream. This is not the same thing due to possible delay and looping.
  // -maybe use a different implementation structure (SVF) for the time-varying filter
  // -try to optimize ram and/or cpu usage by re-ordering
  // -ways to save space (if that should be necessary): use float instead of double, isnteda of
  //  using std::vector, use a hand-rolled dynamic array of indices (could be short int)

};

//===============================================================================================

class SampleBusPlayer : public SamplePlayer
{

public:

  using uchar = unsigned char;

  void setupPlayerSetting(const PlaybackSetting& s, double sampleRate, 
    RegionPlayer* rp, PlayStatus* iv) override;

  bool setGroupOrInstrumToPlay(const SfzInstrument::HierarchyLevel* thingToPlay, 
    uchar key, uchar vel, double sampleRate, RegionPlayer* regionPlayer, bool busMode, 
    PlayStatus* intermediates);
  // busMode is superfluous - when a SampleBusPlayer is invoked, we are in busMode by definition

  virtual void releaseResources()
  {
    disassembleEffectChain();
    grpOrInstr = nullptr;
  }

protected:

  bool assembleEffectChain(bool busMode) override;

  const SfzInstrument::HierarchyLevel* grpOrInstr = nullptr;
  // pointer to the group or isntrument that this player should play

};

//===============================================================================================

/** A class for collecting all the SignalProcessors that apply to a given group. This is used
only when the group's DSP settings should go on top of the region's settings, i.e. in 
"drum-sampler" mode where the groups map to sub-busses and the instrument maps to the master 
bus. */

class GroupPlayer : public SampleBusPlayer
{

public:

  /** Generates one stereo sample frame at a time. */
  void processFrame(float* L, float* R);

  // implement processBlock

  /** Release the resources that were acquired for playback (DSPs, RegionPlayers, etc.). */
  void releaseResources() override;
  // make virtual method in baseclass

  /** Adds a new region player to our regionPlayers array. */
  void addRegionPlayer(RegionPlayer* newPlayer);

  /** Removes the given player from our regionPlayers array. */
  void removeRegionPlayer(RegionPlayer* player);

  /** Returns true, iff the given regionPlayer is part of this GroupPlayer, i.e.  */
  bool contains(RegionPlayer* rp) { return RAPT::rsContains(regionPlayers, rp); }
  // ToDo: make parameter rp const - for some reason, it doesn't compile

  /** Returns true, iff this GroupPlayer has no RegionPlayer objects running. */
  bool hasNoRegionPlayers() { return regionPlayers.empty(); }

  /** Returns a pointer to the group to play. If there is currently no group assigned to this 
  player, this is a nullptr. */
  const SfzInstrument::Group* getGroupToPlay() const 
  { return (const SfzInstrument::Group*) grpOrInstr; }

  /** Sets the group that should be played back by this player. */
  bool setGroupToPlay(const SfzInstrument::Group* groupToPlay, uchar key, uchar vel, 
    double fs, RegionPlayer* rp, bool busMode, PlayStatus* iv)
  { return setGroupOrInstrumToPlay(groupToPlay, key, vel, fs, rp, busMode, iv); }
    // ...it's just a convenience function to make the call site look nicer.

protected:

  std::vector<RegionPlayer*> regionPlayers;
  // Pointers to the players for all the regions in this group.

};

//===============================================================================================

class InstrumPlayer : public SampleBusPlayer
{

public:

  void addRegionPlayer(RegionPlayer* newPlayer);

  void processFrame(float* L, float* R) { effectChain.processFrame(L, R);  }

  // implement processBlock

  bool setInstrumToPlay(const SfzInstrument::Global* instrumToPlay, uchar key, uchar vel, 
    double fs, RegionPlayer* rp, bool busMode, PlayStatus* iv)
  { return setGroupOrInstrumToPlay(instrumToPlay, key, vel, fs, rp, busMode, iv); }
    // Convenience function to make the call site look nicer.

};

}}

#endif