#ifndef rosic_SamplerPlayers_h
#define rosic_SamplerPlayers_h

namespace rosic { namespace Sampler {

//=================================================================================================

class SignalProcessorChain
{
public:

  void processFrame(rsFloat64x2& inOut);
  void processBlock(rsFloat64x2* inOut, int N);
  void prepareToPlay(double fs) { for(auto & p : processors) p->prepareToPlay(fs); }
  //void resetState()    { for(auto & p : processors) p->resetState();    }
  //void resetSettings() { for(auto & p : processors) p->resetSettings(); }
  //void reset() { resetState(); resetSettings(); }

  void reserve(size_t num) { processors.reserve(num); }
  void addProcessor(SignalProcessor* p) { processors.push_back(p); }
  void clear() { processors.clear(); }

  bool isEmpty() const { return processors.empty(); }

  /** Returns the total number of processors in the chain. */
  size_t getNumProcessors() const { return processors.size(); }

  /** Returns the number of processors of given type in the chain. */
  size_t getNumProcessors(DspType type) const;


  SignalProcessor* getProcessor(int i) { return processors[i]; }

  /** Returns the index-th processor of the given type within the chain or nullptr, if there are
  not enough (i.e. less than i+1, counting starts at zero) processors of the given type in the 
  chain. To get the 3rd filter, you would pass type = SignalProcessorType::Filter, index = 2.  */
  SignalProcessor* getProcessor(DspType type, int index);

protected:

  std::vector<SignalProcessor*> processors;

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

/** A class for playing back a given Region object. */

class RegionPlayer 
{

public:

  // For convenience:
  using uchar = unsigned char;
  using Region = rsSamplerData::Region; // todo: make a subclass here that adds the stream field
  using Group  = rsSamplerData::Group;
  using PlaybackSetting = rsSamplerData::PlaybackSetting;

  /** Sets up the region object that this player should play. You need to also pass the output
  sample-rate which is the sample rate at which the player should run (not the sample rate of the
  audio file associated with the region). The groupSettingsOverride parameter determines whether
  the group settings should override the instrument settings (true) or be applied on top of them
  (false). Likewise, regionSettingsOverride = true lets the region settings override the group
  settings. */
  rsReturnCode setRegionToPlay(const Region* regionToPlay, 
    const AudioFileStream<float>* sampleStream, double outputSampleRate,
    bool groupSettingsOverride, bool regionSettingsOverride);
  // todo: later maybe have default values (false) for the settingsOnTop variables for 
  // convenience - but for implementing the signal-flow stuff, it makes sense to enforce the 
  // caller to apps a value
  // change API to take group/regionSettingsAccumulate as parameters

  const Region* getRegionToPlay() const { return region; }

  /** Sets the midi note number for which this player was started. This needs to be set up when
  receiving a noteOn. This information is used later when receiving a noteOff to identify which
  players need to stop. */
  void setKey(uchar newKey) { key = newKey; }

  /** Generates one stereo sample frame at a time. */
  rsFloat64x2 getFrame();

  /** Writes a block of given length into the outBuffer. */
  void processBlock(rsFloat64x2* outBuffer, int length);

  /** Returns true, iff this player has finished its playback job, for example by having reached
  the end of the sample stream and/or amplitude envelope. */
  bool hasFinished(); // should be const?

  /** Retrieves the information about the midi note for which this player was started. Used to
  identify players that need to stop, when a noteOff is received. @see setKey */
  uchar getKey() const { return key; }

  /** Sets up the pool of DSP resource objects (such as filters, modulators, etc.) from which the
  player may grab items in order to prepare itself for playing back a region. This should be set
  up by the engine soon after it has created all its players. */
  void setDspResourcePool(DspResourcePool* newPool) { dspPool = newPool; }

  /** Releases all resources that were acquired for playback such as signal processors, 
  modulators, etc. */
  void releaseDspObjects();

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
  rsReturnCode prepareToPlay(double sampleRate, bool groupSettingsOverride, 
    bool regionSettingsOverride);
  // change API to take group/regionSettingsAccumulate as parameters

  /** Returns a pointer to a processor of given type, if available, otherwise a nullptr. Used in
  buildProcessingChain. */
  SignalProcessor* getProcessor(DspType type);

  bool buildProcessingChain(bool withGroupDsps, bool withInstrumDsps);
  bool setupModulations();

  //void resetDspState();
  void resetPlayerSettings();
  void setupDspSettingsFor(const Region* r, double sampleRate, bool groupSettingsOverride,
    bool regionSettingsOverride);
  void setupDspSettings(const std::vector<PlaybackSetting>& settings,
    double sampleRate, bool overrideOldSetting);
  void setupProcessorSetting(const PlaybackSetting& s);


  // see comment at prepareToPlay - maybe make onTop default to false
  // change API: replace onTop with override

  const Region* region;                 //< The Region object that this object should play
  const AudioFileStream<float>* stream; //< Stream object to get the data from

  rsFloat64x2 amp = 1.0;         //< Amplitude (for both channels)
  //int sampleTime = 0;          //< Elapsed time in samples, negative values used for delay
  double sampleTime = 0.0;       //< Time index in the sample. Negative values used for delay.
  double increment  = 1.0;       //< Increment of sampleTime per sample

  // new, under construction:
  float offset    = 0;  // maybe rename to startTime or startSample
  float endTime   = 0;
  float loopStart = 0;
  float loopEnd   = 0;
  uchar loopMode  = 0;  // use an enum class with None


  uchar key = 0;                 //< Midi note number used for starting this player

  // ToDo: arrange members to avoid padding to minimize memory footprint of this object

  std::vector<Modulator*> modulators;
  std::vector<ModulationConnection*> modMatrix;  // not a literal matrix but conceptually
  SignalProcessorChain dspChain;

  DspResourcePool* dspPool = nullptr;
  /** This should be set by the engine once and for all when it creates it RegionPlayers. */

  // We may need a state, too. Can be attack/decay/sustain/release. Or maybe just play/release?
  // Or maybe no state at all but the triggerRelease just triggers the release of all envelopes?

  friend class rsSamplerEngine; // try to get rid!


  //using Biquad = RAPT::rsBiquadDF1<rsFloat64x2, double>; // todo: use TDF2
  //Biquad flt, eq1, eq2, eq3;     //< Filter and equalizers

  // ToDo: 
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
  // -env-generators need as state variables: stage (int), time-into-stage (float/double),
  //  LFO need just phase
  // -Why are so many functions declared virtual? there are not supposed to be any subclasses.
  //  -> remove the virtual declarations
};

//===============================================================================================

/** A class for collecting all the SignalProcessors that apply to a given group. This is used
only when the group's DSP settings should go on top of the region's settings, i.e. in 
"drum-sampler" mode where the groups map to sub-busses and the instrument maps to the master 
bus. */

class GroupPlayer
{

public:

  /** Generates one stereo sample frame at a time. */
  rsFloat64x2 getFrame();

  /** Resets the internal state. */
  void reset();

  /** Adds a new region player to our regionPlayers array. */
  void addRegionPlayer(RegionPlayer* newPlayer);

  /** Removes the given player from our regionPlayers array. */
  void removeRegionPlayer(RegionPlayer* player);

  /** Returns true, iff the given regionPlayer is part of this GroupPlayer, i.e.  */
  bool contains(RegionPlayer* rp) { return RAPT::rsContains(regionPlayers, rp); }
  // ToDo: make parameter rp const - for some reason, it doesn't compile

  /** Returns true, iff this GroupPlayer has no RegionPlayer objects running. */
  bool hasNoRegionPlayers() { return regionPlayers.empty(); }

  bool buildDspChain(); // maybe rename to assembleDspChain
  void clearDspChain(); // maybe rename to disassembleDspChain, teardown, clearDspChain

protected:

  std::vector<RegionPlayer*> regionPlayers;
  // Pointers to the players for all the regions in this group.

  SignalProcessorChain dspChain;
  // The chain of additional per-group signal processors that apply to the group as a whole.

  const rsSamplerData::Group* group = nullptr;
  // Pointer to the group object which is played back by this player


  //class rsSamplerEngine2;
  friend class rsSamplerEngine2;

  rsSamplerEngine2* engine = nullptr;
  // Needed for communication channel with enclosing sampler-engine..can we get rid of this?


};





}}

#endif