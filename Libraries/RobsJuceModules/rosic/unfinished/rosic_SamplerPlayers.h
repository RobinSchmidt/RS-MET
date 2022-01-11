#ifndef rosic_SamplerPlayers_h  // rename files to SamplePlayers.h/cpp
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

/** Baseclass for RegionPlayer and GroupPlayer to factor out the common stuff. */

class SamplePlayer
{

public:

  /** Sets up the pool of DSP resource objects (such as filters, modulators, etc.) from which the
  player may grab items in order to prepare itself for playing back a region. This should be set
  up by the engine soon after it has created all its players. */
  void setDspResourcePool(DspResourcePool* newPool) { dspPool = newPool; }

protected:

  using PlaybackSetting = rsSamplerData::PlaybackSetting;

  /** Returns a pointer to a processor of given type, if available, otherwise a nullptr. Used in
  buildProcessingChain. */
  SignalProcessor* getProcessor(DspType type)
  {
    RAPT::rsAssert(dspPool, "This pointer should be assigned soon after creation");
    return dspPool->processorPool.grabProcessor(type);
  }

  /** Adds the DSPs of the given types to the chain of actual DSP objects, if needed. Adding a 
  particular DSP is needed, if no suitable such DSP is already there in our dspChain where 
  "suitable" means: "with right type and index". The return value informs, whether or not adding
  the desired DSPs was succesful. */
  bool addDspsIfNeeded(const std::vector<DspType>& dspTypeChain);

  /** Given a playback setting (i.e. opcode, value, possibly index), it finds the processor in our
  dspChain to which this setting applies and sets the corresponding parameter in the DSP. It 
  assumes that a suitable processor exists in our chain - if not, then something went wrong with
  building the assembling the dspChain in a step before and an assert is triggered. */
  void setupProcessorSetting(const PlaybackSetting& s);

  /** This is supposed to be overriden by subclasses to actually assemble the DSP chain they 
  need. The implementation should return true, if assembling the chain was successful and false 
  otherwise (when not enough DSPs are available). In the latter case, it is also the job of the 
  function to clean up any partially built chain if necessary. On return, the dspChain should 
  either be built completely and correctly (and true be returned) or not at all (and false be 
  returned). The post-condition should be: the dspChain is either built fully or it is empty. For
  the meaning of the busMode parameter, see rsSamplerEngine2::setBusMode. For normal sfz-like 
  behavior, it should be set to false. */
  virtual bool assembleDspChain(bool busMode) = 0;
  // -maybe use an int mode parameter later when more flexibility is needed
  // -maybe provide default argument false for busMode

  /** Reposits all the DSP objects back into the dspPool and clears our dspChain. */
  void disassembleDspChain();


  SignalProcessorChain dspChain;
  /** This is the chain of our actual DSP objects. A SignalProcessorChain is basically an array
  of pointers to polymorphic signal processor classes that can be assembled at runtime,  
  typically on noteOn. */

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
    const AudioFileStream<float>* sampleStream, double outputSampleRate, bool busMode);
  // todo: later maybe have default values (false) for the busMode 

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
  rsReturnCode prepareToPlay(double sampleRate, bool busMode);


  //bool assembleDspChain(bool withGroupDsps, bool withInstrumDsps);
  // maybe use only one bool called busMode or something...but it's a 
  // bit more complicated: some group/instrument settings should apply to RegionPlayers even in
  // busMode or drumMode - namely the settings that affect the sample playback source (delay, 
  // pitch, etc.). Maybe we should always call all 3 setup functions and pass a busMode flag 
  // through all the way down

  bool assembleDspChain(bool busMode) override;


  bool setupModulations();

  //void resetDspState();
  void resetPlayerSettings();

  // move to baseclass, if possible and/or maybe have a virtual detupDspSettings function in 
  // baseclass that we override here and in the GroupPlayer:
  void setupDspSettingsFor(const Region* r, double sampleRate, bool busMode);
  void setupDspSettings(const std::vector<PlaybackSetting>& settings,
    double sampleRate, bool busMode);



  const Region* region;                 //< The Region object that this object should play
  const AudioFileStream<float>* stream; //< Stream object to get the data from

  rsFloat64x2 amp = 1.0;         //< Amplitude (for both channels)
  // maybe remove or replace by channel-matrix gains gLL, gLR, gRL, gRR

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

  std::vector<Modulator*> modulators;
  std::vector<ModulationConnection*> modMatrix;  // not a literal matrix but conceptually


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

class GroupPlayer : public SamplePlayer
{

public:

  /** Generates one stereo sample frame at a time. */
  rsFloat64x2 getFrame();

  /** Release the resources that were acquired for playback (DSPs, RegionPlayers, etc.). */
  void releaseResources();
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


  const rsSamplerData::Group* getGroupToPlay() const { return group; }

  bool setGroupToPlay(const rsSamplerData::Group* groupToPlay, bool busMode);


  //void setGroupToPlay(const rsSamplerData::Group* groupToPlay) { group = groupToPlay; }
  // should do more stuff: assmeble and set up the DSP chain, return a bool to report success or
  // failure


protected:

  bool assembleDspChain(bool busMode) override;

  void setupDspChain();
  // maybe make this an override of a baseclass method...if possible

  std::vector<RegionPlayer*> regionPlayers;
  // Pointers to the players for all the regions in this group.


  const rsSamplerData::Group* group = nullptr;
  // Pointer to the group object which is played back by this player


  friend class rsSamplerEngine2;  
  // Try to get rid! We need it because SamplerEmgine2 calls assembleDspChain - maybe we should 
  // move it to unprotected...or: setGroupToPlay should call it - yes that seems cleaner

  //rsSamplerEngine2* engine = nullptr;
  // For communication with enclosing sampler-engine - currently not needed - try to keep it like
  // that (reduce coupling)

};

//===============================================================================================

class InstrumPlayer : public SamplePlayer
{

public:

  bool assembleDspChain(bool busMode) override;

  void releaseResources()
  {
    disassembleDspChain();
  }
  // if all 3 subclasses have such a method, introduce a virtual baseclass method

protected:

  const rsSamplerData::Instrument* instrum = nullptr;
  // Pointer to the instrument object which is played back by this player


};
// maybe the implementations can be moved into the baseclass as default implementations




}}

#endif