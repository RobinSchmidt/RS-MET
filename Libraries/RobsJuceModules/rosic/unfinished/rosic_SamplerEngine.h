#ifndef rosic_SamplerEngine_h
#define rosic_SamplerEngine_h

namespace rosic
{

//=================================================================================================

/** A class for representing musical events such as note-on/off etc. Think of it as a class to 
represent MIDI events but with some of its anachronistic restrictions lifted, such as the abysmal 
resolution of values (typically 7 or 14 bit). The template parameter T is supposed to be either 
float or double. For easy conversion and compatibility with MIDI, we still follow the (now 
historical) convention that values are in the range 0..127, but now with much higher resolution 
due to the floating point representation. So, in a nutshell, this is a class for MIDI events but 
with higher resolution for all the values. */

template<class T>
class rsMusicalEvent
{

public:

  enum class Type
  {
    noteOn,
    noteOff,
    controlChange,
    pitchWheel,
    reset
    // ...tbc...
  };

  rsMusicalEvent(Type eventType, T value1, T value2)
    : type(eventType), val1(value1), val2(value2) {}

  Type getType() const { return type; }

  T getValue1() const { return val1; }

  T getValue2() const { return val2; }

protected:

  Type type;  // e.g. noteOn/Off, controlChange, pitchWheel
  T    val1;  // e.g. key, controller number, pitchWheelMSB
  T    val2;  // e.g. velocity, controller value, pitchWheelLSB

};
// maybe move the class elsewhere for more general use - maybe it should go into rapt due to the
// templatized nature

//=================================================================================================

/** Under Construction. Not yet ready for general use. 

A sampler engine whose feature set roughly resembles the sfz specification. It's not necessarily 
meant to be feature-complete (certainly not yet) and on the other hand, it may introduce additional
features, but sfz is the spec after which this engine is roughly modeled. 

An instrument definition in sfz is organized in 3 levels of hierarchy. At the lowest level is the 
"region" which defines which sample file should be played along with a bunch of performance 
parameters such as the key- and velocity ranges, at which the sample should be played, its volume, 
pan, filter and envelope settings and a bunch of other stuff. One level higher is the "group" which 
defines common settings that apply to all regions within the given group. Groups allow to edit the 
performance parameters of multiple regions at once: If a region does not define a particular 
performance parameter, the value of the enclosing group will be used. Region specific settings, if 
present, override the group settings (i think - verify!). At the highest level is the whole 
"instrument" itself. Just like groups provide fallback settings for regions, the whole instrument 
can provide fallback settings for all the groups it contains. If some performance parameter isn't 
defined anywhere (neither in the instrument, group or region), a neutrally behaving default value 
will be used, which means the corresponding feature is not used, i.e. bypassed.  */

class rsSamplerEngine
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  rsSamplerEngine(int maxNumLayers = 16);
  // todo: use some higher default value - what is reasonable here needs some testing in realistic
  // scenarios

  virtual ~rsSamplerEngine();


  // for convenience:
  using uchar = unsigned char;
  using Region = rsSamplerData::Region; // todo: make a subclass here that adds the stream field
  using Group  = rsSamplerData::Group;
  using PlaybackSetting = rsSamplerData::PlaybackSetting;


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the maximum number of layers/regions that can be played simultaneously. */
  virtual void setMaxNumLayers(int newMax);

  /** Clears the sfz instrument definition and the samplePool */
  void clearInstrument();

  /** Adds a new sample to our pool of samples. After the sample has been added, regions can be 
  defined that make use of it. */
  int addSampleToPool(float** data, int numFrames, int numChannels, float sampleRate, 
    const std::string& path);
  // Maybe rename to addSample, it should return the index of the sample in the sample-pool
  // maybe make a struct SampleMetaData containing: numFrames, numChannels, sampleRate, rootKey
  // todo: take reference to a metaData object

  /** Loads a sample represented by the given path into our samplePool. The path is supposed to be 
  relative to some fixed root directory which is typcally the directory in which sfz file resides,
  but later this may be switched to some user and/or factory content directory, too (this requires
  to introduce a new opcode to sfz...maybe root_dir or sample_directory or sample_folder or 
  something which should be defined once for the whole instrument). It returns the following
  return codes: 
    success:       sample was succesfully loaded into the pool
    nothingToDo:   sample was already in the pool
    fileLoadError: sample could not be loaded (maybe the path was wrong?)
    memAllocFail:  we could not allocate enough memory to load the sample
  ToDo: verify return codes in unit test  */
  int loadSampleToPool(const std::string& path);

  /** All regions that use the sample with given index will be assigned to use no sample at all. */
  int unUseSample(int sampleIndex);
  int unUseSample(const std::string& path);

  /** Adds a new group to the instrument definition and returns the index of the group. */
  int addGroup() { return sfz.addGroup(); }

  /** Adds a new region to the group with the given index and returns the index of the region 
  within the group or rsReturnCode::invalidIndex, if the passed groupIndex was invalid. If the key 
  range is already known, it makes sense to pass it using the optional loKey/hiKey parameters. This
  can also be set up later, but some memory operations can be saved, if it's known in advance. */
  int addRegion(int groupIndex, uchar loKey = 0, uchar hiKey = 127);

  /** Removes the region with given group- and region index from the instrument. If a note is 
  currently playing that makes use of this region. it will continue to play as is (with the region)
  but the next time it's triggered, the region will not be part of it anymore. */
  int removeRegion(int groupIndex, int regionIndex);

  /** Sets the sample to be used for the given region within the given group. Returns either
  rsReturnCode::success or rsReturnCode::invalidIndex, if the pair of group/region indices and/or the
  sample index was invalid. */
  int setRegionSample(int groupIndex, int regionIndex, int sampleIndex); 

  /** Sets a value for a given type of playback setting a region. Returns either
  rsReturnCode::success or rsReturnCode::invalidIndex, if groupIndex and/or regionIndex was invalid. If 
  this happens, it indicates a bug on the call site. */
  int setRegionSetting(int groupIndex, int regionIdex, PlaybackSetting::Type type, float value);

  int setGroupSetting(int groupIndex, PlaybackSetting::Type type, float value);






  // todo: setGroupSetting, setInstrumentSetting, removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,

  /** Sets up the engine from the given sfz data object and returns rsReturnCode::success, if all
  is well or...  */
  int setupFromSFZ(const rsSamplerData& sfz);


  /** Writes the current instrument definition into an sfz file with given path. */
  bool saveToSFZ(const char* path) const { return sfz.saveToSFZ(path); }
  // -document, whether path is absolute or relative and if the latter, what is the root
  // -return a return-code instead of bool
  // -maybe move elsewhere

  /** Loads the instrument definition given by an sfz file with the given path. Returns 
  rsReturnCode::success if all wen well or rsReturnCode::fileLoadError if loading of the sfz or any of 
  the used samples has failed. */
  int loadFromSFZ(const char* path);
  // ToDo: In case of failure, maybe return a more specific error code/object, indicating, which 
  // file(s) exactly failed to load. This is an information that may be eventually displayed to the
  // user on a GUI.


  /** Sets the sample-rate, at which this engine should operate. This change will affect only 
  RegionPlayer objects that were started after calling this function. It's supposed to be called in
  a suspended state anyway, not in the middle of the processing. */
  void setSampleRate(double newRate) { sampleRate = newRate; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of groups in the instrument. */
  int getNumGroups() const { return sfz.getNumGroups(); }

  /** Returns the number of regions in the group with given groupIndex. */
  int getNumRegions(int groupIndex) const { return sfz.getNumRegions(groupIndex); }
  // todo: maybe assert the groupIndex is valid - if not, return invalidIndex

  /** Returns a pointer to the region object with the given group- and region index or a nullptr 
  if the combination of indices is invalid. If the client wants to edit the region, it can do so 
  only by using appropriate region-editing functions of the rsSamplerEngine object from which it 
  has requested the region pointer, for example:

    rsSamplerEngine::Region* r = se.getRegion(gi, ri);   // se is the sampler engine object
    using PST = rsSamplerEngine::PlaybackSetting::Type;
    se.setRegionSetting(r, PST::PitchKeyCenter, 69.f);
  
  The idea is that client code should not modify the settings of a Region object behind the sampler
  engine's back because the engine may have to take additional actions when certain aspects of a 
  region change. This is actually enforced by the fact that Region provides no public setters. */
  Region* getRegion(int groupIndex, int regionIndex);

  const Region* getRegionConst(int gi, int ri) const { return sfz.getRegion(gi, ri); }

  /** Returns the number of regions in the instrument definition that use the sample with the given
  index in our samplePool or rsReturnCode::invalidIndex, if the given sampleIndex is invalid. */
  int getNumRegionsUsing(int sampleIndex) const;
  int getNumRegionsUsing(const std::string& samplePath) const;

  // getGroup, getRegion, getStateAsSFZ, isSampleInPool, getNumGroups, getNumRegionsInGroup(int)
  // getNumRegions(), getSampleIndex(const string& uniqueName) ..or maybe it should take a pointer
  // to a SampleMetaData object

  /** Returns true, iff the given group index is valid, i.e. >= 0 and < numGroups. */
  bool isGroupIndexValid(int i) const { return sfz.isGroupIndexValid(i); }

  /** Returns true, iff the given pair of group- and region index is valid, i.e. a region with this
  pair of indices actually exists in the current instrument definition. */
  bool isIndexPairValid(int groupIndex, int regionIndex) const
  { return sfz.isIndexPairValid(groupIndex, regionIndex); }

  /** Returns true, iff the given sample index is valid, i.e. a sample with this index actually 
  exists our sample pool. */
  bool isSampleIndexValid(int sampleIndex) const
  { return samplePool.isSampleIndexValid(sampleIndex); }

  /** Returns the index of the sample represented by the given string in our sample pool or -1, if
  the sample is not in the pool. */
  int findSampleIndexInPool(const std::string& sample) const;

  /** Returns true, if the sample represented by the given string (as relative path with respect to
  some root directory) is present in our samplePool. */
  bool isSampleInPool(const std::string& sample) const
  { return findSampleIndexInPool(sample) != -1; }

  /** Returns the maximum number of layers that can play simultaneously. */
  int getMaxNumLayers() const { return (int) playerPool.size(); }

  /** Returns the number of currently playing layers. */
  int getNumActiveLayers() const { return (int) activePlayers.size(); }

  /** Returns the number of layers that are currently not playing, i.e. still available for adding
  a new layer to the playback. */
  int getNumIdleLayers() const { return (int) idlePlayers.size(); }

  /** Returns the number of samples that are currently in our samplePool. */
  int getNumSamples() const { return samplePool.getNumSamples(); }

  /** Returns the number of samples that were loaded to the sample pool in the most recent call to
  setupFromSFZ of loadFromSFZ. */
  int getNumSamplesLoaded() const { return numSamplesLoaded; }

  /** Returns the number of samples that were removed from the sample pool in the most recent call 
  to setupFromSFZ of loadFromSFZ. */
  int getNumSamplesRemoved() const { return numSamplesRemoved; }

  /** Returns the number of samples that failed to load to the sample pool in the most recent call 
  to setupFromSFZ of loadFromSFZ. */
  int getNumSamplesFailed() const { return numSamplesFailed; }

  /** Returns a const pointer to the rsSamplerData object that represents the current instrument
  settings. */
  const rsSamplerData& getInstrumentData() const { return sfz; }
    





  //-----------------------------------------------------------------------------------------------
  // \name Processing



  virtual void processFrame(double* left, double* right);

  virtual void processFrame(float* left, float* right);

  virtual void processBlock(float** block, int numFrames);

  // void processFrameVoice, processBlockVoice

  /** Stops the playback of all currently active RegionPlayers immediately. This is a rather hard 
  reset which may be appropriate to call when a midi reset message is received or before loading a
  new patch. It returns the number of players that were affected, i.e. the number of players that 
  were in active state before the call. */
  virtual int stopAllPlayers();



  /** Calls stopAllPlayers. Function is for consistency with the rest of the library. */
  void reset() { stopAllPlayers(); }

  void handleMusicalEvent(const rsMusicalEvent<float>& ev);



  //===============================================================================================

  //-----------------------------------------------------------------------------------------------
  // \name Internal Helper Classes


  class SignalProcessor
  {
  public:
    virtual void processFrame(rsFloat64x2& inOut) = 0;
    virtual void processBlock(rsFloat64x2* inOut, int N) = 0;
    virtual void resetState() = 0;
    virtual void resetSettings() = 0;
  };

  class Modulator
  {
  public:
    virtual double getSample() = 0;
    virtual void resetState() = 0;
    virtual void resetSettings() = 0;
    // todo: processBlock
  };


protected:

  class SignalProcessorChain
  {
  public:
    void processFrame(rsFloat64x2& inOut);
    void processBlock(rsFloat64x2* inOut, int N);
    void resetState();
    void resetSettings();
  protected:
    std::vector<SignalProcessor*> processors;
  };

  class ModulationConnection
  {

  private:
    std::function<void(double)> targetSetter; // target callback that sets some parameter
    double amount = 0.0;  // strength of modulation
    double refVal = 0.0;  // unmodulated reference value
  };

  /** A class for playing back a given Region object. */
  class RegionPlayer
  {

  public:

    /** Sets up the region object that this player should play. You need to also pass the output 
    sample-rate which is the sample rate at which the player should run (not the sample rate of the
    audio file associated with the region). */
    void setRegionToPlay(const Region* regionToPlay, double outputSampleRate);

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


  protected:

    /** A basic sanity check for the given region. Mostly for catching bugs. */
    bool isPlayable(const Region* region);

    /** Sets up the internal values for the playback settings (including DSP objects) according
    to the assigned region and resets all DSP objects. */
    void prepareToPlay(double sampleRate);

    bool buildProcessingChain();
    void resetDspState();
    void resetDspSettings();
    void setupDspSettings(const std::vector<PlaybackSetting>& settings, double sampleRate);

    const Region* region;                 //< The Region object that this object should play
    const AudioFileStream<float>* stream; //< Stream object to get the data from
    rsFloat64x2 amp = 1.0;                //< Amplitude (for both channels)
    //int sampleTime = 0;            //< Elapsed time in samples, negative values used for delay
    double sampleTime = 0.0;       //< Time index in the sample. Negative values used for delay.
    double increment  = 1.0;       //< Increment of sampleTime per sample
    uchar key = 0;                 //< Midi note number used for starting this player

    std::vector<Modulator*> modulators;
    std::vector<ModulationConnection*> modMatrix;  // not a literal matrix but conceptually
    SignalProcessorChain dspChain;

    // We may need a state, too. Can be attack/decay/sustain/release. Or maybe just play/release?
    // Or maybe no state at all but the triggerRelease just triggers the release of all envelopes?


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

  /** Defines a set of regions. Used to handle note-on/off events efficiently. Not to be confused 
  with groups. This class exists for purely technical reasons (i.e. implementation details) and 
  does not map to any user concept. */
  class RegionSet
  {

  public:


    int findRegionIndex(const Region* r) const
    { 
      if(regions.empty()) return -1;
      return RAPT::rsArrayTools::findIndexOf(&regions[0], r, (int)regions.size()); 
    }

    bool containsRegion(const Region* r) const {  return findRegionIndex(r) != -1; }

    void addRegion(const Region* r) 
    { 
      if(containsRegion(r)) {
        RAPT::rsError("Don't add regions twice!"); return; }
      regions.push_back(r); 
    }

    bool removeRegion(int i)
    {
      if(i < 0 || i >= (int)regions.size())  {
        RAPT::rsError("Invalid region index"); return false; }
      RAPT::rsRemove(regions, i);
      return true;
    }

    bool removeRegion(const Region* r)
    {
      return removeRegion(findRegionIndex(r));
    }

    // todo: removeRegion, containsRegion, findRegion

    int getNumRegions() const { return (int)regions.size(); }

    const Region* getRegion(int i) const { return regions[i]; }

    void clear() { regions.clear(); }

  private:

    std::vector<const Region*> regions; // pointers to the regions belonging to this set
  };


  //-----------------------------------------------------------------------------------------------
  // \name Internal functions

  /** Adds the given region to our regionsForKey array at the k-th position iff the region should
  potentially play at the key k, as determined by the region's loKey,hiKey settings. */
  void addRegionForKey(uchar k, const Region* region);
  // maybe rename to addRegionForNoteOn and have a similar functionality for noteOff

  /** Returns true, iff the given region should play when the given key is pressed with given 
  velocity. This will also take into account other playback constraints defined for the region 
  and/or its enclosing group. */
  bool shouldRegionPlay(const Region* r, uchar key, uchar vel);

  /** Finds the group index and region index within the group for the given region and assigns the
  output parameters group/regionIndex accordingly. If the region is not found in any of our groups,
  both will be assigned to -1. This should actually not happen, though: a region should always be 
  found in one and only one group. In sfz files, there can actually be regions that are not in any
  group, but in this implementation, we just put them into an additional group (at group index 0), 
  if necessary. Handling it all uniformly is more elegant and efficient and the organization of the
  sfz file does not need to match all of our implementation details here exactly. The group with 
  index 0 contains either all the free regions or corresponds to the first defined group in the sfz
  file, if the file contains no free regions. If (-1,-1) is returned, it indicates that the caller 
  still holds a pointer to a region that doesn't exist anymore in the instrument, which may 
  indicate a bug at the call site. */
  void findRegion(const Region* region, int* groupIndex, int* regionIndex);

  /** Returns a pointer to a player for the given region by grabbing it from the idlePlayers array.
  This will have the side effects that this player will be removed from idlePlayers and added to 
  activePlayers. If no player is available (i.e. idle), this will return a nullptr. The caller 
  should interpret that as a layerOverload condition and return the appropriate return code to 
  its caller. */
  RegionPlayer* getRegionPlayerFor(const Region* r, uchar key, uchar vel);

  /** Returns true, iff the given sample is used in the instrument definition represented by the 
  given sfz */
  bool isSampleUsedIn(const AudioFileStream<float>* sample, const rsSamplerData& sfz);
  
  /** Stops the player at the given "activeIndex" which is the index into our "activePlayers" 
  array. This results in the removal of the player from "activePlayers" and adding it back to
  "idlePlayers". The return value is either rsReturnCode::success or rsReturnCode::invalidIndex, if
  the activeIndex was not a valid index into our activePlayers array. */
  int stopRegionPlayer(int activeIndex);

  /** Returns the AudioFileStream object that is used to stream the actual sample data for the
  given region. A pointer to this object is supposed to be stored within the region object
  ...tbc... */
  static const AudioFileStream<float>* getSampleStreamFor(const Region* r);

  /** Handles a noteOn event with given key and velocity and returns either rsReturnCode::success, if
  we had enough voices available to serve the request or rsReturnCode::voiceOverload, in case the 
  noteOn could not be handled due to inavailability of a sufficient number of idle voices. If no
  sufficient number of idle voices was available and the noteOn should actually have triggered 
  playback of multiple samples, none of them will be triggered. It's an all-or-nothing thing: we 
  don't ever trigger playback for only a subset of samples for a given noteOn. */
  virtual int handleNoteOn(uchar key, uchar vel);

  /** Analogous to handleNoteOn. It may also return rsReturnCode::voiceOverload in cases where the 
  noteOff is supposed to trigger relase-samples. In such a case, none of the release-samples will 
  be triggered. */
  virtual int handleNoteOff(uchar key, uchar vel);


  /** Removes those samples from our sample pool that are not used in the given sfz instrument 
  specification. Returns the number of samples that were removed. */
  int removeSamplesNotUsedIn(const rsSamplerData& sfz);
  // maybe rename to removeUnusedSamples. But that name is more ambiguous: it could be interpreted
  // as "unused in the current sfz member", so maybe don't

  /** Adds all samples to our sample pool that are used in the given sfz instrument definition, if
  they are not already there. Returns the number of samples that were added or 
  rsReturnCode::fileLoadError if any of the files failed to load. */
  int addSamplesUsedIn(const rsSamplerData& sfz);
  // maybe rename to loadSamples or loadSamplesFor

  /** Sets up all the AudioStream pointers in all the regions in our sfz member. */
  int setupAudioStreams();

  /** Sets up ourregionsForKey array according to our sfz member. */
  void setupRegionsForKey();


  //-----------------------------------------------------------------------------------------------
  // \name Data

  rsSamplerData sfz;
  /**< The data structure that defines the sfz instrument. */
 
  static const int numKeys = 128;
  RegionSet regionsForKey[numKeys];
  /**< For each key, we store a set of regions that *may* need to be played, when the key is 
  pressed. Whether or not a region is a candidate for playback for a given key is determined by the
  loKey, hiKey settings of that region. If the playback candidate region then *really* needs to be 
  played in a particular situation is determined by other constraints as well, such as velocity 
  range, last received controller and/or pitch-wheel values, etc. The key is the first and primary 
  filter for which regions need to be played when a noteOn is received and the purpose of this 
  array is to optimize this primary filter to avoid having to loop through all regions in the 
  instrument on each received noteOn. Secondary, tertiary, etc. filters may follow and are 
  implemented by indeed looping through all candidate regions for a given key. It is assumed that 
  the number of candidate regions for each key is typically much smaller than the total number of 
  regions in the instrument - like a few instead of a few hundred. */
  // maybe use a std::vector, maybe we need a similar array for note-off samples

  SamplePool<float> samplePool;
  /**< The pool of samples that are in use for the currently loaded instrument. The samples are 
  pooled to avoid redundant storage in memory when multiple regions use the same sample. */
  // Maybe the pool should be a pointer, so it can be shared between multiple instances and/or with
  // some other object that also uses the same sample pool (for example, a DAW that uses the sampler
  // as plugin)


  //std::vector<PlaybackSetting> settings;
  /**< Playback settings that apply to all groups within this instrument, unless a group (or 
  region) overrides a setting with its own value. **/
  // get rid - should go into sfz.instrument.settings

  std::vector<RegionPlayer*> activePlayers;
  /**< Array of pointers to region players that are currently active, i.e. playing. */
  // rename to activeRegionPlayers

  std::vector<RegionPlayer*> idlePlayers;
  /**< Array of pointers to region players that are currently idle, i.e. not playing, and therefore
  available to be used for new incoming notes. */

  std::vector<RegionPlayer> playerPool;
  /**< This is our pool of the actual RegionPlayer objects from which we grab a player when we need
  to trigger the playback of a region. The invariant is that at any given time, all players in this 
  pool are either pointed to by some element in the activePlayers or by some element in the 
  idlePlayers. The size of this array determines our maximum number of layers and the memory usage 
  without taking memory for the samples into account. */
  // todo: use a sort of multi-threaded, speculatively pre-allocating, non-deallocating dynamic 
  // array data structure for that later. The same strategy should then later be used for DSP 
  // objects as well


  double sampleRate = 44100.0;
  /**< Sample rate at which this object runs. */

  // Some info that can be inquired from client code after loading a new sfz file:
  int numSamplesRemoved = 0;  /**< Number of samples that were unloaded. */
  int numSamplesLoaded  = 0;  /**< Number of samples that were loaded. */
  int numSamplesFailed  = 0;  /**< Number of samples that failed to load. */

  // These are not yet used - currently, both are assumed to be the project directory (at least
  // for the unit tests):
  std::string sfzDir;         /**< Root directory for .sfz files */
  std::string wavDir;         /**< Root directory for .wav files */


  //int numChannels = 2;
  /**< The number of output channels. By default, we have two channels, i.e. a stereo output. */
  // maybe that should be determined by TSig? multi-channel output should be realized by using
  // a multichannel (simd) type ...maybe get rid and support only stereo output


  //float midiCC[128];     // most recently received controller values in 0...127
  //float midiPitchWheel;  // most recently received pitch-wheel value in -8192...+8291

  // Maybe have buffers for the outputs
  // TSig* outBuffer;
  // maybe for the pitch-env, we can render the pitch-env itself into the buffer first and then 
  // overwrite it with the actual output data

  // we may need a pool of filter objects, eq-objects, etc. for the different voices
};

//=================================================================================================

/** A subclass of rsSamplerEngine that adds a couple of features. In particular, it's meant for
adding those features that are not part of the original SFZ specification. These include more 
flexible signal routing capabilities like the ability to apply the group- and instrument-wide 
settings on top of the region settings instead of using them as fallback values. */

class rsSamplerEngine2 : public rsSamplerEngine
{

public:

  // for convenience:
  using uchar = unsigned char;
  //using Region = rsSamplerData::Region; // todo: make a subclass here that adds the stream field
  //using Group  = rsSamplerData::Group;
  //using PlaybackSetting = rsSamplerData::PlaybackSetting;


  rsSamplerEngine2(int maxNumLayers);


  void setMaxNumLayers(int newMax) override;


  /** Decides if the group settings should be applied on top of the region settings (true) or if 
  they should just act as fallback values for when a region doesn't define them (false). The 
  "on top" mode means that if a region defines a gain of -6dB and its enclosing group defines a 
  gain of -3dB, the total gain will be -9dB. The "fallback" mode means that the region will just 
  use it's defined -6dB gain and only if that would not be defined, it would fall back to the 
  group's -3dB setting. The latter behavior is the default in sfz (verify!) but the the former is 
  also often convenient. 
  This doesn't do anything yet...this feature is not yet implemented - for the time being, it's 
  just the infrastructure. */
  void setGroupSettingsOnTop(bool onTop) { groupSettingsOnTop = onTop; }

  /** Decides if the group modulations should be applied on top of the region modulations (true) or
  if their settings should just act as fallback values for when a region doesn't define them 
  (false). Note that technically, we don't have a set of modulators per group. Instead, if "on top"
  mode is selected all RegionPlayers must actually duplicate all their modulators, one using the
  region's own settings and one using the enclosing group's settings and add them up before 
  applying. It doesn't seem to make a lot of sense to run common modulators per-group. Think of an
  envelope: it gets triggered with the note that starts the RegionPlayer, so it must be part of the
  RegionPlayer. ...or actually, it could make sense to trigger the modulators with the first 
  RegionPlayer for each group - this is not a useful behavior for envelopes, but it could be for 
  LFOs, sequencers, etc. ...maybe the modulators should have 3 modes: override/fallback, 
  accumulate/duplicate, accumulate...actually, it would be useful, if this could be set for each 
  modulator individually. */
  void setGroupModulationsOnTop(bool onTop) { groupModulationsOnTop = onTop; }

  /** Like setGroupSettingsOnTop, but for the instrument settings. */
  void setInstrumentSettingsOnTop(bool onTop) { instrumentSettingsOnTop = onTop; }

  /** Like setGroupModulationsOnTop, but for the instrument modulations. */
  void setInstrumentModulationsOnTop(bool onTop) { instrumentModulationsOnTop = onTop; }




  void processFrame(double* left, double* right) override;

  //void processFrame(float* left, float* right) override;

  //void processBlock(float** block, int numFrames) override;

  // void processFrameVoice, processBlockVoice

  /** Stops the playback of all currently active RegionPlayers immediately. This is a rather hard 
  reset which may be appropriate to call when a midi reset message is received or before loading a
  new patch. It returns the number of players that were affected, i.e. the number of players that 
  were in active state before the call. */
  int stopAllPlayers() override;


protected:


  int handleNoteOn(uchar key, uchar vel) override;

  int handleNoteOff(uchar key, uchar vel) override;



  /** A class for collecting all the SignalProcessors that apply to a given group. This is used 
  only when the group's DSP settings should go on top of the region's settings */
  class GroupPlayer  // maybe rename to GroupPlayer
  {

  public:

    /** Generates one stereo sample frame at a time. */
    rsFloat64x2 getFrame();

  protected:

    std::vector<RegionPlayer*> regionPlayers;
    SignalProcessorChain dspChain;

    rsSamplerEngine2* engine = nullptr;
    // We need a communication channel to the enclosing sampler-engine.

  };


  // Flags to decide if the group- and/or instrument settings and/or modulations should be applied 
  // on top of the region settings/modulations:
  bool groupSettingsOnTop         = false;
  bool instrumentSettingsOnTop    = false;
  bool groupModulationsOnTop      = false;
  bool instrumentModulationsOnTop = false;
  // maybe move this into a subclass rsSamplerEngine2 - this is an added non-sfz feature


  // under construction:
  std::vector<GroupPlayer*> activeGroupPlayers;
  std::vector<GroupPlayer*> idleGroupPlayers;
  std::vector<GroupPlayer>  groupPlayerPool;
  // I think, we need as many GroupPlayers as there are RegionPlayers because in the worst case,
  // each region could be within its own group...although that's probably really uncommon

};

//=================================================================================================

/** Subclass that contains some extra functions that facilitate testing which should not go into 
the production code. */
class rsSamplerEngineTest : public rsSamplerEngine
{

public:

  using rsSamplerEngine::rsSamplerEngine;  // inherit constructors

  /** Returns true, iff this object is in the same state (with regard to content of the sample pool
  and instrument definition) as the given other engine */
  bool isInSameStateAs(const rsSamplerEngineTest& other) const
  {
    return sfz == other.sfz; // && samplePool.hasSameContentAs(other.samplePool);
    // todo: compare also content of samplePool in both objects..the function may optionally allow
    // a different order of the samples in both pools
    // maybe compare also the regionsForKey arrays
  }
  // maybe rename to hasSameInstrument - isInSameState may also compare the state with regard to 
  // activeLayers, etc.

  /** Returns the byte size of the RegionPlayer class. We want to keep this small so we use this 
  function to keep track of its size in the tests. */
  static int getRegionPlayerSize() { return sizeof(rsSamplerEngine::RegionPlayer); }

};
// maybe move into the test project or the rs_testing juce module



}
#endif