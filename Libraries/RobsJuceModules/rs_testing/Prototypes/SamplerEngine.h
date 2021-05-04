#ifndef RAPT_SAMPLERENGINE_H
#define RAPT_SAMPLERENGINE_H

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

protected:

  Type type;  // e.g. noteOn/Off, controlChange, pitchWheel
  T    val1;  // e.g. key, controller number, pitchWheelMSB
  T    val2;  // e.g. velocity, controller value, pitchWheelLSB

};

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

  rsSamplerEngine(int maxPolyphony = 16);

  virtual ~rsSamplerEngine();


  //-----------------------------------------------------------------------------------------------
  // \name Helper Classes

  using uchar = unsigned char;
  class Region;
  class AudioFileStream;

  /** A class to represent various additional (and optional) playback settings of a region, group 
  or instrument. Such additional settings include additional constraints for the circumstances 
  under which a particular sample should be played. Key- and velocity ranges are the obvious 
  primary constraints (and they therefore are directly baked into the Region class below), but sfz
  defines many more. But the settings doesn't need to be playback constraints - that's only one 
  type of setting. Other types are things like envelope settings, filter frequencies, etc. */
  class PlaybackSetting
  {

  public:

    /** Enumeration of possible types of settings. These types correspond to the opcodes defined
    in the sfz specification. */
    enum Type
    {
      ControllerRangeLo, ControllerRangeHi, PitchWheelRange,  // 

      PitchKeyCenter,

      AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,

      FilterCutoff, FilterResonance, FilterType,

      Unknown,
      NumTypes

      //...tbc...
    };

    PlaybackSetting(Type type, float value)
    { this->type = type; this->value = value; }

    Type getType() const { return type; }

    /** Returns the stored value. Values are always stored as floats and it is understood that in
    cases, where the corresponding parameter in the sfz spec is defined to be an integer, we just
    represent it as that integer with all zeros after the decimal dot. */
    float getValue() const { return value; }
    // todo: (decide and) document, how choice parameters like filter-type are represented. In sfz,
    // they are just represented as text. Maybe for each such choice parameter, we need another 
    // enum to represent its allowed values...

  private:

    Type  type  = Type::Unknown;
    float value = 0.f;
    // hmm - it seems, for the controllers, we need 2 values: controller number and value - but it
    // would be wasteful to store two values for all other settings as well...hmmm...maybe 
    // groups/regions need to maintain 2 arrays with settings, 1 for the 1-valued settings and 
    // another for the 2-valued settings - maybe have classes PlaybackSetting, PlaybackSetting2Val
    // or: have indeed all the ccN as different type - but that would blow up the enum excessively
  };

  /** A group organizes a bunch of regions into a single entity for which performance settings can 
  be set up which will be applicable in cases where the region does not itself define these 
  settings, so they act as fallback values. */
  class Group
  {

  public:


    ~Group() { clearRegions(); }

    int addRegion();   // make private!

    // todo: removeRegion, etc.

    /** Returns the index of the given region within this group, if present or -1 if the region is
    not present in this group. */
    int getRegionIndex(const Region* region); // make const


    /** Returns true, if the given index i refers toa valid region within this group. */
    bool isRegionIndexValid(int i) const { return i >= 0 && i < (int)regions.size(); }




    const Region* getRegion(int i) const
    {
      if(i < 0 || i >= (int)regions.size()) {
        rsError("Invalid region index");
        return nullptr; 
      }
      //return &regions[i];
      return regions[i];
    }
    // for some reason, i get compiler errors when trying to put this into the cpp file 
    // -> figure out

    Region* getRegionNonConst(int i)
    {
      if(i < 0 || i >= (int)regions.size()) {
        rsError("Invalid region index");
        return nullptr; 
      }
      //return &regions[i];
      return regions[i];
    }
    // move to cpp


    void clearRegions();   // make private
    void clearSettings();  // ditto

  private:

    std::vector<Region*> regions;
    /**< Pointers to the regions belonging to this group. */

    std::vector<PlaybackSetting> settings;
    /**< Settings that apply to all regions within this group, unless a region overrides them with
    its own value for a particular setting. */

    // may be add these later:
    //std::string name;  
  };

  /** A region contains a sample along with performance settings including information for which 
  keys and velocities the sample should be played and optionally other constraints for when the the
  sample should be played and also settings for pitch, volume, pan, filter, envelopes, etc. */
  class Region
  {

  public:

    /** Returns a (const) pointer the audio stream object that should be used for this region. */
    const AudioFileStream* getSampleStream() const { return sampleStream; }

    /** Returns a const reference to our playback settings. */
    const std::vector<PlaybackSetting>& getSettings() const { return settings; }

  private:


    /** Sets the audio stream object that should be used for this region. */
    void setSampleStream(const AudioFileStream* newStream) { sampleStream = newStream; }

    const AudioFileStream* sampleStream = nullptr;  
    Group* group = nullptr;             // pointer to the group to which this region belongs
    uchar loKey = 0, hiKey = 127;
    uchar loVel = 0, hiVel = 127;
    // todo: maybe package loKey/hiKey, loVel/hiVel into a single uchar to save memory

    std::vector<PlaybackSetting> settings;
    // for more restrictions (optional) restrictions - sfz can restrict the playback of samples
    // also based on other state variables such as the last received controller of some number,
    // last received pitchwheel, etc. ...but maybe a subclass RestrictedRegion should be used
    // for that - i don't think, it will be used a lot and will just eat up memory when it's
    // present in the baseclass...or maybe it should have a more general array of RegionFeatures
    // which may also include loop-settings and the like

    //std::string name;

    friend class Group;
    friend class rsSamplerEngine;
    // The Region class shall not provide any public functions that can modify the region because
    // those could be used by client code to modify the region behind the back of the 
    // rsSamplerEngine which could mess things up. Client code can modify regions only through the
    // appropriate functions of rsSamplerEngine. It acts as man-in-the-middle and can the call the
    // private setters of the Region (by virtue of being a friend class) and it may also trigger 
    // additional actions, if necessary. The same should probably apply to the Group class as well.
    // Is this a known pattern? -> figure out
  };


  /** A class to represent a pool of audio samples...tbc...
  ToDo:
  Maybe factor out, i.e. move out of rsSamplerEngine. It may be useful in other contexts as well. 
  The same goes for the AudioStream class and its subclasses. But in a more general context, we 
  would probably need a more complex referencing system for the AudioStream objects - regions
  would need to be something like AudioStreamClients, that register/deregister themselves, etc.
  Maybe AudioStream should be a subclass of some DataStream baseclass. We'll see... */
  class SamplePool
  {

  public:


    ~SamplePool() { clear();  }


    int addSample(const AudioFileStream* newSample)
    {
      // rsAssert(!contains(newSample))
      samples.push_back(newSample);
      return ((int) samples.size()) - 1;
    }
    // should the pool take ownership? ...i think so


    /** Returns true, if the given index i refers to a valid sample index. */
    bool isSampleIndexValid(int i) const { return i >= 0 && i < (int)samples.size(); }


    const AudioFileStream* getSampleStream(int i)
    {
      if(!isSampleIndexValid(i)) {
        rsError("Invalid sample index");
        return nullptr; }
      return samples[i];
    }

    void clear();


    // todo:
    // setup: removeSample...but if regions refer to it, we need to update them, too by 
    // invalidating their pointers. We either need an observer mechanism (complex and general) or 
    // we allow reomval of samples only via a member function of the outlying rsSamplerEngine 
    // class, which also takes care of resetting the sample-streams in all regions that use it
    // (simpler but less general). Maybe to make it safer, we could also introduce a reference
    // counter and check, if it is zero, before a stream objects gets removed
    // inquiry: hasSample


  protected:

    std::vector<const AudioFileStream*> samples;

  };


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Return codes for the setup functions. We use encodings as negative integers so we can use 
  them also for functions which use positive integers as valid return values. */
  enum ReturnCode
  {
    success       = -1,  //< Operation completed successfully. 
    nothingToDo   = -2,  //< There was nothing to actually do. State was already as desired.
    memAllocFail  = -3,  //< Memory allocation failure.
    invalidIndex  = -4,  //< An invalid index was passed.
    voiceOverload = -5,  //< Not enough free voices available (in e.g. new noteOn).
    notFound      = -6   //< A region, group or whatever was not found.
  };
  // todo: make it an enum class, maybe include also return codes for inquiry functions such as for
  // "unknown", etc. ...but maybe that's no good idea when we want to use it for functions which
  // need to return valid integers (like, for numChannels, etc. - we could use negative numbers to
  // encode such things)
  // maybe rename "success" to "completed" because "success" has actually a more general meaning:
  // "nothingToDo" is also a kind of "success" (or maybe "workDone" or "workCompleted"

  /** Adds a new sample to our pool of samples. After the sample has been added, regions can be 
  defined that make use of it. */
  int addSampleToPool(float** data, int numFrames, int numChannels, float sampleRate, 
    const std::string& uniqueName);
  // Maybe rename to addSample, it should return the index of the sample in the sample-pool
  // maybe make a struct SampleMetaData containing: numFrames, numChannels, sampleRate, rootKey
  // todo: take reference to a metaData object

  /** Adds a new group to the instrument definition and returns the index of the group. */
  int addGroup();

  /** Adds a new region to the group with the given index and returns the index of the region 
  within the group or ReturnCode::invalidIndex, if the passed groupIndex was invalid. If the key 
  range is already known, it makes sense to pass it using the optional loKey/hiKey parameters. This
  can also be set up later, but some memory operations can be saved, if it's known in advance. */
  int addRegion(int groupIndex, uchar loKey = 0, uchar hiKey = 127);

  /** Sets the sample to be used for the given region within the given group. Returns either
  ReturnCode::success or ReturnCode::invalidIndex, if the pair of group/region indices and/or the
  sample index was invalid. */
  int setRegionSample(int groupIndex, int regionIndex, int sampleIndex); 

  /** Sets a value for a given type of playback setting for the given region. Return either
  ReturnCode::success
  */
  int setRegionSetting(const Region* region, PlaybackSetting::Type type, float value);

  // todo: setGroupSetting, setInstrumentSetting

  //int setRegionSetting(Region

  // todo:
  // setRegionLoKey, setRegionHiKey, setRegionLoVel, setRegionHiVel


  // todo: removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns a pointer to the (const) region object with the given group- and region index or a 
  nullptr if the combination of indices is invalid. */
  const Region* getRegion(int groupIndex, int regionIndex)
  {
    int gi = groupIndex, ri = regionIndex;
    if(gi < 0 || gi >= (int)groups.size()) {
      rsError("Invalid group index");
      return nullptr; 
    }
    return groups[gi].getRegion(ri);
  }
  // move to cpp
  // Maybe it should be non-const - but no: the caller should not be able to change the loKey/hiKey
  // settings because that would require a change to the regionsForKey array...but maybe we should 
  // get rid of that anyway - it might be a pointless attempt to optimization -> benchmark!
  // for some reason, i get compiler errors when trying to put this into the cpp file 
  // -> figure out

  // getGroup, getRegion, getStateAsSFZ, isSampleInPool, getNumGroups, getNumRegionsInGroup(int)
  // getNumRegions(), getSampleIndex(const string& uniqueName) ..or maybe it should take a pointer
  // to a SampleMetaData object

  /** Returns true, iff the given pair of group- and region index is valid, i.e. a region with this
  pair of indices actually exists in the current instrument definition. */
  bool isIndexPairValid(int groupIndex, int regionIndex) const
  {
    int gi = groupIndex, ri = regionIndex;
    return gi >= 0 && gi < (int)groups.size() && groups[gi].isRegionIndexValid(ri);
  }

  /** Returns true, iff the given sample index is valid, i.e. a sample with this index actually 
  exists our sample pool. */
  bool isSampleIndexValid(int sampleIndex) const
  { return samplePool.isSampleIndexValid(sampleIndex); }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(float* frame);
  // maybe have frameL, frameR inputs

  void processBlock(float** block, int numFrames);

  void handleMusicalEvent(const rsMusicalEvent<float>& ev);

  // void processFrameVoice, processBlockVoice


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internal Helper Classes

  class AudioStream
  {

  public:

    virtual ~AudioStream() {}

    /** For random access. Writes the sample frame with given index into the given destination. */
    virtual void getFrame(int sampleIndex, float* destination) = 0;

    /** Subclasses may want to override this for optimizing the blockwise access. */
    /*
    virtual void getBlock(int startIndex, int length, TSmp** destination)
    {
      for(int i = 0; i < length; i++)
        getFrame(startIndex + i, destination[i]);
        // assumes interleaved storage in memory
    }
    */

  protected:

    float sampleRate  = 44100.f;
    int   numChannels = 0;
    int   numFrames   = 0;         // maybe use -1 to encode "unknown"? would that be useful?

  };
  // maybe move out of this class - this may be useful in other contexts, too - maybe templatize

  class AudioFileStream : public AudioStream
  {

    // ToDo: add == operator based on fileName, extension, path (and maybe (meta)data such as 
    // sampleRate, numChannels, numFrames, ...)

    // inquiry: existsOnDisk (idea: we should allow to work with samples that are not stored on 
    // disk but rather pre-rendered programmatically into RAM at runtime)

    //virtual void getFrame(int sampleIndex, TSmp* destination) {}

    // todo: getBlock


  protected:

    std::string fileName;  // without filename extension (e.g. Piano_A4)
    std::string extension; // filename extension (e.g. wav, flac)
    std::string path;      // relative path from a predefined root directory
    int rootDirIndex = 0;  // index of root directory (among a couple of predefined choices)
    // rootDirIndex stores the root-directory to which the path is relative, but just as an integer
    // index that selects between various pre-defined root-directories that should exist at the 
    // rsSamplerEngine level. For example: 0: instrument directory (containing e.g. Piano.sfz), 
    // 1: factory sample directory, 2: user sample directory, etc. (but not hardcoded - meanings of
    // the indices should be flexible). This makes it more reasonably possible to uniquely identify 
    // samples. It's totally possible to have samples in an instrument with same relative paths and 
    // filenames but with respect to different root directories. Yes - that would be weird, but the 
    // engine should neverless be able to handle such situations.
  };

  // maybe rename to AudioFileStreamRAM, another subclass can be named AudioFileStreamDFD
  class AudioFileStreamPreloaded : public AudioFileStream 
  {

  public:


    virtual ~AudioFileStreamPreloaded() { clear(); }


    int setData(float** newData, int numFrames, int numChannels, float sampleRate, 
      const std::string& uniqueName);
    // todo: include fileName etc. later, too


    void clear();


    void getFrame(int sampleIndex, float* destination) override
    {
      int n = sampleIndex;
      rsAssert(n >= 0 && n < numFrames, "sampleIndex out of range");
      for(int c = 0; c < this->numChannels; c++)
        destination[c] = channelPointers[c][n];
        // What, if the number of output channels shall be different than the number of of channels
        // in the data? Maybe it's best (most efficient) to ensure that this cannot happen on a 
        // higher level. We can just let our channelPointers all point to the same actual data, for
        // example. When the number of output channels of the sampler engine is changed (which 
        // should happen rarely, if ever), we need to update all channelPointer arrays in all 
        // AudioFileStreamPreloaded objects.
    }

  protected:

    float*  flatData = nullptr;         // pointer to the sample data
    float** channelPointers = nullptr;  // pointers to the channels
    // If we store the data in interleaved format, the channelPointers will be not needed and 
    // getFrame must be implemented differently. Maybe that's better (more efficient)

  };



  /** A class for playing back a given Region object. */
  class RegionPlayer
  {

  public:

    /** Sets up the region object that this player should play. */
    virtual void setRegionToPlay(const Region* regionToPlay);

    /** Generates one stereo sample frame at a time. */
    virtual rsFloat64x2 getFrame();

    /** Writes a block of given length into the outBuffer. */
    virtual void processBlock(rsFloat64x2* outBuffer, int length);

  protected:

    /** Sets up the internal values for the playback settings (including DSP objects) according
    to the assigned region and resets all DSP objects. */
    virtual void prepareToPlay();

    using Biquad = RAPT::rsBiquadDF1<rsFloat64x2, double>; // todo: use TDF2

    const Region* region;          //< The Region object that this object should play
    const AudioFileStream* stream; //< Stream object to get the data from
    rsFloat64x2 amp = 1.0;         //< Amplitude (for both channels)
    int sampleTime = 0;            //< Elapsed time in samples, negative values used for delay
    Biquad flt, eq1, eq2, eq3;     //< Filter and equalizers

    // ToDo: 
    // -add more DSP objects: envelope generators, LFOs, we also need to somehow take care
    //  of the effect sends (reverb, chorus)
    // -Maybe in addition to elapsed time, we also need to maintain position inside the sample 
    //  stream. This is not the same thing due to possible delay and looping.
    // -maybe use a different implementation structure (SVF) for the time-varying filter
    // -try to optimize ram and/or cpu usage by re-ordering
    // -env-generators need as state variables: stage (int), time-into-stage (float/double),
    //  LFO need just phase
  };

  /** Defines a set of regions. Used to handle note-on/off events efficiently. Not to be confused 
  with groups. This class exists for purely technical reasons (i.e. implementation details) and 
  does not map to any user concept. */
  class RegionSet
  {

  public:

    void addRegion(const Region* r) { regions.push_back(r); }
    // todo: removeRegion, containsRegion

    size_t getNumRegions() const { return regions.size(); }

    const Region* getRegion(size_t i) const { return regions[i]; }

  private:

    std::vector<const Region*> regions; // pointers to the regions belonging to this set
  };


  //-----------------------------------------------------------------------------------------------
  // \name Internal functions


  /** Adds the given region to our regionsForKey array at the k-th position. */
  void addRegionForKey(uchar k, const Region* region);

  /** Returns true, iff the given region should play when the given key is pressed with given 
  velocity. This will also take into account other playback constraints defined for the region 
  and/or its enclosing group. */
  bool shouldRegionPlay(const Region* r, const char key, const char vel);

  /** Returns a non-constant point to the region with given index pair, so this allow modifying the
  Region object via the pointer. Use it with care. In particular, don't change the loKey, hiKey 
  settings because such changes require an update of the regionsForKey array. */
  Region* getRegionNonConst(int groupIndex, int regionIndex)
  {
    int gi = groupIndex, ri = regionIndex;
    if(gi < 0 || gi >= (int)groups.size()) {
      rsError("Invalid group index");
      return nullptr; 
    }
    return groups[gi].getRegionNonConst(ri);
  }
  // Move to cpp file or better: get rid and keep only getRegion which shall return a non-const 
  // pointer. Client code actually should be able to modify the region, but only indirectly, 
  // mediated through the sampler engine.


  /** Finds the group index and region index within the group for the given region and assigns the
  output parameters group/regionIndex accordingly. If the region is not found in any of our groups,
  both will be assigned to -1. This should actually not happen, though: a region should always be 
  found in one and only one group. In sfz files, there can actually be regions that are not in any
  group, but in this implementation, we just put them into an additional group (at group index 0), 
  if necessary. The organization of the sfz file does not need to match the implementation. The 
  group with index 0 contains either all the free regions or corresponds to the first defined group 
  in the sfz file, if the file contains no free regions. If -1,-1 is returned, it indicates that 
  the caller still holds a pointer to a region that doesn't exist anymore in the instrument, which 
  may indicate a bug at the call site. */
  void findRegion(const Region* region, int* groupIndex, int* regionIndex);


  /** Returns a pointer to a player for the given region by grabbing it from the idlePlayers array.
  This will have the side effects that this player will be removed from idlePlayers and added to 
  activePlayers. If no player is available (i.e. idle), this will return a nullptr. */
  RegionPlayer* getRegionPlayerFor(const Region* r);

  /** Handles a noteOn event with given key and velocity and returns the number of voices that
  were triggred or ReturnCode::voiceOverload, in case the noteOn could not be handled due to 
  inavailability of a sufficient number of idle voices. */
  int handleNoteOn(uchar key, uchar vel);

  int handleNoteOff(uchar key, uchar vel);

  // return code should inform, whether a region was triggered (maybe how many) and should also
  // return an error code for when not enough idle voices are available



  //-----------------------------------------------------------------------------------------------
  // \name Data
 
  RegionSet regionsForKey[128];
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
  // maybe use a std::vector

  SamplePool samplePool;
  /**< The pool of samples that are in use for the currently loaded instrument. The samples are 
  pooled to avoid redundant storage in memory when multiple regions use the same sample. */

  std::vector<Group> groups;
  /**< The groups contained in this instrument. Each group may contain set of regions. */

  std::vector<PlaybackSetting> settings;
  /**< Playback settings that apply to all groups within this instrument, unless a group (or 
  region) overrides a setting with its own value. **/

  std::vector<RegionPlayer*> activePlayers;
  /**< Array of pointers to region players that are currently active, i.e. playing. */

  std::vector<RegionPlayer*> idlePlayers;
  /**< Array of pointers to region players that are currently idle, i.e. not playing, and therefore
  available to be used for new incoming notes. */

  std::vector<RegionPlayer> playerPool;
  /**< This is our pool of the actual RegionPlayer objects for which we grab a player when we need
  to trigger the playback of a region. The invariant is that at any given time, all players in this 
  pool are either pointed to by some element in the activePlayers or by some element in the 
  idlePlayers. The size of this array determines our maximum polyphony and the memory usage 
  without taking memory for the samples into account. */


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




/** Subclass that contains some extra functions that facilitate testing which should not go into 
the production code. */
class rsSamplerEngineTest : public rsSamplerEngine
{

public:

  static int getRegionPlayerSize() { return sizeof(rsSamplerEngine::RegionPlayer); }

};


#endif