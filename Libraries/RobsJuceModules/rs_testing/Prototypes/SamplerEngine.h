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

//=================================================================================================

class AudioStream
{

public:

  virtual ~AudioStream() {}

  /** For random access. Writes the sample frame with given index into the given destination. */
  virtual void getFrame(int sampleIndex, float* destination) const = 0;

  /** Function for speficically handling stereo signals to allow handling that importnat, common 
  special case more efficiently than with the more general implementation. */
  virtual void getFrameStereo(int sampleIndex, float* left, float* right) const = 0;

  /** Sets the number of output channels for this object. By default, this number will be equal to
  the number of channels in the data, as set by the setData call. However, the number of desired 
  output channels to be produced may actually be different from that. For example, we may want to
  produce stereo output from mono data by just copying the single output into both channels. Let
  M be the number of channels in the data and N be the number of output channels of the stream, 
  then we will produce: out[n] = data[m % N] where n is the output channel index and m % N the data
  channel index. */
  //virtual void setNumOutputChannels(int newNumChannels) = 0;



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


// maybe rename to AudioFileStreamRAM, another subclass can be named AudioFileStreamDFD,
// or ...StreamFromMemory/FromDisk
class AudioFileStreamPreloaded : public AudioFileStream 
{

public:


  virtual ~AudioFileStreamPreloaded() { clear(); }


  /** Returns true, iff everything went alright and false if it failed to allocate the required
  memory. */
  bool setData(float** newData, int numFrames, int numDataChannels, float sampleRate, 
    int numStreamChannels, const std::string& uniqueName);
  // todo: 
  // -we need to distiguish between the number of channels in the data and the desired number of 
  //  output channels
  // -include fileName etc. later, too

  //void setNumOutputChannels(int newNumChannels) override;


  void clear();


  void getFrame(int sampleIndex, float* destination) const override
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
  // this api sucks! pass a float** destinations - but for the sfz player, this is too general
  // we really need to support only mono and stereo. Maybe make a function getFrameStereo


  void getFrameStereo(int sampleIndex, float* left, float* right) const override
  {
    rsAssert(numChannels == 2); // Can be used only for stereo signals
    int n = sampleIndex;
    *left  = channelPointers[0][n];
    *right = channelPointers[1][n];
    //*left  = 0.0;
    //*right = 0.0;
  }

  //void getFrameStereo(int sampleIndex, float* destination) const




protected:

  float*  flatData = nullptr;         // pointer to the sample data
  float** channelPointers = nullptr;  // pointers to the channels
  // If we store the data in interleaved format, the channelPointers will be not needed and 
  // getFrame must be implemented differently. Maybe that's better (more efficient)
};

//=================================================================================================

/** A class to represent a pool of audio samples...tbc...
ToDo: in a more general context, we would probably need a more complex referencing system for the 
AudioStream objects - regions would need to be something like AudioStreamClients, that 
register/deregister themselves, etc. Maybe AudioStream should be a subclass of some DataStream 
baseclass. We'll see... */
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
// maybe templatize, rename to AudioFileStreamPool


//=================================================================================================

/** Data structure to define sample based instruments conforming to the sfz specification. */

class rsDataSFZ // todo: move into its own pair of .h/.cpp files
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Helper Classes

  using uchar = unsigned char;
  class Region;
  class Group;
  class Instrument;


  //-----------------------------------------------------------------------------------------------
  /** A class to represent various playback settings of a region, group or instrument. Such 
  settings include constraints for the circumstances under which a particular sample should be 
  played. Key- and velocity ranges are the obvious primary constraints (and they therefore are 
  directly baked into the Region class below), but sfz defines many more. But settings don't 
  need to be playback constraints. Other types are things like envelope settings, filter 
  frequencies, etc. */
  class PlaybackSetting  // rename to Setting
  {

  public:

    /** Enumeration of possible types of settings. These types correspond to the opcodes defined
    in the sfz specification. */
    enum Type
    {
      ControllerRangeLo, ControllerRangeHi, PitchWheelRange,  // 

      PitchKeyCenter,

      Volume,
      AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,

      FilterCutoff, FilterResonance, FilterType,

      Unknown,
      NumTypes

      //...tbc...
    };
    // maybe don't capitalize first letter - make it conistent with other (newer) enums in the 
    // library

    enum FilterType
    {
      Off, Lowpass_6, Lowpass_12, Highpass_6, Highpass_12, Bandpass_6_6, Bandreject_6_6,
      
      NumFilterTypes
    };

    PlaybackSetting(Type type, float value)
    { this->type = type; this->value = value; }

    Type getType() const { return type; }

    /** Returns the stored value for this setting. Values are always stored as floats and it is 
    understood that in cases, where the corresponding parameter in the sfz spec is defined to be an
    integer, we just represent it as that integer with all zeros after the decimal dot and the 
    caller is supposed to convert to int by writing e.g.:

      int intValue = (int)setting.getValue();

    For settings whose value is represented by text that indciates a particular choice, the caller
    has to look up the int in an enumeration corresponding to the type of the parameter...tbc... */
    float getValue() const { return value; }
    // todo: (decide and) document, how choice parameters like filter-type are represented. In sfz,
    // they are just represented as text. Maybe for each such choice parameter, we need another 
    // enum to represent its allowed values.

    /** Some settings need to specify an index in addtion to the value. An example is a setting 
    involving midi control changes. In the sfz file, they are written as e.g. loccN where the N is
    replaced by the actual controller number, like locc74=20 to indicate, that the sample should
    only play, if the last received value for CC#74 was >= 20. If indexing not applicable to a 
    particular setting/opcode, this will return -1. */
    int getIndex() const { return index; }

  private:

    Type  type  = Type::Unknown;
    float value = 0.f;
    int   index = -1;  //< Used e.g. for conrol-change settings. Is -1, if not applicable.
  };
  // Maybe rename to Opcode - but no: "opcodes" are the strings that appear in the sfz file, such
  // "lokey". They map to the Type of the playback setting. Maybe this class should provide the
  // mapping (maybe std::map or some selfmade class for a 2-way associative array)


  //-----------------------------------------------------------------------------------------------
  /** Baseclass for the 3 organizational levels of the sfz specification, factoring out their 
  commonalities. Subclasses are Region, Group, Instrument. */
  class OrganizationLevel
  {

  public:

    /** Returns a const reference to our playback settings. */
    const std::vector<PlaybackSetting>& getSettings() const { return settings; }

    /** Returns a pointer to the parent level which encloses this level. In a Region, this would 
    point to its enclosing Group, in a Group to its enclosing Instrument and in an Instrument, it
    would remain nullptr (unless we introduce an even higher level such as an "Ensemble"). */
    const OrganizationLevel* getParent() const { return parent; }

    // todo: float getSetting(PlaybackSetting::Type, int index) this should loop through the 
    // settings to see, if it finds it and if not call the same method on the parent or return the
    // default value, if the parent is nullptr.

  protected:

    OrganizationLevel* parent = nullptr;  

    void clearSettings() { settings.clear(); }

    std::vector<PlaybackSetting> settings;

  };


  //-----------------------------------------------------------------------------------------------
  /** A region is the lowest organizational level ins sfz. It contains a sample along with 
  performance settings including information for which keys and velocities the sample should be 
  played and optionally other constraints for when the the sample should be played and also 
  settings for pitch, volume, pan, filter, envelopes, etc. */
  class Region : public OrganizationLevel
  {

  public:

    // replaces getSampleStream

    /** Returns the generic pointer for custom satellite data or objects that are associated with
    this region. This pointer is intended to be used for some sort of audio stream object that is 
    used for accessing the sample data. It has been made a generic void pointer to decouple 
    rsDataSFZ from the AudioFileStream class that is used in rsSamplerEngine. The sampler-engine 
    assigns this pointer with appropriate stream object and when retriveing them, does an 
    appropriate type cast. ToDo: try to find a better design, maybe move up into baseclass */
    const void* getCustomPointer() const { return custom; }

    /** Return a pointer to the group to which this region belongs. */
    const Group* getGroup() const { return (const Group*) getParent(); }

    //const Group* getGroup() const { return group; }
    // todo: return (const Group*) getParentLevel();


    /** Returns the lowest key at which this region will be played. */
    uchar getLoKey() const { return loKey; }

    /** Returns the higest key at which this region will be played. */
    uchar getHiKey() const { return hiKey; }

    /** Returns the lowest velocity at which this region will be played. */
    uchar getLoVel() const { return loVel; }

    /** Returns the highest velocity at which this region will be played. */
    uchar getHiVel() const { return hiVel; }

  private:

    //Group* group = nullptr;  //< Pointer to the group to which this region belongs
    // mayb get rid by having a general parentLevel pointer defined in the baseclass


    // todo: setters for loKey,...

    /** Sets the audio stream object that should be used for this region. */
    //void setSampleStream(const AudioFileStream* newStream) { sampleStream = newStream; }

    void setCustomPointer(const void* newPointer) { custom = newPointer; }


    //const AudioFileStream* sampleStream = nullptr;
    // try to get rid - that member should be added by rsSamplerEngine::Region which should be
    // a subclass of rsInstrumentDataSFZ::Region, and/or move up into baseclass. maybe to decouple
    // rsDataSFZ from AudioFileStream, just keep it as pointer-to-void which client code may 
    // typecast to any sort of stream...or maybe that coupling makes sense?..hmm - not really.
    // maybe a pointer-to-void named customData should be stored in OrganizationLevel

    const void* custom = nullptr;





    uchar loKey = 0, hiKey = 127;
    uchar loVel = 0, hiVel = 127;
    // todo: maybe package loKey/hiKey, loVel/hiVel into a single uchar to save memory.
    // To prepare for this, provide get/setLoKey() etc. accessors and use them consistently in
    // rsSamplerEngine. Should these be moved into the baseclass, meaning that groups and 
    // instruments can also restrict the keyrange additionally? Or is this opcode really 
    // specifically applicable to regions only? Test with SFZPlayer and replicate its behavior.
    // I think, it could be useful to restrict keyranges of groups and even instruments, when
    // they are part of an enseble - for example, for keyboard splits.


    //std::string name; sample

    friend class Group;  // do we need this? if not, get rid.
    friend class rsSamplerEngine;  // try to get rid
    // The Region class shall not provide any public functions that can modify the region because
    // those could be used by client code to modify the region behind the back of the 
    // rsSamplerEngine which could mess things up. Client code can modify regions only through the
    // appropriate functions of rsSamplerEngine. It acts as man-in-the-middle and can the call the
    // private setters of the Region (by virtue of being a friend class) and it may also trigger 
    // additional actions, if necessary. The same should probably apply to the Group class as well.
    // Is this a known pattern? -> figure out
  };

  //-----------------------------------------------------------------------------------------------
  /** A group organizes a bunch of regions into a single entity for which performance settings can 
  be set up which will be applicable in cases where the region does not itself define these 
  settings, so they act as fallback values. It's the mid-level of organization in sfz. */
  class Group : public OrganizationLevel
  {

  public:

    ~Group() { clearRegions(); }

    /** Returns the index of the given region within this group, if present or -1 if the region is
    not present in this group. */
    int getRegionIndex(const Region* region) const;

    /** Returns true, if the given index i refers toa valid region within this group. */
    bool isRegionIndexValid(int i) const { return i >= 0 && i < (int)regions.size(); }


    /** Return a pointer to the instrument to which this group belongs. */
    const Instrument* getInstrument() const { return (const Instrument*) getParent(); }


    /** Returns a pointer to the region with the given index within the group. */
    Region* getRegion(int i) const;


  private:

    //Group* Instrument = nullptr;  //< Pointer to the instrument to which this group belongs

    int addRegion();      // todo: removeRegion, etc.
    void clearRegions();

    std::vector<Region*> regions;
    /**< Pointers to the regions belonging to this group. */

    // may be add these later:
    //std::string name;  

    friend class rsSamplerEngine;  // try to get rid
  };

  // todo: class Instrument - group should have a pointer to its enclosing instrument

  //-----------------------------------------------------------------------------------------------
  /** The instrument is the highest organizational level in sfz. There is actually no section 
  header in the .sfz file format for the whole instrument that corresponds to this class. The whole
  content of the sfz file *is* the instrument. But for consistency, we represent it by a class as 
  well. Maybe later an additional "Ensemble" or "Orchestra" level can be added on top. */
  class Instrument : public OrganizationLevel
  {

  public:

  private:

    std::vector<Group> groups;
    // Should that be an array of pointers, too? Like the regions array in Group? That would make
    // the implementations of Group and Instrument more consistent but is actually technically not 
    // necessary. So, for the time being, let's keep it an array of direct value objects.

  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  // todo: factor out some code from rsSamplerEngine - the names of the correponding setters here 
  // and there should match and rsSamplerEngine should call functions from here and perhaps do
  // additional stuff, if necessary


protected:

  Instrument instrument; // Maybe we could maintain an array of such isntruments

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

  rsSamplerEngine(int maxNumLayers = 16);
  // todo: use some higher default value - what is reasonable here needs some testing in realistic
  // scenarios

  virtual ~rsSamplerEngine();


  // for convenience:
  using uchar = unsigned char;
  using Region = rsDataSFZ::Region; // todo: make a subclass here that adds the stream field
  using Group  = rsDataSFZ::Group;
  using PlaybackSetting = rsDataSFZ::PlaybackSetting;



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
  // rename voiceOverload to layerOverload. in sfz, a voice refers to asingle key which can contain
  // multiple layers/regions
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

  /** Sets a value for a given type of playback setting for the given region. Returns either
  ReturnCode::success or ReturnCode::notFound, if the region was not found in this instrument. If 
  this happens, it indicates a bug on the call site. */
  int setRegionSetting(Region* region, PlaybackSetting::Type type, float value);

  // todo: setGroupSetting, setInstrumentSetting, removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

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


  int getMaxNumLayers() const { return (int) playerPool.size(); }

  int getNumActiveLayers() const { return (int) activePlayers.size(); }

  int getNumIdleLayers() const { return (int) idlePlayers.size(); }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(float* left, float* right);

  void processBlock(float** block, int numFrames);

  void handleMusicalEvent(const rsMusicalEvent<float>& ev);

  // void processFrameVoice, processBlockVoice


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internal Helper Classes

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

    /** A basic sanity check for the given region. Mostly for catching bugs. */
    virtual bool isPlayable(const Region* region);

    /** Sets up the internal values for the playback settings (including DSP objects) according
    to the assigned region and resets all DSP objects. */
    virtual void prepareToPlay();

    virtual void resetDspState();
    virtual void resetDspSettings();
    virtual void setupDspSettings(const std::vector<PlaybackSetting>& settings);

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

  /** Returns the AudioFileStream object that is used to stream the actual sample data for the
  given region. A pointer to this object is supposed to be stored within the region object
  ...tbc... */
  static const AudioFileStream* getSampleStreamFor(const Region* r);

  /** Handles a noteOn event with given key and velocity and returns either ReturnCode::success, if
  we had enough voices available to serve the request or ReturnCode::voiceOverload, in case the 
  noteOn could not be handled due to inavailability of a sufficient number of idle voices. If no
  sufficient number of idle voices was available and the noteOn should actually have triggered 
  playback of multiple samples, none of them will be triggered. It's an all-or-nothing thing: we 
  don't ever trigger playback for only a subset of samples for a given noteOn. */
  int handleNoteOn(uchar key, uchar vel);

  /** Analogous to handleNoteOn. It may also return ReturnCode::voiceOverload in cases where the 
  noteOff is supposed to trigger relase-samples. In such a case, none of the release-samples will 
  be triggered. */
  int handleNoteOff(uchar key, uchar vel);

  //-----------------------------------------------------------------------------------------------
  // \name Data

  rsDataSFZ sfz;
  /**< The data structure that defines the sfz instrument. */
 
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
  // maybe use a std::vector, maybe we need a similar array for note-off samples

  SamplePool samplePool;
  /**< The pool of samples that are in use for the currently loaded instrument. The samples are 
  pooled to avoid redundant storage in memory when multiple regions use the same sample. */

  std::vector<Group> groups;
  /**< The groups contained in this instrument. Each group may contain set of regions. */
  // get rid - should go into rsInstrumentDataSFZ::Instrument

  std::vector<PlaybackSetting> settings;
  /**< Playback settings that apply to all groups within this instrument, unless a group (or 
  region) overrides a setting with its own value. **/
  // get rid - should go into sfz.instrument.settings

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

  using rsSamplerEngine::rsSamplerEngine;  // inherit constructors

  static int getRegionPlayerSize() { return sizeof(rsSamplerEngine::RegionPlayer); }

};


#endif