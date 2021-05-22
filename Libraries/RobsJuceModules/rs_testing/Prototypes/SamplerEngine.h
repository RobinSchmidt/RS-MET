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

template<class T>
class AudioStream
{

public:

  virtual ~AudioStream() {}

  /** For random access. Writes the sample frame with given index into the given destination. */
  virtual void getFrame(int sampleIndex, T* destination) const = 0;

  /** Function for speficically handling stereo signals to allow handling that importnat, common 
  special case more efficiently than with the more general implementation. */
  virtual void getFrameStereo(int sampleIndex, T* left, T* right) const = 0;

  /** Sets the number of output channels for this object. By default, this number will be equal to
  the number of channels in the data, as set by the setData call. However, the number of desired 
  output channels to be produced may actually be different from that. For example, we may want to
  produce stereo output from mono data by just copying the single output into both channels. Let
  M be the number of channels in the data and N be the number of output channels of the stream, 
  then we will produce: out[n] = data[m % N] where n is the output channel index and m % N the data
  channel index. */
  //virtual void setNumOutputChannels(int newNumChannels) = 0;

  /** Returns the number of sample frames in this stream. */
  int getNumFrames() const { return numFrames; }


  T getSampleRate() const { return sampleRate; }

  //virtual bool hasFinished(int sampleTime)

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

  T   sampleRate  = T(44100);
  int numChannels = 0;
  int numFrames   = 0;         // maybe use -1 to encode "unknown"? would that be useful?

};
// maybe move out of this class - this may be useful in other contexts, too - maybe templatize


/*
class rsAudioFileInfo  // maybe rename to rsFileInfo for use in more general context
{

public:

  bool operator==(const rsAudioFileInfo& rhs) const 
  { 
    bool same = true;
    same &= fileName     == rhs.fileName;
    same &= extension    == rhs.extension;
    same &= path         == rhs.path;
    same &= rootDirIndex == rhs.rootDirIndex;
    return same;
  }

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
*/


template<class T>
class AudioFileStream : public AudioStream<T>
{

public:

  // ToDo: add == operator based on fileName, extension, path (and maybe (meta)data such as 
  // sampleRate, numChannels, numFrames, ...)

  // inquiry: existsOnDisk (idea: we should allow to work with samples that are not stored on 
  // disk but rather pre-rendered programmatically into RAM at runtime)

  //virtual void getFrame(int sampleIndex, TSmp* destination) {}

  // todo: getBlock

  /** Returns true, iff the given other stream object uses the same audio file as this. */
  bool usesSameFile(const AudioFileStream<T>* otherStream) const
  {
    //return fileInfo == otherStream->fileInfo;

    bool same = true;
    //same &= fileName     == otherStream->fileName;
    //same &= extension    == otherStream->extension;
    same &= path         == otherStream->path;
    same &= rootDirIndex == otherStream->rootDirIndex;
    return same;
  }

  const std::string& getPath() const { return path; }
  //std::string getPath() const { return path + fileName + extension; }
  // todo: 
  // -don't split the path into 3 fields -> have only a path field which is the full path
  // -return a const reference to that here



protected:

  //rsAudioFileInfo fileInfo;


   // factor out into class rsFileInfo
  //std::string fileName;  // without filename extension (e.g. Piano_A4)
  //std::string extension; // filename extension (e.g. wav, flac)
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

  // todo: maybe keep just the path which should represents the full relative path. It doesn't seem
  // to be a good idea to split it into 3 parts
};


// maybe rename to AudioFileStreamRAM, another subclass can be named AudioFileStreamDFD,
// or ...StreamFromMemory/FromDisk

template<class T>
class AudioFileStreamPreloaded : public AudioFileStream<T>
{

public:


  virtual ~AudioFileStreamPreloaded() { clear(); }


  /** Returns true, iff everything went alright and false if it failed to allocate the required
  memory. */
  bool setData(T** newData, int numFrames, int numDataChannels, T sampleRate, 
    int numStreamChannels, const std::string& path);
  // todo: 
  // -we need to distiguish between the number of channels in the data and the desired number of 
  //  output channels
  // -include fileName etc. later, too

  //void setNumOutputChannels(int newNumChannels) override;


  void clear();


  void getFrame(int sampleIndex, T* destination) const override
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


  void getFrameStereo(int sampleIndex, T* left, T* right) const override
  {
    rsAssert(numChannels == 2); // Can be used only for stereo signals
    int n = sampleIndex;
    *left  = channelPointers[0][n];
    *right = channelPointers[1][n];
  }

  //void getFrameStereo(int sampleIndex, float* destination) const

protected:

  T*  flatData = nullptr;         // pointer to the sample data
  T** channelPointers = nullptr;  // pointers to the channels
  // If we store the data in interleaved format, the channelPointers will be not needed and 
  // getFrame must be implemented differently. Maybe that's better (more efficient)
};

//=================================================================================================

/** A class to represent a pool of audio samples...tbc...
ToDo: in a more general context, we would probably need a more complex referencing system for the 
AudioStream objects - regions would need to be something like AudioStreamClients, that 
register/deregister themselves, etc. Maybe AudioStream should be a subclass of some DataStream 
baseclass. We'll see... */

template<class T>
class SamplePool
{

public:


  ~SamplePool() { clear();  }


  int addSample(const AudioFileStream<T>* newSample)
  {
    // rsAssert(!contains(newSample))
    samples.push_back(newSample);
    return ((int) samples.size()) - 1;
  }
  // should the pool take ownership? ...i think so

  void removeSample(int i)
  {
    rsAssert(i >= 0 && i < (int)samples.size());
    delete samples[i];
    rsRemove(samples, (size_t)i);
  }


  /** Returns the number of samples in this pool. */
  int getNumSamples() const { return (int)samples.size(); }


  /** Returns true, if the given index i refers to a valid sample index. */
  bool isSampleIndexValid(int i) const { return i >= 0 && i < (int)samples.size(); }


  const AudioFileStream<T>* getSampleStream(int i) const
  {
    if(!isSampleIndexValid(i)) {
      rsError("Invalid sample index");
      return nullptr; }
    return samples[i];
  }


  const std::string& getSamplePath(int i) const { return getSampleStream(i)->getPath();  }

  void clear();

  // todo:
  // setup: removeSample...but if regions refer to it, we need to update them, too by 
  // invalidating their pointers. We either need an observer mechanism (complex and general) or 
  // we allow removal of samples only via a member function of the outlying rsSamplerEngine 
  // class, which also takes care of resetting the sample-streams in all regions that use it
  // (simpler but less general). Maybe to make it safer, we could also introduce a reference
  // counter and check, if it is zero, before a stream objects gets removed
  // inquiry: hasSample

protected:

  std::vector<const AudioFileStream<T>*> samples;

};
// maybe templatize, rename to AudioFileStreamPool


//=================================================================================================

/** Data structure to define sample based instruments conforming to the sfz specification. 

ToDo: 
-use size_t or int consistently for indexing groups and regions
-use pointers or references consistently for returning sub-levels
*/

class rsSamplerData // todo: move into its own pair of .h/.cpp files, rename to rsSamplerData
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Default constructor. */
  rsSamplerData() { }

  ~rsSamplerData() { clearInstrument(); }

  /** Copy constructor.   */
  rsSamplerData(const rsSamplerData& d) { copy(d, *this); }

  /** Copy assignment operator. */
  rsSamplerData& operator=(const rsSamplerData& d) { if(this != &d) copy(d, *this); return *this; }






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

      Volume, Pan, PanRule,
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
      off, lp_6, lp_12, hp_6, hp_12, bp_6_6, br_6_6,
      
      numFilterTypes
    };


    enum PanRule  
    {
      linear, sinCos,

      numPanRules
    };
    // Or maybe it should be called PanLaw? Maybe have different variations with respect to total
    // gain - for linear: either factor 2 for hard left/right setting or a factor of 0.5 for a 
    // center setting. The former would imply that with neutral default settings, stereo samples 
    // are played as is. The latter would imply that hard-panned mono samples would be played as is
    // on their respective channel. Both behaviors may be useful, although, it would be a bit 
    // redundant because we also have an overall gain parameter as well which can always be used to 
    // compensate...although a factor of exactly 2 or 0.5 may be hard to achieve because gain is 
    // given in dB, so the sfz file would have to specify +-6.0205999132796239....., which is 
    // inconvenient


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


    bool operator==(const PlaybackSetting& rhs) const
    { return type == rhs.type && value == rhs.value && index == rhs.index; }


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


    /** Sets the sample to be used for this region. This should be a string that represents the 
    path of the sample relative to some fixed root directory (typically the directory of the sfz 
    file). */
    void setSamplePath(const std::string& newPath) { samplePath = newPath; }


    void copyDataFrom(const OrganizationLevel* lvl);



    /** @see setSample */
    const std::string& getSamplePath() const { return samplePath; }

    /** Returns a const reference to our playback settings. */
    const std::vector<PlaybackSetting>& getSettings() const { return settings; }

    /** Returns a pointer to the parent level which encloses this level. In a Region, this would 
    point to its enclosing Group, in a Group to its enclosing Instrument and in an Instrument, it
    would remain nullptr (unless we introduce an even higher level such as an "Ensemble"). */
    const OrganizationLevel* getParent() const { return parent; }



    void addSetting(const PlaybackSetting& s) { settings.push_back(s); }
    // Maybe we should have a function setSetting that either adds a new setting or overwrites
    // an existing one


    /** Returns the generic pointer for custom satellite data or objects that are associated with
    this region. This pointer is intended to be used for some sort of audio stream object that is 
    used for accessing the sample data. It has been made a generic void pointer to decouple 
    rsSamplerData from the AudioFileStream class that is used in rsSamplerEngine. The sampler-engine 
    assigns this pointer with appropriate stream object and when retriveing them, does an 
    appropriate type cast. ToDo: try to find a better design, maybe move up into baseclass */
    const void* getCustomPointer() const { return custom; }
    // replaces getSampleStream


    // todo: float getSetting(PlaybackSetting::Type, int index) this should loop through the 
    // settings to see, if it finds it and if not call the same method on the parent or return the
    // default value, if the parent is nullptr.

  protected:


    /** Sets the audio stream object that should be used for this region. */
    //void setSampleStream(const AudioFileStream* newStream) { sampleStream = newStream; }

    void setCustomPointer(const void* newPointer) { custom = newPointer; }

    //const AudioFileStream* sampleStream = nullptr;
    // try to get rid - that member should be added by rsSamplerEngine::Region which should be
    // a subclass of rsInstrumentDataSFZ::Region, and/or move up into baseclass. maybe to decouple
    // rsSamplerData from AudioFileStream, just keep it as pointer-to-void which client code may 
    // typecast to any sort of stream...or maybe that coupling makes sense?..hmm - not really.
    // maybe a pointer-to-void named customData should be stored in OrganizationLevel

    std::string samplePath; 
    // This is the full (relative) path. ToDo: maybe this should be moved into the baseclass. We'll
    // need to figure out, if the sfz player allows samples to be defined also for groups. The same
    // goes for the loKey,etc. stuff as well. Actually, the string is redundant here when the 
    // "custom" pointer actually points to an AudioFileStream, because that stream object also 
    // stores the path. Maybe revert the custom void pointer to a pointer-to-AudioFileStream again. 
    // Yes, this will introduce coupling but it gets rid of the redundant storage. Maybe trying to 
    // decouple it amounts to the "speculative generality" antipattern here - we'll see...
    // https://refactoring.guru/smells/speculative-generality

    const void* custom = nullptr;

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



    void setLoKey(uchar newKey) { loKey = newKey; }

    void setHiKey(uchar newKey) { hiKey = newKey; }


    void copyDataFrom(const Region* src);


    /** Return a pointer to the group to which this region belongs. */
    const Group* getGroup() const { return (const Group*) getParent(); }
    // rename to getParentGroup or getEnclosingGroup

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

    bool shouldPlayForKey(uchar key) const { return key >= loKey && key <= hiKey; }
    // todo: later maybe also take loKey/hiKey of group (and instrument) into account, requires to
    // move loKey/hiKey members into baseclass

    bool operator==(const Region& rhs) const;


  private:

    //Group* group = nullptr;  //< Pointer to the group to which this region belongs
    // mayb get rid by having a general parentLevel pointer defined in the baseclass

    // todo: setters for loKey,...

    uchar loKey = 0, hiKey = 127;
    uchar loVel = 0, hiVel = 127;
    // todo: maybe package loKey/hiKey, loVel/hiVel into a single uchar to save memory.
    // To prepare for this, provide get/setLoKey() etc. accessors and use them consistently in
    // rsSamplerEngine. Should these be moved into the baseclass, meaning that groups and 
    // instruments can also restrict the keyrange additionally? Or is this opcode really 
    // specifically applicable to regions only? Test with SFZPlayer and replicate its behavior.
    // I think, it could be useful to restrict keyranges of groups and even instruments, when
    // they are part of an enseble - for example, for keyboard splits.
    // -maybe these should go into the baseclass as well
    // -maybe the getters should use min/max with the stored settings and those of the parent such
    //  that regions can only further restrict the the range...or maybe they should override the
    //  group setting -> check, how sfzPlayer behaves

    friend class Group;  // do we need this? if not, get rid.
    friend class rsSamplerData;
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


    int addRegion(uchar loKey = 0, uchar hiKey = 127);      // todo: removeRegion, etc.

    int addRegion(Region* newRegion); 

    void copyDataFrom(const Group* scrGroup);

    void clearRegions();



    /** Returns the index of the given region within this group, if present or -1 if the region is
    not present in this group. */
    int getRegionIndex(const Region* region) const;

    /** Returns true, if the given index i refers toa valid region within this group. */
    bool isRegionIndexValid(int i) const { return i >= 0 && i < (int)regions.size(); }

    size_t getNumRegions() const { return regions.size(); }


    /** Return a pointer to the instrument to which this group belongs. */
    const Instrument* getInstrument() const { return (const Instrument*) getParent(); }


    /** Returns a pointer to the region with the given index within the group. */
    Region* getRegion(int i) const;
    // rename to getRegionPointer or Ptr


    bool operator==(const Group& rhs) const;
    // the comparison is quite strict in the sense that the settings must occur in the same order



  private:




    //std::vector<Region> regions;
    std::vector<Region*> regions;
    /**< Pointers to the regions belonging to this group. */

    friend class rsSamplerData;
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

    size_t getNumGroups() const { return groups.size(); }


    /** Returns a pointer to the region with the given index within the group. */
    Group* getGroup(int i) { return groups[i]; }


    const std::vector<PlaybackSetting>& getGroupSettings(size_t groupIndex) const
    { return groups[groupIndex]->getSettings(); }


    bool operator==(const Instrument& rhs) const;
    //{ return settings == rhs.settings && groups == rhs.groups; }


  //private:  // make protected later

    // maybe move to public
    int addGroup();      // todo: removeGroup, etc.

    int addGroup(Group* g);


    void clearGroups();


    std::vector<Group*> groups;
    // Should that be an array of pointers, too? Like the regions array in Group? That would make
    // the implementations of Group and Instrument more consistent but is actually technically not 
    // necessary. So, for the time being, let's keep it an array of direct value objects.

    friend class rsSamplerData;
  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  // todo: factor out some code from rsSamplerEngine - the names of the correponding setters here 
  // and there should match and rsSamplerEngine should call functions from here and perhaps do
  // additional stuff, if necessary

  int addGroup() { return instrument.addGroup(); }

  int addGroup(Group* newGroup) { return instrument.addGroup(newGroup); }


  int addRegion(int gi, uchar loKey = 0, uchar hiKey = 127);

  void setRegionCustomPointer(int gi, int ri, void* ptr)
  { instrument.groups[gi]->regions[ri]->setCustomPointer(ptr); }

  void setRegionSample(int gi, int ri, const std::string& samplePath)
  { instrument.groups[gi]->regions[ri]->setSamplePath(samplePath); }

  void setRegionSetting(int gi, int ri, PlaybackSetting::Type type, float value)
  {
    instrument.groups[gi]->regions[ri]->settings.push_back(PlaybackSetting(type, value));
    // Preliminary. We need to figure out, if that setting already exists and if so, just change 
    // its value instead of pushing another value for the same parameter
  }



  /** Clears the whole instrument definition. */
  void clearInstrument() 
  { instrument.clearGroups(); }
  //{ instrument.groups.clear(); }
  // todo: may wrap into instrument.clear() - don't access the groups array directly


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  size_t getNumGroups() const { return instrument.getNumGroups(); }


  //const Group& getGroupRef(size_t i) const { return instrument.groups[i]; }

  const Group& getGroupRef(size_t i) const { return *instrument.groups[i]; }
  // replace by getGroupPtr ..maybe renma to just getGroup


  const Group* getGroupPtr(size_t i) const { return instrument.groups[i]; }

  const Region* getRegionPtr(size_t gi, size_t ri) const 
  { return instrument.groups[gi]->regions[ri]; }

  /** Returns a const reference to the playback settings if the i-th group. */
  //const std::vector<PlaybackSetting>& getGroupSettings(size_t i) const 
  //{ return instrument.getGroupSettings(i); }


  bool operator==(const rsSamplerData& rhs) const { return instrument == rhs.instrument; }

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Produces the string that represents the settings in an sfz-file compliant format, i.e. a 
  string that can be written into an .sfz file. */
  std::string getAsSFZ() const;


  /** Sets up this data object according to the given string which is supposed to represent the 
  contents of an .sfz file. */
  void setFromSFZ(const std::string& sfzFileContents);
  // todo: return a return-code, including unknownOpcode, invalidValue, invalidIndex, ...


  bool saveToSFZ(const char* path) const;
  // todo: return a return-code, including fileWriteError

  bool loadFromSFZ(const char* path);
  // todo: return a return-code, including sfzFileNotFound, sampleFileNotFound



//protected:  // preliminarily commented - make protected again later

  Instrument instrument; 
  // Maybe we could maintain an array of such instruments that define an ensmeble

protected:

  static void writeSettingToString(const PlaybackSetting& setting, std::string& str);

  static PlaybackSetting getSettingFromString(
    const std::string& opcode, const std::string& value);

  static void copy(const rsSamplerData& src, rsSamplerData& dst);

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
  using Region = rsSamplerData::Region; // todo: make a subclass here that adds the stream field
  using Group  = rsSamplerData::Group;
  using PlaybackSetting = rsSamplerData::PlaybackSetting;



  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Return codes for the setup functions. We use encodings as negative integers so we can use 
  them also for functions which use positive integers as valid return values. */
  enum ReturnCode
  {
    success        = -1,  //< Operation completed successfully. 
    nothingToDo    = -2,  //< There was nothing to actually do. State was already as desired.
    memAllocFail   = -3,  //< Memory allocation failure.
    invalidIndex   = -4,  //< An invalid index was passed.
    voiceOverload  = -5,  //< Not enough free voices available (in e.g. new noteOn).
    notFound       = -6,  //< A region, group or whatever was not found.
    fileLoadError  = -7,  //< A file could not be loaded (reasons: not found or failed alloc).

    notImplemented = -8   //< Feature not yet implemented (relevant during development).
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
    const std::string& path);
  // Maybe rename to addSample, it should return the index of the sample in the sample-pool
  // maybe make a struct SampleMetaData containing: numFrames, numChannels, sampleRate, rootKey
  // todo: take reference to a metaData object

  /** Loads a sample represented by the given path into our samplePool. The path is supposed to be 
  relative to some fixed root directory which is typcally the directory in which sfz file resides,
  but later this may be switched to some user and/or factory content directory, too (this requires
  to introduce a new opcode to sfz...maybe root_dir or sample_directory or sample_folder or 
  something which should be defined once for the whole instrument). The returns the following
  return codes: 
    success:       sample was succesfully loaded into the pool
    nothingToDo:   sample was already in the pool
    fileLoadError: sample could not be loaded (maybe the path was wrong?)
    memAllocFail:  we could not allocate enough memory to load the sample
  ToDo: verify return codes in unit test  */
  int loadSampleToPool(const std::string& path);



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

  /** Sets a value for a given type of playback setting a region. Returns either
  ReturnCode::success or ReturnCode::invalidIndex, if groupIndex and/or regionIndex was invalid. If 
  this happens, it indicates a bug on the call site. */
  int setRegionSetting(int groupIndex, int regionIdex, PlaybackSetting::Type type, float value);


  // todo: setGroupSetting, setInstrumentSetting, removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,

  /** Sets up the engine from the given sfz data object and returns ReturnCode::success, if all
  is well or...  */
  int setupFromSFZ(const rsSamplerData& sfz);


  /** Writes the current instrument definition into an sfz file with given path. */
  bool saveToSFZ(const char* path) const { return sfz.saveToSFZ(path); }
  // -document, whether path is absolute or relative and if the latter, what is the root
  // -return a return-code instead of bool
  // -maybe move elsewhere

  /** Loads the instrument definition given by an sfz file with the given path. Returns 
  ReturnCode::success if all wen well or ReturnCode::fileLoadError if loading of the sfz or any of 
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
    return gi >= 0 && gi < (int)sfz.instrument.groups.size() 
      && sfz.instrument.groups[gi]->isRegionIndexValid(ri);
  }

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


  int getMaxNumLayers() const { return (int) playerPool.size(); }

  int getNumActiveLayers() const { return (int) activePlayers.size(); }

  int getNumIdleLayers() const { return (int) idlePlayers.size(); }

  /** Returns a const pointer to the rsSamplerData object that represents the current instrument
  settings. */
  const rsSamplerData& getInstrumentData() const { return sfz; }
    



  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(double* left, double* right);

  void processFrame(float* left, float* right);

  void processBlock(float** block, int numFrames);

  void handleMusicalEvent(const rsMusicalEvent<float>& ev);

  // void processFrameVoice, processBlockVoice

  /** Stops the playback of all currently active RegionPlayers immediately. This is a rather hard 
  reset which may be appropriate to call when a midi reset message is received or before loading a
  new patch. It returns the number of players that were affected, i.e. the number of players that 
  were in active state before the call. */
  int stopAllPlayers();

  /** Calls stopAllPlayers. Function is for consistency with the rest of the library. */
  void reset() { stopAllPlayers(); }




protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internal Helper Classes

  class SignalProcessor
  {
  public:
    virtual void processFrame(rsFloat64x2& inOut) = 0;
    virtual void resetState() = 0;
    virtual void resetSettings() = 0;
    // todo: processBlock
  };

  class Modulator
  {
  public:
    virtual double getSample() = 0;
    virtual void resetState() = 0;
    virtual void resetSettings() = 0;
    // todo: processBlock
  };

  class ModulationConnection
  {

  private:
    std::function<void(double)> targetSetter; 
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
    virtual void setRegionToPlay(const Region* regionToPlay, double outputSampleRate);

    /** Sets the midi note number for which this player was started. This needs to be set up when 
    receiving a noteOn. This information is used later when receiving a noteOff to identify which 
    players need to stop. */
    void setKey(uchar newKey) { key = newKey; }

    /** Generates one stereo sample frame at a time. */
    virtual rsFloat64x2 getFrame();

    /** Writes a block of given length into the outBuffer. */
    virtual void processBlock(rsFloat64x2* outBuffer, int length);

    /** Returns true, iff this player has finished its playback job, for example by having reached
    the end of the sample stream and/or amplitude envelope. */
    virtual bool hasFinished(); // should be const?

    /** Retrieves the information about the midi note for which this player was started. Used to 
    identify players that need to stop, when a noteOff is received. @see setKey */
    uchar getKey() const { return key; }


  protected:

    /** A basic sanity check for the given region. Mostly for catching bugs. */
    virtual bool isPlayable(const Region* region);

    /** Sets up the internal values for the playback settings (including DSP objects) according
    to the assigned region and resets all DSP objects. */
    virtual void prepareToPlay(double sampleRate);

    virtual bool buildProcessingChain();
    virtual void resetDspState();
    virtual void resetDspSettings();
    virtual void setupDspSettings(const std::vector<PlaybackSetting>& settings, double sampleRate);

    const Region* region;                 //< The Region object that this object should play
    const AudioFileStream<float>* stream; //< Stream object to get the data from
    rsFloat64x2 amp = 1.0;                //< Amplitude (for both channels)
    //int sampleTime = 0;            //< Elapsed time in samples, negative values used for delay
    double sampleTime = 0.0;       //< Time index in the sample. Negative values used for delay.
    double increment  = 1.0;       //< Increment of sampleTime per sample
    uchar key = 0;                 //< Midi note number used for starting this player

    std::vector<Modulator*> modulators;
    std::vector<SignalProcessor*> dspChain;
    std::vector<ModulationConnection*> modMatrix;  // not a literal matrix but conceptually

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
  "idlePlayers". The return value is either ReturnCode::success or ReturnCode::invalidIndex, if
  the activeIndex was not a valid index into our activePlayers array. */
  int deactivateRegionPlayer(size_t activeIndex);

  /** Returns the AudioFileStream object that is used to stream the actual sample data for the
  given region. A pointer to this object is supposed to be stored within the region object
  ...tbc... */
  static const AudioFileStream<float>* getSampleStreamFor(const Region* r);

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


  /** Removes those samples from our sample pool that are not used in the given sfz instrument 
  specification. Returns the number of samples that were removed. */
  int removeSamplesNotUsedIn(const rsSamplerData& sfz);
  // maybe rename to removeUnusedSamples. But that name is more ambiguous: it could be interpreted
  // as "unused in the current sfz member", so maybe don't

  /** Adds all samples to our sample pool that are used in the given sfz instrument definition, if
  they are not already there. Returns the number of samples that were added or 
  ReturnCode::fileLoadError if any of the files failed to load. */
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

/** A class for setting up an object of class rsSamplerEngine according to a string representing
an .sfz file. */

class rsSamplerEngineLoaderSFZ  // maybe rename to Serializer
{

public:

  /** Sets up the given rsSamplerEngine object according to the given string which is supposed to
  represent the contents of an .sfz file. */
  static void setFromSFZ(rsSamplerEngine* engine, const std::string& sfzFileContents);
  // todo: return a return-code, including unknownOpcode, invalidValue, invalidIndex, ...

  /** Given an rsSamplerData object, this function produces the string that represents the settings in
  an sfz-file compliant format, i.e. a string that can be written into an .sfz file. */
  static std::string getAsSFZ(const rsSamplerData& sfz); 

};

// hmmm...maybe this additional class will not be needed after all - provide de/serialize in
// rsSamplerData and maybe also in rsSamplerEngine, where serialize just forwards to sfz.serialize and
// deserialize may have to take additional actions
// Maybe this should not be a separate class and rsSamplerEngine should just have a pair of 
// functions getAsSFZ/setFromSFZ. We'll see, how complex the code gets. If it's not too complex,
// integrate it into rsSamplerEngine. Or maybe the code should go into rsSamplerData...at least, the
// code related to generating and parsing sfz strings, maybe not the code related to actually
// setting up the engine object such that the rsSamplerData remains independent from rsSamplerEngine.




//=================================================================================================

/** Subclass that contains some extra functions that facilitate testing which should not go into 
the production code. */
class rsSamplerEngineTest : public rsSamplerEngine
{

public:

  using rsSamplerEngine::rsSamplerEngine;  // inherit constructors

  static int getRegionPlayerSize() { return sizeof(rsSamplerEngine::RegionPlayer); }

};


#endif