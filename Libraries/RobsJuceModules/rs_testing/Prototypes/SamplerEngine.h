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
"instrument" itself. Just like groups provide fallback settings for region, the whole instrument 
can provide fallback settings for all the groups it contains. If some performance parameter isn't 
defined anywhere (neither in the instrument, group or region), a neutrally behaving default value 
will be used.  */

class rsSamplerEngine
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Helper Classes

  using uchar = unsigned char;
  class Region;
  class AudioFileStream;

  /*
  class SampleMetaData
  {

  public:

    
  private:

    std::string uniqueName;  // Unique within the context of the instrument definition.
    std::string fileName;    // Without extension for the format.
    std::string extension;   // e.g. wav, flac, etc.
    std::string filePath;    // Relative path with respect to instrument definition file or some
                             // global sample content directory.

    int  numFrames   = 0;
    int  numChannels = 0;
    TPar sampleRate  = TPar(-1);  // in Hz, -1 is code for unknown
    TPar rootKey     = TPar(-1);  // as MIDI key in 0..127, not necessarily integer

    // Maybe add: peakAmplitude, rmsAmplitude, etc., maybe Category (harmonic, inharmonic, noisy, 
    // transient, speech, music, ...)
  };
  */

  /** A class to represent various additional (and optional) playback settings of a region, group 
  or instrument. Such additional settings include additional constraints for the circumstances 
  under which a particular sample should be played. Key- and velocity ranges are the obvious 
  primary constraints (and they therefore are directly baked into the Region class below), but sfz
  defines many more. But the settings doesn't need to be playback constraints - that's only one 
  type of setting. Other types are things like envelope settings, filter frequencies, etc. */
  class PlaybackSetting
  {

    enum class Type
    {
      ControllerRangeLo, ControllerRangeHi, PitchWheelRange,  // 

      AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,

      FilterCutoff, FilterResonance, FilterType

      //...tbc...
    };

  private:

    float value = 0.f;
    // hmm - it seems, for the controllers, we need 2 values: controller number and value - but it
    // would be wasteful to store two values for all other settings as well...hmmm...maybe 
    // groups/regions need to maintain 2 arrays with settings, 1 for the 1-valued settings and 
    // another for the 2-valued settings - maybe have classes PlaybackSetting, PlaybackSetting2Val

  };

  /** A group organizes a bunch of regions into a single entity for which performance settings can 
  be set up which will be applicable in cases where the region does not itself define these 
  settings, so they act as fallback values. */
  class Group
  {

  public:

    int addRegion();

    // todo: removeRegion, etc.

    const Region* getRegion(int i)
    {
      if(i < 0 || i >= (int)regions.size()) {
        rsError("Invalid region index");
        return nullptr; 
      }
      return &regions[i];
    }
    // for some reason, i get compiler errors when trying to put this into the cpp file 
    // -> figure out

  private:

    std::vector<Region> regions;
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

  private:

    AudioFileStream* sample = nullptr;  
    Group* group = nullptr;             // pointer to the group to which this region belongs

    char loKey = 0, hiKey = 127;
    char loVel = 0, hiVel = 127;
    // todo: maybe package loKey/hiKey, loVel/hiVel into a single char to save memory


    std::vector<PlaybackSetting> settings;
    // for more restrictions (optional) restrictions - sfz can restrict the playback of samples
    // also based on other state variables such as the last received controller of some number,
    // last received pitchwheel, etc. ...but maybe a subclass RestrictedRegion should be used
    // for that - i don't think, it will be used a lot and will just eat up memory when it's
    // present in the baseclass...or maybe it should have a more general array of RegionFeatures
    // which may also include loop-settings and the like

    //std::string name;

    friend class Group;
  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Return codes for the setup functions. We use encodings as negative integers so we can use 
  them also for functions which use positive integers as valid return values. */
  enum ReturnCode
  {
    success      = -1,  //< Operation completed successfully. 
    nothingToDo  = -2,  //< There was nothing to actually do. State was already as desired.
    memAllocFail = -3,  //< Memory allocation failure.
    invalidIndex = -4   //< An invalid index was passed.
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

  // todo:
  // int setRegionSample(int groupIndex, int regionIndex, int sampleIndex); 
  // setRegionLoKey, setRegionHiKey, setRegionLoVel, setRegionHiVel


  // todo: addGroup, addRegion(int group, ..), removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, addSampleToPool, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,



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
  // Maybe it should be non-const - but no: the caller should not be able to change the loKey/hiKey
  // settings because that would require a change to the regionsForKey array...but maybe we should 
  // get rid of that anyway - it might be a pointless attempt to optimization -> benchmark!
  // for some reason, i get compiler errors when trying to put this into the cpp file 
  // -> figure out

  // getGroup, getRegion, getStateAsSFZ, isSampleInPool, getNumGroups, getNumRegionsInGroup(int)
  // getNumRegions()


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

  class AudioFileStream : public AudioStream
  {

    // ToDo: add == operator based on fileName, extension, path (and maybe (meta)data such as 
    // sampleRate, numChannels, numFrames, ...)

    // inquiry: existsOnDisk (idea: we should allow to work with samples that are not stored on 
    // disk but rather pre-rendered programmatically into RAM at runtime)

    //virtual void getFrame(int sampleIndex, TSmp* destination) {}

    // todo: getBlock


  protected:

    std::string fileName;  // without filename extension
    std::string extension; // filename extension
    std::string path;      // relative path from instrument definition file (e.g. Piano.sfz)

    // Maybe we should also store the root-directory to which the path is relative, but maybe just
    // as an integer that selects between various pre-defined root-directories that should exist at
    // the rsSamplerEngine level. For example: 0: instrument directory, 1: factory sample 
    // directory, 2: user sample directory, etc. (but not hardcoded - meanings of the indices 
    // should be flexible). This makes it more reasonably possible to uniquely identify samples. 
    // It's totally possible to have samples in an instrument with same relative paths and 
    // filenames but with respect to different root directories. Yes - that would be weird,
    // but the engine should neverless be able to handle such situations.
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


    void clear();


    // todo:
    // setup:   removeSample
    // inquiry: hasSample


  protected:

    std::vector<const AudioFileStream*> samples;

  };




  /** Defines a set of regions. Used to handle note-on/off events efficiently. Not to be confused 
  with groups. This class exists for purely technical reasons (i.e. implementation details) and 
  does not map to any user concept. */
  class RegionSet
  {

  public:

    void addRegion(const Region* r) { regions.push_back(r); }
    // todo: removeRegion, containsRegion

  private:

    std::vector<const Region*> regions; // pointers to the regions belonging to this set
  };


  //-----------------------------------------------------------------------------------------------
  // \name Internal functions

  void addRegionForKey(uchar k, const Region* region);


  /** Returns true, iff the given region should play when the given key is pressed with given 
  velocity. This will also take into account other playback constraints defined for the region 
  and/or its enclosing group. */
  bool shouldRegionPlay(const Region* r, const char key, const char vel);



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

  int numChannels = 2;
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



#endif