#ifndef RAPT_SAMPLERENGINE_H
#define RAPT_SAMPLERENGINE_H

//=================================================================================================

/** A class for representing musical events such as note-on/off etc. Think of it as a class to 
represent MIDI events but with some of its anachronistic restrictions lifted, such as the abysmal 
resolution of values. The template parameter T is supposed to be either float or double. For easy 
conversion and compatibility with MIDI, we still follow the (now historical) convention that 
values are in the range 0..127, but now with much higher resolution due to the floating point 
representation. So, in a nutshell, this is a class for MIDI events but with higher resolution for 
all the values. */

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
will be used.

*/

template<class TSig, class TPar, class TSmp> 
// TSig: type for signals during processing (typically float, double or maybe a SIMD type)
// TPar: type for continuous numeric parameters (typically float or double)
// TSmp: type of sample values as they are stored in memory (typically float)

class rsSamplerEngine
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  // todo: addGroup, addRegion(int group, ..), removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, addSampleToPool, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  // getGroup, getRegion, getStateAsSFZ, isSampleInPool, 


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(TSig* frame, int numChannels);
  // maybe numChannles should be a member

  void processBlock(TSig** block, int numFrames, int numChannels);

  void handleMusicalEvent(const rsMusicalEvent<TPar>& ev);

  // void processFrameVoice, processBlockVoice


protected:

  //-----------------------------------------------------------------------------------------------
  // \name Helper Classes (may be factored out at some point)

  class AudioStream
  {

  public:

    /** For random access. Writes the sample frame with given index into the given destination. */
    virtual void getSampleFrame(int sampleIndex, TSmp* destination) = 0;

    /** Subclasses may want to override this for optimizing the blockwise access. */
    virtual void getBlock(int startIndex, int length, TSmp** destination)
    {
      for(int i = 0; i < length; i++)
        getSampleFrame(startIndex + i, destination[i]);
        // assumes interleaved storage in memory
    }

  protected:

    int sampleRate  = 44100;
    int numChannels = 1;
    int numFrames   = 0;         // maybe use -1 to encode "unknown"? would that be useful?

  };

  class AudioFileStream : public AudioStream
  {

    // ToDo: add == operator based on fileName, extension, path (and maybe (meta)data such as 
    // sampleRate, numChannels, numFrames, ...)

    // inquiry: existsOnDisk (idea: we should allow to work with samples that are not stored on 
    // disk but rather pre-rendered programmatically into RAM at runtime)

  protected:

    std::string fileName;  // without filename extension
    std::string extension; // filename extension
    std::string path;      // relative path from instrument definition file (e.g. Piano.sfz)
  };

  class AudioFileStreamPreloaded
  {

    void getSampleFrame(int sampleIndex, TSmp* destination) override
    {
      // ...something to do...
    }

  protected:

    TSmp** data;           // pointer to the sample data

  };


  class SamplePool
  {

    // todo: 
    // setup:   addSample, removeSample, clearAllSamples
    // inquiry: hasSample


  protected:

    std::vector<AudioFileStream*> samples;

  };


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
      ControllerRange, PitchWheelRange,

      AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,

      FilterCutoff, FilterResoance, FilterType

      //...tbc...
    };

  private:

  };





  class Region;

  class Group
  {


  private:

    std::vector<Region*> regions; 
    /**< Pointers to the regions belonging to this group. */

    std::vector<PlaybackSetting> settings;
    /**< Settings that apply to all regions within this group, unless a region overrides them with
    its own value for a particular setting. */


    // may be add these later:
    //std::string name;  
  };





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
  };

  /** Defines a set of regions. Used to handle note-on/off events efficiently. Not to be confused 
  with groups. This class exists for purely technical reasons (i.e. implementation details) and 
  does not map to any user concept. */
  class RegionSet
  {

    std::vector<Region*> regions; // pointers to the regions belonging to this set
  };


  //-----------------------------------------------------------------------------------------------
  // \name Internal functions


  /** Returns true, iff the given region should play when the given key is pressed with given 
  velocity. This will also take into account other playback constraints defined for the region 
  and/or its enclosing group. */
  bool shouldRegionPlay(const Region* r, const char key, const char vel);



  //-----------------------------------------------------------------------------------------------
  // \name Data
 
  RegionSet regionsForKey[128];
  /**< For each key, we store a set of regions that may need to be played, when the key is pressed.
  If they actually need to be played indeed is determined by other constraints as well, such as 
  velocity, last received controller and/or pitch-wheel values, etc. */

  SamplePool samplePool;
  /**< The pool of samples that are in use for the currently loaded instrument. The samples are 
  pooled to avoid redundant storage in memory when multiple regions use the same sample. */

  std::vector<Group> groups;
  /**< The groups contained in this instrument. Each group may contain set of regions. */

  std::vector<PlaybackSetting> settings;
  /**< Playback settings that apply to all groups within this instrument, unless a group (or 
  region) overrides a setting with its own value. **/


};



#endif