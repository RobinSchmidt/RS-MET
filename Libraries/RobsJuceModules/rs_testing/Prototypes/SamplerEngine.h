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
features, but sfz is the spec after which this engine is roughly modeled. */

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

  // getGroup, getRegion, getStateAsSFZ


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(TSig* frame, int numChannels);

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

  };

  class AudioFileStream : public AudioStream
  {

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


  class Region;

  class Group
  {


  private:

    std::vector<Region*> regions; // pointers to the regions belonging to this group

    // may be add these later:
    //std::string name;  
  };


  /** A class to represent various additional (and optional) features of a region, such as 
  additional constraints for the circumstances under which a particular sample should be played. 
  Key- and velocity ranges are the obvious primary constraints (and they therefore are directly 
  baked into the region class below), but sfz defines many more. But the "feature" doesn't need to 
  be a playback constraint - that's only one type of feature. Other types are things like envelope 
  settings, filter frequencies, etc. */
  class RegionFeature
  {

    enum class Type
    {
      ControllerRange,
      PitchWheelRange
      //...tbc...
    };

  private:

  };
  // maybe rename - the same features could be applied to groups and the whole instrument, too

  class Region
  {

  private:

    AudioFileStream* sample = nullptr;  
    Group* group = nullptr;             // pointer to the group to which this region belongs

    char loKey = 0, hiKey = 127;
    char loVel = 0, hiVel = 127;
    // todo: maybe package loKey/hiKey, loVel/hiVel into a single char to save memory


    // may be add these later:
    //std::vector<Restriction> restrictions; 
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
  // \name Data
 
  RegionSet regionsForKey[128];
  // For each key, we store a set of regions that may need to be played, when the key is pressed.
  // If they actually need to be played indeed is determined by other constraints as well, such as 
  // velocity, last received controller and/or pitch-wheel values, etc.

};



#endif