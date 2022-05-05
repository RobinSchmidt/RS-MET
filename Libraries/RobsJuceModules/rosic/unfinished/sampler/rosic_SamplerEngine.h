#ifndef rosic_SamplerEngine_h
#define rosic_SamplerEngine_h

namespace rosic {
namespace Sampler {

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
will be used, which means the corresponding feature is not used, i.e. bypassed.  

SFZ-Links:
  https://sfzformat.com/  https://sfzformat.com/opcodes/?v=1  https://sfzformat.com/legacy/
  https://github.com/sfzformat
  https://github.com/sfzformat/sfzformat.github.io
  https://github.com/sfz  https://github.com/sfz/opcode-suggestions/issues
  https://www.kvraudio.com/forum/viewtopic.php?f=42&t=508861

*/

class rsSamplerEngine
{

public:


  // For convenience:
  using uchar = unsigned char;
  using Region = SfzInstrument::Region; // todo: make a subclass here that adds the stream field
  using Group  = SfzInstrument::Group;

  //-----------------------------------------------------------------------------------------------
  // \name Internal Helper Classes

  /** A struct that can be returned from midi event handling functions to inform the caller, how
  the event has changed the playback status of the engine. For example, a noteOn event will
  typically result in the start of playback for one or more layers/regions. This will be reflected
  in the numLayersStarted field. */
  struct PlayStatusChange
  {
    int numLayersStarted = 0;
    int numLayersStopped = 0;
  };
  // ToDo: have also fields for numProcessorsAdded/Removed, numModulatorsAdded/Removed, etc.


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  rsSamplerEngine(int maxNumLayers = 16);
  // todo: use some higher default value - what is reasonable here needs some testing in realistic
  // scenarios

  virtual ~rsSamplerEngine();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the maximum number of layers/regions that can be played simultaneously. */
  virtual void setMaxNumLayers(int newMax);

  /** Clears the sfz instrument definition and the samplePool. */
  void clearInstrument();

  void clearAllSfzSettings() { sfz.clearAllSettings(); }
  // what's the differece to clearInstrument? seems like clearInstrument() also clears the 
  // samplePool? -> document this. maybe try to get rid of one of them - it's a bit confusing to 
  // have both. does it also clear the sample opcode? i don't think so. document this. maybe rename
  // the function to something that conveys this information


  /** Adds a new sample to our pool of samples. After the sample has been added, regions can be
  defined that make use of it. */
  int addSampleToPool(float** data, int numFrames, int numChannels, double sampleRate,
    const std::string& path);
  // Maybe rename to addSample, it should return the index of the sample in the sample-pool
  // maybe make a struct SampleMetaData containing: numFrames, numChannels, sampleRate, rootKey
  // todo: take reference to a metaData object
  // make const-correct

  /** Conveninece function for mono samples stored in a std::vector. */
  int addSampleToPool(const std::vector<float>& data, double sampleRate, const std::string& path);

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
  int setRegionSetting(int groupIndex, int regionIndex, Opcode type, float value, int index)
  { return sfz.setRegionSetting(groupIndex, regionIndex, type, value, index); }

  int setGroupSetting(int groupIndex, Opcode type, float value, int index)
  { return sfz.setGroupSetting(groupIndex, type, value, index); }

  int setInstrumentSetting(Opcode type, float value, int index)
  { return sfz.setInstrumentSetting(type, value, index); }

  // under construction:
  int setRegionModulation(int groupIndex, int regionIdex, OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode)
  {
    return sfz.setRegionModulation(groupIndex, regionIdex, 
      modSrcType, modSrcIndex, modTarget, modTargetIndex, modDepth, modMode);
  }

  int setGroupModulation(int groupIndex, OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode)
  {
    return sfz.setGroupModulation(groupIndex, 
      modSrcType, modSrcIndex, modTarget, modTargetIndex, modDepth, modMode);
  }

  int setInstrumentModulation(OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode)
  {
    return sfz.setInstrumentModulation(
      modSrcType, modSrcIndex, modTarget, modTargetIndex, modDepth, modMode);
  }


  // ToDo: These functions are currently used only in test code. I think, they are not yet ready
  // for production because for certain removals, we may have to take additional actions such as 
  // updating our regionsForKey array

  /** NOT YET READY FOR PRODUCTION. */
  rsReturnCode removeRegionSetting(int groupIndex, int regionIdex, Opcode type, int index)
  {
    return sfz.removeRegionSetting(groupIndex, regionIdex, type, index);
  }
  /** NOT YET READY FOR PRODUCTION. */
  rsReturnCode removeGroupSetting(int groupIndex, Opcode type, int index)
  {
    return sfz.removeGroupSetting(groupIndex, type, index);
  }
  /** NOT YET READY FOR PRODUCTION. */
  rsReturnCode removeInstrumentSetting(Opcode type, int index)
  {
    return sfz.removeInstrumentSetting(type, index);
  }

  rsReturnCode removeRegionModulation(int groupIndex, int regionIdex, OpcodeType modSrcType, 
    int modSrcIndex,  Opcode modTarget, int modTargetIndex)
  {
    return sfz.removeRegionModulation(groupIndex, regionIdex, 
      modSrcType, modSrcIndex, modTarget, modTargetIndex);
  }
  // Maybe instead of a return code, it should return the number of removed connections? But no:
  // that's supposed to be 1 or 0 anyway because we don't allow multiple connections between the 
  // same pair of pins. I think, that would be technically possible though, but it doesn't seem 
  // to make much sense from a user perspective.
  // todo: removeRegionModulation, removeInstrumentModulation...but maybe these function will not
  // be needed in the context of loading sfz files...hmm...but on an interactive GUI, we probably
  // want to be able to remove modulations one by one



  /** NOT YET READY FOR PRODUCTION. */
  rsReturnCode clearRegionSettings(int groupIndex, int regionIdex)
  {
    return sfz.clearRegionSettings(groupIndex, regionIdex);
  }









  // todo: setGroupSetting, setInstrumentSetting, removeRegion/Group, clearGroup, clearRegion, 
  // clearInstrument, removeSampleFromPool, replaceSampleInPool, setupFromSFZ,

  /** Sets up the engine from the given sfz data object and returns rsReturnCode::success, if all
  is well or...  */
  int setupFromSFZ(const SfzInstrument& sfz);

  /** Sets up the root directory, with respect to which sfz file paths are interpreted in saveToSfz
  and loadFromSfz. This is supposed to be an absolute path. Returns false in case, the directory
  doesn't exist. */
  bool setSfzRootDir(const char* path);

  /** Writes the current instrument definition into an sfz file with given path which can be given
  as absolute path or relative with respect what was set via setSfzRootDir. */
  bool saveToSFZ(const char* path, bool pathIsAbsolute = false) const;
  // -don't use all caps in SFZ
  // -return a return-code instead of bool

  /** Loads the instrument definition given by an sfz file with the given path which can be given
  as absolute path or relative with respect what was set via setSfzRootDir. Returns
  rsReturnCode::success if all wen well or rsReturnCode::fileLoadError if loading of the sfz or any
  of the used samples has failed. */
  int loadFromSFZ(const char* path, bool pathIsAbsolute = false);
  // ToDo: In case of failure, maybe return a more specific error code/object, indicating, which 
  // file(s) exactly failed to load. This is an information that may be eventually displayed to the
  // user on a GUI.

  /** Sets the engine up according to a string representing the contents of an sfz file. It doesn't
  matter whether or not such an sfz file actually exists. The function is mostly meant for unit
  testing the sfz-parser with programmatically created sfz-data which doesn't necessarily have to
  come from a file on disk. */
  int setFromSFZ(const std::string& sfzFileContents);



  /** Sets the sample-rate, at which this engine should operate. This change will affect only
  RegionPlayer objects that were started after calling this function. It's supposed to be called in
  a suspended state anyway, not in the middle of the processing. */
  void setSampleRate(double newRate) // use float
  { 
    playStatus.sampleRate = (float) newRate;
    //sampleRate = newRate;   // redundant -> remove
  }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of regions in the group with given groupIndex. */
  int getNumRegions(int groupIndex) const { return sfz.getNumRegions(groupIndex); }
  // todo: maybe assert the groupIndex is valid - if not, return invalidIndex, maybe rename to 
  // getNumRegionsInGroup, getNumRegions should have no parameter and return the total number

  /** Returns const a pointer to the region object with the given group- and region index or a 
  nullptr if the combination of indices is invalid. The client can use the pointer to inquire the
  region in some more detail, perhaps for GUI purposes, but it can't edit the region's state 
  because that should never happen behind the back of the sampler engine because the engine may 
  have to take additional actions when certain aspects of a region change. If the client wants to 
  edit the region, it needs to use the appropriate region-editing functions of the rsSamplerEngine 
  object. */
  const Region* getRegion(int gi, int ri) const  { return sfz.getRegion(gi, ri); }

  /** Returns the number of groups in the instrument. */
  int getNumGroups() const { return sfz.getNumGroups(); }

  /** Returns a pointer to the index-th group or a nullptr if the index is invalid. */ 
  const Group* getGroup(int index) { return sfz.global.getGroup(index); }

  /** Returns the number of regions in the instrument definition that use the sample with the given
  path. The path is interpreted to be relative to the sample-path which is by default the directory
  where the sfz file resides. */
  int getNumRegionsUsing(const std::string& samplePath) const;

  /** Returns the number of regions in the instrument definition that use the sample with the given
  index in our samplePool or rsReturnCode::invalidIndex, if the given sampleIndex is invalid. */
  int getNumRegionsUsing(int sampleIndex) const;
  // ToDo: move to protected - client code has nothing to do with the innards of our sample pool


  // getStateAsSFZ, getSampleIndex(const string& uniqueName) ..or maybe it should take a pointer
  // to a SampleMetaData object

  /** Returns true, iff the given group index is valid, i.e. >= 0 and < numGroups. */
  bool isGroupIndexValid(int i) const { return sfz.isGroupIndexValid(i); }

  /** Returns true, iff the given pair of group- and region index is valid, i.e. a region with this
  pair of indices actually exists in the current instrument definition. */
  bool isIndexPairValid(int groupIndex, int regionIndex) const
  {
    return sfz.isIndexPairValid(groupIndex, regionIndex);
  }

  /** Returns true, iff the given sample index is valid, i.e. a sample with this index actually
  exists our sample pool. */
  bool isSampleIndexValid(int sampleIndex) const
  {
    return samplePool.isSampleIndexValid(sampleIndex);
  }

  /** Returns the index of the sample represented by the given string in our sample pool or -1, if
  the sample is not in the pool. */
  int findSampleIndexInPool(const std::string& sample) const;

  /** Returns true, if the sample represented by the given string (as relative path with respect to
  some root directory) is present in our samplePool. */
  bool isSampleInPool(const std::string& sample) const
  {
    return findSampleIndexInPool(sample) != -1;
  }

  /** Returns the maximum number of layers that can play simultaneously. */
  int getMaxNumLayers() const { return (int)playerPool.size(); }

  /** Returns the number of currently playing layers. */
  int getNumActiveLayers() const { return (int)activePlayers.size(); }
  // todo: rename to getNumActiveRegions or getNumPlayingRegions

  /** Returns the number of layers that are currently not playing, i.e. still available for adding
  a new layer to the playback. */
  int getNumIdleLayers() const { return (int)idlePlayers.size(); }

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

  /** Returns a const pointer to the SfzInstrument object that represents the current instrument
  settings. */
  const SfzInstrument& getInstrumentData() const { return sfz; }

  /** Given a path which can be either relative to our sfzDir or absolute, this function returns
  the corresponding absolute path as std::string. That means, if pathIsAbsolute is true, it just
  converts the given char-array to a std::string as is and pathIsAbsolute is false, it assumes that
  the given path is relative and prepends the sfzDir in the returned string */
  std::string getAbsolutePath(const char* path, bool pathIsAbsolute = false) const;

  /** Returns the sample rate at which this engine runs, i.e. produces its output. */
  float getOutputSampleRate() const { return playStatus.sampleRate; }

  //-----------------------------------------------------------------------------------------------
  // \name Processing



  //virtual void processFrame(double* left, double* right);

  /** Produces one stereo sample frame at a time. ...tbc... */
  virtual void processFrame(float* left, float* right);

  /** Produces one block of stereo samples at a time. Note that the sampler's output will be added
  to whatever audio data is already present in the block. The rationale behind this is that in the
  context of ToolChain or on a DAW track, modules or plugins are supposed to be in a series 
  connection but for instrument plugins, such a series connection would mean that we always only 
  hear the last instrument in the chain whereas the desired behavior is most probably layering.
  If you don't like that behavior, you need to clear the buffer before calling processBlock.  */
  virtual void processBlock(float** block, int numFrames);

  // void processFrameVoice, processBlockVoice

  /** Stops the playback of all currently active RegionPlayers immediately. This is a rather hard
  reset which may be appropriate to call when a midi reset message is received or before loading a
  new patch. It returns the number of players that were affected, i.e. the number of players that
  were in active state before the call. */
  int stopAllRegionPlayers();
  // todo: return a PlayStatusChange

  /** Calls stopAllRegionPlayers. Function is for consistency with the rest of the library. */
  virtual void reset() { stopAllRegionPlayers(); }

  /** Handles a musical (i.e. midi) event. This will typically change the playback status of the
  engine, for example by triggering the playback of one or more layers/region. If this status
  change is relevant for the caller, it can inspect the returned PlayStatusChange after the call,
  otherwise, the caller may ignore it. */
  PlayStatusChange handleMusicalEvent(const rsMusicalEvent<float>& ev);

  //===============================================================================================

protected:

  /** Returns a pointer to a region that allows modification of the region settings. */
  Region* getRegionMutable(int gi, int ri) const  { return sfz.getRegionMutable(gi, ri); }


  /** Defines a set of regions. Used to handle note-on/off events efficiently. Not to be confused
  with groups. This class exists for purely technical reasons (i.e. implementation details) and
  does not map to any user concept. It holds pointers to all regions that should be played when
  a particluar key is pressed so we can quickly find them. */
  class RegionSet
  {

  public:


    int findRegionIndex(const Region* r) const
    {
      if(regions.empty()) return -1;
      return RAPT::rsArrayTools::findIndexOf(&regions[0], r, (int)regions.size());
    }

    bool containsRegion(const Region* r) const { return findRegionIndex(r) != -1; }

    void addRegion(const Region* r)
    {
      if(containsRegion(r)) {
        RAPT::rsError("Don't add regions twice!"); return;
      }
      regions.push_back(r);
    }

    bool removeRegion(int i)
    {
      if(i < 0 || i >= (int)regions.size())  {
        RAPT::rsError("Invalid region index"); return false;
      }
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
  bool isSampleUsedIn(const AudioFileStream<float>* sample, const SfzInstrument& sfz);

  /** Stops the player at the given "activeIndex" which is the index into our "activePlayers"
  array. This results in the removal of the player from "activePlayers" and adding it back to
  "idlePlayers". The return value is either rsReturnCode::success or rsReturnCode::invalidIndex, if
  the activeIndex was not a valid index into our activePlayers array. */
  virtual rsReturnCode stopRegionPlayer(int activeIndex);
  // maybe rename to stopLayer

  /** Stops the given number of most recently started layers (i.e. RegionPlayers). This is used
  to roll back RegionPlayers that have been started in noteOn when it is later discovered that
  some more RegionPlayers that would be needed for the same not could not be made available due
  to lack of resources like RegionPlayers themselves, required DSP objects, etc. In such a case
  we need to roll back the starts of all players started so far for the note. We don't want 
  notes to play with an incomplete set of layers. It's all or nothing - either all lyers for the
  given note will play or none. We either do it right or not at all! */
  virtual void stopMostRecentLayers(int number);

  /** Returns the AudioFileStream object that is used to stream the actual sample data for the
  given region. A pointer to this object is supposed to be stored within the region object
  ...tbc... */
  static const AudioFileStream<float>* getSampleStreamFor(const Region* r);

  /** Handles a noteOn event with given key and velocity and returns

  obsolete:
  either rsReturnCode::success,
  if we had enough voices available to serve the request or rsReturnCode::voiceOverload, in case
  the noteOn could not be handled due to inavailability of a sufficient number of idle voices. If
  no sufficient number of idle voices was available and the noteOn should actually have triggered
  playback of multiple samples, none of them will be triggered. It's an all-or-nothing thing: we
  don't ever trigger playback for only a subset of samples for a given noteOn. */
  virtual PlayStatusChange handleNoteOn(uchar key, uchar vel);

  /** Analogous to handleNoteOn. It may also return rsReturnCode::voiceOverload in cases where the
  noteOff is supposed to trigger relase-samples. In such a case, none of the release-samples will
  be triggered. */
  virtual PlayStatusChange handleNoteOff(uchar key, uchar vel);


  /** Removes those samples from our sample pool that are not used in the given sfz instrument
  specification. Returns the number of samples that were removed. */
  int removeSamplesNotUsedIn(const SfzInstrument& sfz);
  // maybe rename to removeUnusedSamples. But that name is more ambiguous: it could be interpreted
  // as "unused in the current sfz member", so maybe don't

  /** Adds all samples to our sample pool that are used in the given sfz instrument definition, if
  they are not already there. Returns the number of samples that were added or
  rsReturnCode::fileLoadError if any of the files failed to load. */
  int addSamplesUsedIn(const SfzInstrument& sfz);
  // maybe rename to loadSamples or loadSamplesFor

  /** Sets up all the AudioStream pointers in all the regions in our sfz member. */
  int setupAudioStreams();

  /** Sets up ourregionsForKey array according to our sfz member. */
  void setupRegionsForKey();

  /** Pre-allocates the memory needed for the DSP objects (filters, equalizers, etc). It's virtual
  so our subclass rsSamplerEngine2 can override, because it may have different requirement for how 
  much to allocate and where because it has a different architecture of the DSP chain. */
  virtual void preAllocateDspMemory();


  //-----------------------------------------------------------------------------------------------
  // \name Data

  SfzInstrument sfz;
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

  DspResourcePool dspPool;
  /**< The pool of DSP processors, modulators and modulation connection. When a region starts 
  playing, the respective RegionPlayer object may grab objects from this pool such as filters, 
  waveshapers, envelopes, etc. - whatever it needs to play the region appropriately. */

  //SignalProcessorPool processorPool;
  /**< The pool of DSP processors. When a region starts playing, the respective RegionPlayer object
  may grab objects from this pool such as filters, waveshapers, etc. - whatever it needs to play 
  the region appropriately. */

  //std::vector<Processor*> masterEffects; 
  // not yet used..their oscs etc should probably retriggered only when there's currently no note
  // playing *and* the output is silent


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

  PlayStatus playStatus;
  /**< Intermediate variables used for the computation of things like per-sample increment, final 
  amplitude etc. according to key, vel, keytrack, veltrack, tune, transpose, etc. We need to pass 
  around such a struct from noteOn to have a place into which we can accumulate all the modifiers
  such that the RegionPlayer itself needs ot store only the final values. */

  //double sampleRate = 44100.0; // redundant with PlayStatus.sampleRate..obsolete
  /**< Sample rate at which this object runs. */

  // Some info that can be inquired from client code after loading a new sfz file:
  int numSamplesRemoved = 0;  /**< Number of samples that were unloaded. */
  int numSamplesLoaded  = 0;  /**< Number of samples that were loaded. */
  int numSamplesFailed  = 0;  /**< Number of samples that failed to load. */

  // These are not yet used - currently, both are assumed to be the project directory (at least
  // for the unit tests):
  std::string sfzDir;         /**< Root directory for .sfz files */
  //std::string wavDir;         /**< Root directory for .wav files */



  // Some member variables for features whose implementation is actually deferred to the subclass
  // rsSamplerEngine2. For convenience of implementation, we already need them to be members of 
  // this baseclass here. In this class here, they are always false and we don't even have setters 
  // for them. Having them as baseclass members avoids the necessity to pass around an excessive 
  // amount of parameters to certain member functions. Their existence here is a provision to 
  // enable features in the subclass.
  //
  // Flags to decide if the group- and/or instrument settings and/or modulations should be applied 
  // on top of the region settings/modulations. This is a feature not present in the sfz spec. 
  //bool groupSettingsOverride  = true;
  //bool regionSettingsOverride = true;
  bool busMode = false;

  // This object is used for various translation tasks:
  //SfzCodeBook sfzTranslator;

  static int instanceCounter;
  /** An instance counter which is needed to clean up the SfzCodeBook singleton object when the
  last instance of the sampler engine is deleted. */


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

class rsSamplerEngine2 : public rsSamplerEngine // maybe rename to SamplerMixEngine
{

public:

  // for convenience:
  using uchar = unsigned char;
  //using Region = SfzInstrument::Region; // todo: make a subclass here that adds the stream field
  //using Group  = SfzInstrument::Group;
  //using PlaybackSetting = SfzInstrument::PlaybackSetting;


  rsSamplerEngine2(int maxNumLayers = 16);


  void setMaxNumLayers(int newMax) override;


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
  //void setGroupModulationsOnTop(bool onTop) { groupModulationsOnTop = onTop; }
  // deprecated ...the comment may still be relevant for implementing the override vs accumulate
  // modes for modulations, so it has not yet been deleted

  /** Decides whether sfz opcodes work accumulatively or in override mode where the latter is the
  default, as specified by the sfz spec and the former is another kind of behavior which is also
  sometimes useful. "Override", in this context, means that a region opcode, if present, will
  override the corresponding opcode of the enclosing group and/or instrument. Likewise, a group
  opcode, if present, will  override the corresponding opcode of the instrument. "Accumulate", on
  the other hand, means that group opcodes will be applied on top of instrument opcodes and region
  opcodes will be applied on top of group (and/or instrument) opcodes. For example, if a group
  specifies a volume of -6dB and a region within the group specifies a volume of -3dB, the
  resulting volume of the region in playback will be -9dB in accumulative mode and -3dB in override
  mode. -3dB is what you should get as per the sfz spec [verify!]. */
  /*
  void setOpcodesAccumulate(bool shouldAccumulate)
  {
    regionSettingsOverride = !shouldAccumulate;
    groupSettingsOverride  = !shouldAccumulate;
    reset();  // changing this setting is disruptive - we need a reset
  }
  */
  // maybe call it setDrumBusMode, maybe call the whole class rsDrumSamplerEngine. It's typical for
  // Drum patches that each sample is mapped to a particular key and then the drums are combined
  // to a drum bus to which more FX are applied etc. ...maybe call it groupsAreBusses or something
  // like that

  /** The availability of a bus-mode is the main feature that distinguishes this subclass from the
  rsSamplerEngine baseclass. The baseclass provides behavior according to the sfz specification 
  where opcodes defined in an enclosing group and/or instrument define fallback values for region 
  opcodes. In bus-mode, these group- or instrument-wide settings are not interpreted as fallback 
  values but instead as values that should be applied additionally on top of the region settings
  after regions within a group have been mixed together (to form a sub-bus) and after the group
  sub-busses have been mixed together to form the instrument's master-bus. This function switches
  to this optional other kind of global behavior. */
  void setBusMode(bool groupsAreSubBussesInstrumentIsMasterBus) 
  { 
    reset();  // changing this setting is disruptive - we need a reset
    busMode = groupsAreSubBussesInstrumentIsMasterBus; 
  }


  /** Returns the number of currently playing group players. */
  int getNumActiveGroupPlayers() const { return (int)activeGroupPlayers.size(); }


  void processFrame(float* left, float* right) override;

  //void processFrame(float* left, float* right) override;

  //void processBlock(float** block, int numFrames) override;

  // void processFrameVoice, processBlockVoice

  /** Stops the playback of all currently active RegionPlayers, GroupPlayers and the 
  instrumPlayer and takes care of rolling back all acquired resources inot idle state. The 
  returned integer is the number of stopped region players. */
  int stopAllPlayers();

  virtual void reset() override { stopAllPlayers(); }


protected:


  PlayStatusChange handleNoteOn(uchar key, uchar vel) override;

  PlayStatusChange handleNoteOff(uchar key, uchar vel) override;

  rsReturnCode stopRegionPlayer(int activeIndex) override;

  /** Returns true, iff the settings are such that the additional features of this subclass are not
  used so we can fall back to the (simpler) baseclass implementations in our overrides. */
  bool canFallBackToBaseclass() const { return !busMode; }

  /** Updates our active/idleGroupPlayer arrays according to a status change in the
  active/idleRegionPlayer arrays...tbc... */
  void updateGroupPlayers(PlayStatusChange psc, uchar key, uchar vel);

  /** Returns the index within our activeGroupPlayers array at which the group player for the given
  group is located or -1 if there is no currently active player for the given group. */
  int getActiveGroupPlayerIndexFor(const SfzInstrument::Group* group);

  /** This starts a new group player for the given region player. This is used in e.g. handleNoteOn
  when region players were triggered for which we do not already have an active group player in
  use. Return true, if the player could be started. If it could not be started, the caller may 
  need to take care of some cleanup (roll back the regionPlayer). */
  bool startGroupPlayerFor(RegionPlayer* regionPlayer, uchar key, uchar vel);

  /** Stops the groupPlayer with the given activeIndex, i.e. moves it from the activeGroupPlayers
  array to the idleGroupPlayers array. */
  int stopGroupPlayer(int activeIndex);

  // under construction
  bool startInstrumPlayerFor(RegionPlayer* regionPlayer, uchar key, uchar vel); 
  void stopInstrumPlayer();



  // under construction:
  InstrumPlayer instrumPlayer;  // we only need one of them
  std::vector<GroupPlayer*> activeGroupPlayers;
  std::vector<GroupPlayer*> idleGroupPlayers;
  std::vector<GroupPlayer>  groupPlayerPool;
  // I think, we need as many GroupPlayers as there are RegionPlayers because in the worst case,
  // each region could be within its own group...although that's probably really uncommon

};

//=================================================================================================

/** Subclass that contains some extra functions that facilitate testing which should not go into
the production code.

ToDo: maybe move into the test project or the rs_testing juce module */

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
  //static int getRegionPlayerSize() { return sizeof(RegionPlayer); }


  void preAllocateDspMemory() override
  {
    rsSamplerEngine::preAllocateDspMemory() ;
  }
  // Function is protected in baseclass because normally, client code should not need to deal with
  // it but in some unit tests, we call it manually so we un-protect it here by overriding it with
  // a public function that just calls the basclass method.

  PlayStatusChange handleNoteOn(uchar key, uchar vel) override
  {
    return rsSamplerEngine::handleNoteOn(key, vel);
    // The purpose of this override is to make it public. It's protected in the baseclass because
    // in production, client code is supposed to use the more general midi-event handler. But in
    // tests, it's sometimes more convenient to be able to call this directly.
  }


};

//=================================================================================================



class rsSamplerEngine2Test : public rsSamplerEngine2
{

public:

  using rsSamplerEngine2::rsSamplerEngine2;  // inherit constructors

  void setMaxNumFilters(int newMax) 
  { 
    reset();
    dspPool.setMaxNumFilters(newMax); 
  }

  // move into baseclass - jura::Sampler will need it:
  void setMaxNumLayers(int newMax) 
  {
    reset();
    playerPool.resize(newMax); // maybe rename to regionPlayerPool
    for(auto& p : playerPool) {
      p.setDspResourcePool(&dspPool);
      p.setPlayStatusPointer(&playStatus); }
    rsSetupPointers(playerPool, idlePlayers);
  }

  void setMaxNumGroups(int newMax) 
  {
    reset();
    groupPlayerPool.resize(newMax);
    for(auto& p : groupPlayerPool)  {
      p.setDspResourcePool(&dspPool);
      p.setPlayStatusPointer(&playStatus); }
    rsSetupPointers(groupPlayerPool, idleGroupPlayers);
  }


  void preAllocateDspMemory() override
  {
    rsSamplerEngine2::preAllocateDspMemory() ;
  }


  int getNumUsedFilters() const { return dspPool.getNumUsedFilters(); }




};


} // namespace Sampler
} // namespace rosic

#endif