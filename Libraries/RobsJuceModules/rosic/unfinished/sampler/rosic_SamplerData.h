#ifndef rosic_SamplerData_h
#define rosic_SamplerData_h

namespace rosic {
namespace Sampler {




//=================================================================================================

/** Data structure to define sample based instruments conforming to the sfz specification. From
here: https://sfzformat.com/legacy/ we quote:

"The basic component of an instrument is a region. An instrument then, is defined by one or more
regions. Multiple regions can be arranged in a group. Groups allow entering common parameters for
multiple regions." 

The parameters are represented as key/value pairs where the key is called an opcode in sfz 
terminology. So, regions can define playback parameters via opcodes, groups can define fallback 
values for these parameters for all regions belonging the group and the global section can define
fallback values for all its groups. We represent the 3 hierarchy levels of sfz (region, group, 
global) by the 3 internal classes Region, Group, Global which have a common basclass 
HierarchyLevel. The common baseclass stores the opcodes defined on the respective level. 
...tbc...   */

class SfzInstrument
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Default constructor. */
  SfzInstrument() {}

  ~SfzInstrument() { clearInstrument(); }

  /** Copy constructor.   */
  SfzInstrument(const SfzInstrument& d) { copy(d, *this); }

  /** Copy assignment operator. */
  SfzInstrument& operator=(const SfzInstrument& d) 
  { if(this != &d)  copy(d, *this); return *this; }

  //-----------------------------------------------------------------------------------------------
  // Forward declarations and abbreviations

  using uchar = unsigned char;
  class Region;
  class Group;
  class Global;


  //-----------------------------------------------------------------------------------------------
  // \name Helper classes

  //-----------------------------------------------------------------------------------------------
  /** Baseclass for the 3 organizational levels of the sfz specification, factoring out their
  commonalities. Subclasses are Region, Group, Global. 
  
  maybe drag out of the class, name it rsSamplerLevel or HierarchyLevel  */
  class HierarchyLevel  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // \name Setup

    /** Sets the sample to be used for this region. This should be a string that represents the
    path of the sample relative to some fixed root directory (typically the directory of the sfz
    file). */
    void setSamplePath(const std::string& newPath) { samplePath = newPath; }

    /** Sets up the given setting. This may either add it to our array of settings or overwrite the
    value, if there is already a value stored for the given type. If there are already more values
    stored for that type (a situation that may arise from badly written sfz files), it will
    overwrite the last one, because that's the one that actually counts. */
    void setSetting(const PlaybackSetting& s);
    // Maybe the function should take a second parameter that specifies the unit? Or maybe the 
    // PlaybackSetting should have a field for that? this will be needed to support a syntax
    // like fileg_depth=1000Hz or fileg_depth=70%

    /** Adds the given setting to our array of settings. */
    //void addSetting(const PlaybackSetting& s) { settings.push_back(s); }
    // Maybe we should have a function setSetting that either adds a new setting or overwrites
    // an existing one, we may also need a function cleanUpSettings that keeps only the last 
    // setting of a particular kind in our list
    // todo: deprecate this in favor of setSetting - this should treat loKey, etc specially

    /** Removes all the settings of given type and index. If the index is -1, it will be ignored, 
    so if it is an indexed opcode, the function can be used to remove all of them. Typically, 
    there should be at most one of any given type, but this may not be the case for poorly 
    written sfz files. Returns true, if there was at least one item actually removed, i.e. it 
    returns false only if no such setting was found in the array. */
    bool removeSetting(Opcode type, int index);

    /** Clears all the settings. */
    void clearSettings()
    {
      settings.clear();
      dspTypes.clear();
      modRoutings.clear();
    }
    // should this also set loKey/hiKey and loVel/hiVel back to their default values of 0/127? 
    // i actually think so...also the signalProcessors


    void setModulation(OpcodeType modSrcType, int modSrcIndex, 
      Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode);
    // rename to setModRouting or setModulationRouting

    void setModulation(const ModulationRouting& newRouting);

    bool removeModulation(OpcodeType modSrcType, int modSrcIndex, 
      Opcode modTarget, int modTargetIndex);

    /** Establishes a modualtion connection from the AmpEnvGen (if one exists) to an Amplifier unit
    at the end of the effect chain. If the last module in the effect chain already is an Amplifier
    with a nominal "amplitudeN" of zero, it will be used. Otherwise, another Amplifier will be 
    inserted at the end. */
    void connectAmpEnv();



    void setLoKey(uchar newKey) { loKey = newKey; }

    void setHiKey(uchar newKey) { hiKey = newKey; }

    void setLoVel(uchar newVel) { loVel = newVel; }

    void setHiVel(uchar newVel) { hiVel = newVel; }


    /** Sets a custom pointer that can be stored in the region object. @see getCustomPointer. */
    void setCustomPointer(const void* newPointer) { custom = newPointer; }


    /** Sets the parent level, i.e. the enclosing group for regions or the enclosing instrument for
    groups. */
    void setParent(HierarchyLevel* newParent) { parent = newParent; }



    void copyDataFrom(const HierarchyLevel* lvl);


    //---------------------------------------------------------------------------------------------
    // \name Inquiry

    /** @see setSample */
    const std::string& getSamplePath() const { return samplePath; }

    /** Returns a const reference to our playback settings. */
    const std::vector<PlaybackSetting>& getSettings() const { return settings; }

    /** Returns a const reference to our modulation routings. */
    const std::vector<ModulationRouting>& getModRoutings() const { return modRoutings; }

    const std::vector<ModulationRouting>& getModulationSettings() const { return modRoutings; }
    // replica! get rid! use only getModRoutings


    /** Returns the value of the given setting, if present. If not present, it will try to figure
    out the parent's setting and so on all the way up the (3-level) hierarchy. If such a setting is
    found in none of the levels, the default value for that setting will be returned. If accumulate
    is true, the settings of the different hierarchy levels will be added up, otherwise, the
    setting in the lower level will override the settimg in the enclosing higher level. */
    float getSettingValue(Opcode type, int index = -1, bool accumulate = false) const;

    /** Tries to find a setting of the given type in our settings array and returns the index of
    the last stored setting of given type (and with given index, if applicable - like for
    controllers), if a setting of the given type was found or -1 if such a setting wasn't found.
    Typically, we want to have only one entry of a particular type in the array, but if for some
    reason (like a poorly written sfz file) there is more than one, it's the last one that counts
    (it will overwrite anything that came before). That's why we return the last index, i.e. we do
    a linear search starting at the end of the array. */
    int findSetting(Opcode type, int index = -1) const;

    // todo: float getSetting(PlaybackSetting::Type, int index) this should loop through the 
    // settings to see, if it finds it and if not call the same method on the parent or return the
    // default value, if the parent is nullptr.

    /** Returns a pointer to the parent level which encloses this level. In a Region, this would
    point to its enclosing Group, in a Group to its enclosing Global and in an Global, it
    would remain nullptr (unless we introduce an even higher level such as an "Ensemble"). */
    const HierarchyLevel* getParent() const { return parent; }

    /** Returns the lowest key at which this region will be played. */
    uchar getLoKey() const { return loKey; }

    /** Returns the higest key at which this region will be played. */
    uchar getHiKey() const { return hiKey; }

    /** Returns the lowest velocity at which this region will be played. */
    uchar getLoVel() const { return loVel; }

    /** Returns the highest velocity at which this region will be played. */
    uchar getHiVel() const { return hiVel; }

    /** Returns the number of signal processors needed to play this region/group/etc. This includes
    all processors regardless of their kind, i.e. all the effects *and* modulators. */
    int getNumProcessors() const { return (int) dspTypes.size(); }

    /** Returns the number of effect processors needed to play this region/group/etc.  */
    int getNumEffects() const;
    // implement also getNumModulators. These count only the number of effects
    // and modulators

    /** Returns the number of processors of a given specific type. */
    int getNumProcessorsOfType(OpcodeType type) const;

    /** Returns true, iff the last effect in the effect chain is an amplifier. For an empty chain,
    it will return false. */
    bool isLastEffectAmplifier() const;


    /** Returns the generic pointer for custom satellite data or objects that are associated with
    this region. This pointer is intended to be used for some sort of audio stream object that is
    used for accessing the sample data. It has been made a generic void pointer to decouple
    SfzInstrument from the AudioFileStream class that is used in rsSamplerEngine. The
    sampler-engine assigns this pointer with appropriate stream object and when retrieving them,
    does an appropriate type cast. ToDo: try to find a better design, maybe move up into
    baseclass */
    const void* getCustomPointer() const { return custom; }
    // replaces getSampleStream..

    /** Returns true, iff this region should be played when the given key is pressed. */
    bool shouldPlayForKey(uchar key) const { return key >= loKey && key <= hiKey; }
    // ToDo: later maybe also take loKey/hiKey parent into account, if not nulltpr

    /** Returns a (const) reference to an array of processor types that is used by the engine to 
    build the dsp chain when this region should be played.*/
    const std::vector<OpcodeType>& getOpcodeTypeChain() const { return dspTypes; }



  protected:

    //---------------------------------------------------------------------------------------------
    // \name Misc

    /** For an opcode of given type, this function makes sure, that a corresponding signal 
    processor type is present in our signalProcessors array. It checks, if the right kind of 
    processor is already there and adds it, if not. You can also pass a parameter to require than
    more than one such processor must be present. */
    void ensureDspsPresent(Opcode opcodeType, int howMany = 1);

    /** Updates our dspTypes array from scratch from the settings array. */
    void updateDspsArray();

    //---------------------------------------------------------------------------------------------
    // \name Hardwired modulations. The depths of the hardwired modulation connections need some 
    // special case handling. The affected opcodes are: fileg_depth, ampeg_depth, pitcheg_depth, 
    // fillfo_depth, amplfo_depth, pitchlfo_depth

    /** Establishes modulation connections from the FilterEnv to the cutoff parameters of all Filter 
    modules. Called from setSetting as a special case handling. */
    void setFilterEnvDepth(float depthInCents);

    void setFilterLfoDepth(float depthInCents);


    /** Establishes a modulation connection from the AmpEnv to the final (last) Amplifier in the 
    effect chain subject to the constraint that this final Amplifier has an amplitudeN=0 setting. 
    If this is not the case or there aren't any Amplifiers present at all or the final effect is 
    something other than an Amplifier, an additional Amplifier with an amplitude=0 setting will be
    appended and the connection will be routed to this new Amplifier. */
    void setAmpEnvDepth(float depthInPercent);


    void setAmpLfoDepth(float depthInDecibels);



    //---------------------------------------------------------------------------------------------
    // \name Data


    /** Sets the audio stream object that should be used for this region. */
    //void setSampleStream(const AudioFileStream* newStream) { sampleStream = newStream; }

    //const AudioFileStream* sampleStream = nullptr;
    // try to get rid - that member should be added by rsSamplerEngine::Region which should be
    // a subclass of rsInstrumentDataSFZ::Region, and/or move up into baseclass. maybe to decouple
    // SfzInstrument from AudioFileStream, just keep it as pointer-to-void which client code may 
    // typecast to any sort of stream...or maybe that coupling makes sense?..hmm - not really.
    // maybe a pointer-to-void named customData should be stored in HierarchyLevel

    std::string samplePath;
    // This is the full (relative) path with respect to either the .sfz file or some predefined
    // sample directory. Actually, the string is redundant here when the "custom" pointer actually
    // points to an AudioFileStream, because that stream object also stores the path. Maybe revert
    // the custom void pointer to a pointer-to-AudioFileStream again. Yes, this will introduce 
    // coupling but it gets rid of the redundant storage and the unelegant (and potentially 
    //  unsafe) typecast of the void pointer. Maybe trying to decouple it amounts to the 
    // "speculative generality" antipattern here - we'll see...
    //   https://refactoring.guru/smells/speculative-generality
    // moreover, this part of the code is not the "hot" part anyway so optimizing it for size may
    // be pointless.

    const void* custom = nullptr;
    // Try to get rid of that - or at least, make it a typed pointer. It's supposed to hold the 
    // pointer to ths AudioStream. But we actually do not want to directly associate regions with
    // samples anymore anyway but instead delegate sample playback to some SamplePlayer subclass
    // of processor. If we do that, then *that* object needs to hold the pointer to the AudioStream


    HierarchyLevel* parent = nullptr;



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

    std::vector<PlaybackSetting> settings;
    /**< The settings which apply to this region/group/instrument, i.e. the opcodes along with 
    their values and, if appliable, index */

    std::vector<ModulationRouting> modRoutings;
    /**< Holds the modulation routings like egN_cutoffX, lfoN_amplitudeX, etc.. These are handled 
    separately from the other settings because they have a different format and do different 
    things. */

    std::vector<OpcodeType> dspTypes;
    /**< Listing of the types of signal processors used in this instrument in the same order like 
    how they should be applied (we assume a serial connection). But the list also contains the 
    modulation sources. The distinction between these two sorts of DSPs is not relevant here but
    will become relevant when actually assmebling the Player objects. */
    // maybe rename - it contains also the settings for the modulators like lfoN_freq



    // ModulationConnection is for the objects that are used in sample processing whereas 
    // ModulationRouting/Routing is used in regions/groups/etc



  };


  //-----------------------------------------------------------------------------------------------
  /** A region is the lowest hierarchy level ins sfz. It contains a sample along with
  performance settings including information for which keys and velocities the sample should be
  played and optionally other constraints for when the the sample should be played and also
  settings for pitch, volume, pan, filter, envelopes, etc. */

  class Region : public HierarchyLevel
  {

  public:

    void copyDataFrom(const Region* src);


    /** Return a pointer to the group to which this region belongs. */
    const Group* getGroup() const { return (const Group*)getParent(); }
    // maybe rename to getParentGroup or getEnclosingGroup


    bool operator==(const Region& rhs) const;


  private:

    //friend class Group;  // do we need this? if not, get rid.
    //friend class SfzInstrument;  // get rid!


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
  class Group : public HierarchyLevel
  {

  public:

    ~Group() { clearRegions(); }


    int addRegion(uchar loKey = 0, uchar hiKey = 127);      // todo: removeRegion, etc.

    int addRegion(Region* newRegion);

    bool removeRegion(int i);

    void copyDataFrom(const Group* scrGroup);

    void clearRegions();


    /** Returns the index of the given region within this group, if present or -1 if the region is
    not present in this group. */
    int getRegionIndex(const Region* region) const;

    /** Returns true, if the given index i refers to a valid region within this group. */
    bool isRegionIndexValid(int i) const { return i >= 0 && i < (int)regions.size(); }

    /** Returns the number of regions in this group. */
    int getNumRegions() const { return (int)regions.size(); }

    /** Return a pointer to the instrument to which this group belongs. */
    const Global* getInstrument() const { return (const Global*)getParent(); }
    // mybe rename to getGlobal ..or rename Global back to Instrument

    /** Returns a pointer to the region with the given index within the group. */
    Region* getRegion(int i) const;
    // maybe return a const pointer

    bool operator==(const Group& rhs) const;
    // the comparison is quite strict in the sense that the settings must occur in the same order

  private:






    std::vector<Region*> regions;
    /**< Pointers to the regions belonging to this group. */

    friend class SfzInstrument;
    //friend class rsSamplerEngine;  // try to get rid
  };

  //-----------------------------------------------------------------------------------------------
  /** The global instrument level is the highest hierarchy level in sfz. There is actually no 
  section header in the .sfz file format for the whole instrument that corresponds to this class. 
  The whole content of the sfz file *is* the instrument. But for consistency, we represent it by a 
  class as well. Maybe later an additional "Ensemble" or "Orchestra" level can be added on top. But 
  actually the .sfz already can apparently represent an ensemble. The spec says: "Each .sfz 
  definition file represents one or a collection of instruments.". I guess, what they have in mind
  is to use a dispatch based on the midi channel. Soo - i guess we are good with the 3 levels. */
  class Global : public HierarchyLevel
  {

  public:

    int getNumGroups() const { return (int)groups.size(); }


    /** Returns a pointer to the region with the given index i within the group or a nullptr if the
    index is invalid. */
    Group* getGroup(int i) 
    { 
      if(i < 0 && i >= (int)groups.size())
        return nullptr;
      return groups[i]; 
    }

    const std::vector<PlaybackSetting>& getGroupSettings(int groupIndex) const
    { return groups[groupIndex]->getSettings(); }


    /** Returns true, if the given index i refers to a valid group within this instrument. */
    bool isGroupIndexValid(int i) const { return i >= 0 && i < (int)groups.size(); }


    bool operator==(const Global& rhs) const;
    //{ return settings == rhs.settings && groups == rhs.groups; }


  //private:  // make protected later

    // maybe move to public
    int addGroup();      // todo: removeGroup, etc.

    int addGroup(Group* g);


    void clearGroups();


    std::vector<Group*> groups;
    // Should that be an array of pointers, too? Like the regions array in Group? That would make
    // the implementations of Group and Global more consistent but is actually technically not 
    // necessary. So, for the time being, let's keep it an array of direct value objects.

    friend class SfzInstrument;  // get rid
  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  // todo: factor out some code from rsSamplerEngine - the names of the correponding setters here 
  // and there should match and rsSamplerEngine should call functions from here and perhaps do
  // additional stuff, if necessary

  int addGroup() { return global.addGroup(); }

  int addGroup(Group* newGroup) { return global.addGroup(newGroup); }

  /** Adds a region to the group with the given groupIndex and returns the region index within the
  group of the newly added region or -1 if the grouIndex was invalid. */
  int addRegion(int groupIndex, uchar loKey = 0, uchar hiKey = 127);

  /** Removes the region with given index from the group with given index and returns true, if this
  was successful and false if it wasn't, which happens when the index pair is invalid. */
  bool removeRegion(int groupIndex, int regionIndex);


  void setRegionCustomPointer(int gi, int ri, void* ptr)
  {
    global.groups[gi]->regions[ri]->setCustomPointer(ptr);
  }

  void setRegionSample(int gi, int ri, const std::string& samplePath)
  {
    global.groups[gi]->regions[ri]->setSamplePath(samplePath);
  }

  rsReturnCode setRegionSetting(int gi, int ri, Opcode type, float value, int index = -1);
  rsReturnCode setGroupSetting(int gi, Opcode type, float value, int index = -1);
  rsReturnCode setInstrumentSetting(Opcode type, float value, int index = -1);  // rename to setGlobalSetting for sfz compliance

  // for setting up a modulation connection
  rsReturnCode setRegionModulation(int gi, int ri, OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode);

  rsReturnCode setGroupModulation(int gi, OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode);

  rsReturnCode setInstrumentModulation(OpcodeType modSrcType, int modSrcIndex,
    Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode);
  // It actually always returns rsReturnCode::success. I'm not sure, if we should keep the 
  // return-value or declare the function void. The former is nicer from a consistency perspective 
  // but the value has actually no information. On the other hand, the returned success report 
  // isn't false either - it's just that the function is always supposed to succeed...hmmm...


  rsReturnCode removeRegionSetting(int gi, int ri, Opcode type, int index);
  rsReturnCode removeGroupSetting(int gi, Opcode type, int index);
  rsReturnCode removeInstrumentSetting(Opcode type, int index);

  rsReturnCode removeRegionModulation(int gi, int ri, OpcodeType modSrcType, 
    int modSrcIndex,  Opcode modTarget, int modTargetIndex);


  rsReturnCode clearRegionSettings(int gi, int ri);

  void clearAllRegionSettings();
  void clearAllGroupSettings();
  void clearAllInstrumentSettings();
  void clearAllSettings();



  /** Clears the whole instrument definition. */
  void clearInstrument()
  {
    global.clearGroups();
    global.clearSettings();
    //signalProcessors.clear();
  }
  //{ instrument.groups.clear(); }
  // todo: may wrap into instrument.clear() - don't access the groups array directly


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  bool isIndexPairValid(int groupIndex, int regionIndex) const
  {
    int gi = groupIndex, ri = regionIndex;
    return global.isGroupIndexValid(gi) && global.groups[gi]->isRegionIndexValid(ri);
  }

  bool isGroupIndexValid(int i) const { return global.isGroupIndexValid(i); }


  int getNumGroups() const { return global.getNumGroups(); }

  int getNumRegions(int i) const { return global.groups[i]->getNumRegions(); }

  //const Group& getGroupRef(size_t i) const { return instrument.groups[i]; }

  //const Group& getGroupRef(int i) const { return *instrument.groups[i]; }
  // replace by getGroupPtr ..maybe rename to just getGroup

  const Global* getGlobal() const { return &global; }

  const Group* getGroup(int i) const { return global.groups[i]; }

  const Region* getRegion(int gi, int ri) const { return global.groups[gi]->regions[ri]; }

  Region* getRegionMutable(int gi, int ri) const { return global.groups[gi]->regions[ri]; }

/** Returns a const reference to the playback settings if the i-th group. */
//const std::vector<PlaybackSetting>& getGroupSettings(size_t i) const 
//{ return instrument.getGroupSettings(i); }

  bool operator==(const SfzInstrument& rhs) const { return global == rhs.global; }




  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Produces the string that represents the settings in an sfz-file compliant format, i.e. a
  string that can be written into an .sfz file. */
  std::string getAsSFZ() const;


  /** Sets up this data object according to the given string which is supposed to represent the
  contents of an .sfz file. */
  rsReturnCode setFromSFZ(const std::string& sfzFileContents);
  // todo: return a return-code, including unknownOpcode, invalidValue, invalidIndex, ...




  /** Writes the data represented by this object into an .sfz file with given path. Warning: This
  function has no safeguards against overwriting an existing file - it will just do it! */
  bool saveToSFZ(const char* path) const;
  // todo: return a return-code, including fileWriteError

  /** Sets up this object according to a given .sfz file.   */
  rsReturnCode loadFromSFZ(const char* path);
  // todo: document return values - should probably be one of: success, loadError, parseError


//protected:  // preliminarily commented - make protected again later

  Global global;

protected:

  static void writeSettingToString(const PlaybackSetting& setting, std::string& str);


  static void writeModRoutingToString(const ModulationRouting& routing, std::string& str);


  static PlaybackSetting getSettingFromString(
    const std::string& opcode, const std::string& value);

  static ModulationRouting getModRoutingFromString(
    const std::string& opStr, const std::string& valStr);


  static void copy(const SfzInstrument& src, SfzInstrument& dst);

};

}
}

#endif