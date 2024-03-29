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
global) by the 3 internal classes Region, Group, Global which havea common basclass 
HierarchyLevel. The common baseclass stores the opcodes defined on the respective level. 
...tbc...   */

class SfzInstrument // todo: move into its own pair of .h/.cpp files, rename to SfzInstrument
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
  // \name Helper classes (maybe drag out - we'll see)

  //-----------------------------------------------------------------------------------------------
  /** A class to represent various playback settings of a region, group or instrument. Such
  settings include constraints for the circumstances under which a particular sample should be
  played. Key- and velocity ranges are the obvious primary constraints (and they therefore are
  directly baked into the Region class below), but sfz defines many more. But settings don't
  need to be playback constraints. Other types are things like envelope settings, filter
  frequencies, etc. */
  class PlaybackSetting  // rename to Setting...or maybe not - it's too generic
  {

  public:


    PlaybackSetting(Opcode type, float value, int index = -1) 
    { 
      this->type  = type; 
      this->value = value; 
      this->index = index;
    }
    // maybe remove the -1 default argument. callers should be explicit...maybe...or: check with an
    // assert that no index applies to the givne Opcode when index == -1

    Opcode getOpcode() const { return type; }
    // rename to getOpcode


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

    /** Returns the default value for a given type of setting. */
    static float getDefaultValue(Opcode op, int index);

    /** For a given type of opcode, this function returns the type of signal processor to which the
    opcode applies, e.g. SignalProcessorType::Filter for Type::FilterCutoff. */
    static DspType getTargetDspType(Opcode op);
    // rename to getTargetDspType


    DspType getTargetDspType() const { return getTargetDspType(type); }


    bool operator==(const PlaybackSetting& rhs) const
    {
      return type == rhs.type && value == rhs.value && index == rhs.index;
    }


  private:

    Opcode type  = Opcode::Unknown;  // rename type to opcode
    float  value = 0.f;
    int    index = -1;  //< Used e.g. for conrol-change settings. Is -1, if not applicable.
    // Maybe use 1 as default - if there's only one such setting anyway, that seems appropriate
    // index should always be a positive real number. But maybe that's not such a good idea - see
    // comment in SfzInstrument::writeSettingToString in the cpp file. For certain things, we need
    // a code for "not applicable".
  };

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
    }
    // should this also set loKey/hiKey and loVel/hiVel back to their default values of 0/127? 
    // i actually think so...also the signalProcessors


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

    /** Returns the number of signal processors needed to play this region/group/etc. */
    int getNumProcessors() const { return (int) dspTypes.size(); }



    /** Returns the generic pointer for custom satellite data or objects that are associated with
    this region. This pointer is intended to be used for some sort of audio stream object that is
    used for accessing the sample data. It has been made a generic void pointer to decouple
    SfzInstrument from the AudioFileStream class that is used in rsSamplerEngine. The
    sampler-engine assigns this pointer with appropriate stream object and when retrieving them,
    does an appropriate type cast. ToDo: try to find a better design, maybe move up into
    baseclass */
    const void* getCustomPointer() const { return custom; }
    // replaces getSampleStream

    /** Returns true, iff this region should be played when the given key is pressed. */
    bool shouldPlayForKey(uchar key) const { return key >= loKey && key <= hiKey; }
    // ToDo: later maybe also take loKey/hiKey parent into account, if not nulltpr

    /** Returns a (const) reference to an array of processor types that is used by the engine to 
    build the dsp chain when this region should be played.*/
    const std::vector<DspType>& getDspTypeChain() const { return dspTypes; }
    // rename to getDspTypeChain


  protected:

    /** For an opcode of given type, this function makes sure, that a corresponding signal 
    processor type is present in our signalProcessors array. It checks, if the right kind of 
    processor is already there and adds it, if not. You can also pass a parameter to require than
    more than one such processor must be present. */
    void ensureDspsPresent(Opcode opcodeType, int howMany = 1);

    /** Updates our dspTypes array from scratch from the settings array. */
    void updateDspsArray();




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
    /**< The settings which applay to this region/group/instrument, i.e. the opcodes along with 
    their values and, if appliable, index */

    std::vector<DspType> dspTypes;
    /**< Listing of the types of signal processors used in this instrument in the same order like 
    how they should be applied (we assume a serial connection). */

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

  rsReturnCode setInstrumentSetting(Opcode type, float value, int index = -1);
  // rename to setGlobalSetting for sfz compliance

  rsReturnCode removeRegionSetting(int gi, int ri, Opcode type, int index);

  rsReturnCode clearRegionSettings(int gi, int ri);

  rsReturnCode removeGroupSetting(int gi, Opcode type, int index);

  rsReturnCode removeInstrumentSetting(Opcode type, int index);


  void clearAllRegionSettings();
  void clearAllGroupSettings();
  void clearAllInstrumentSettings();
  void clearAllSettings();



  /** Clears the whole instrument definition. */
  void clearInstrument()
  {
    global.clearGroups();
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

  static PlaybackSetting getSettingFromString(
    const std::string& opcode, const std::string& value);

  static void copy(const SfzInstrument& src, SfzInstrument& dst);

};

}
}

#endif