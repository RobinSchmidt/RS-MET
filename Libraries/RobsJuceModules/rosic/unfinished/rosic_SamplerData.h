#ifndef rosic_SamplerData_h
#define rosic_SamplerData_h

namespace rosic
{


/** Return codes for the setup functions. We use encodings as negative integers so we can use 
them also for functions which use positive integers as valid return values. */
enum rsReturnCode
{
  success        = -1,  //< Operation completed successfully. 
  nothingToDo    = -2,  //< There was nothing to actually do. State was already as desired.
  memAllocFail   = -3,  //< Memory allocation failure.
  invalidIndex   = -4,  //< An invalid index was passed.
  layerOverload  = -5,  //< Not enough free layers available (in e.g. new noteOn).
  notFound       = -6,  //< A region, group, sample or whatever was not found.
  fileLoadError  = -7,  //< A file could not be loaded (reasons: not found or failed alloc).
  notImplemented = -8   //< Feature not yet implemented (relevant during development).
};
// todo: make it an enum class, maybe include also return codes for inquiry functions such as for
// "unknown", etc. ...but maybe that's no good idea when we want to use it for functions which
// need to return valid integers (like, for numChannels, etc. - we could use negative numbers to
// encode such things)
// -maybe rename "success" to "completed" or "done" because "success" has actually a more general 
//   meaning: "nothingToDo" is also a kind of "success" (or maybe "workDone" or "workCompleted"
// -maybe include a general code for "failed" 
// -rename the layerOverload to a general "overload" or ressourcesUsedUp or something - to make the
//  enum more genrally useful




/** Data structure to define sample based instruments conforming to the sfz specification. */

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
  class PlaybackSetting  // rename to Setting...or maybe not - it's too generic
  {

  public:

    /** Enumeration of possible types of settings. These types correspond to the opcodes defined
    in the sfz specification. */
    enum Type
    {
      // Input controls:
      LoKey, HiKey, LoVel, HiVel,
      ControllerRangeLo, ControllerRangeHi, PitchWheelRange,  // 

      // Muted: convenient to switch regiosn or groups off wihthout removing them - check if 
      // sfz has such a thing

      // Pitch:
      PitchKeyCenter,

      // Amplitude:
      Volume, Pan, PanRule,
      AmpEnvAttack, AmpEnvDecay, AmpEnvSustain, AmpEnvRelease,
      // ToDo: Width, Position

      // Filtering:
      Delay,
      FilterCutoff, FilterResonance, FilterType,

      // Looping:


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

    /** Returns the default value for a given type of setting. */
    static float getDefaultValue(Type type);

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
    void addSetting(const PlaybackSetting& s) { settings.push_back(s); }
    // Maybe we should have a function setSetting that either adds a new setting or overwrites
    // an existing one, we may also need a function cleanUpSettings that keeps only the last 
    // setting of a particular kind in our list
    // todo: deprecate this in favor of setSetting - this should treat loKey, etc specially


    void clearSettings() { settings.clear(); }
    // should this also set loKey/hiKey and loVel/hiVel to 0/127?


    void setLoKey(uchar newKey) { loKey = newKey; }

    void setHiKey(uchar newKey) { hiKey = newKey; }

    void setLoVel(uchar newVel) { loVel = newVel; }

    void setHiVel(uchar newVel) { hiVel = newVel; }


    /** Sets a custom pointer that can be stored in the region object. @see getCustomPointer. */
    void setCustomPointer(const void* newPointer) { custom = newPointer; }


    void copyDataFrom(const OrganizationLevel* lvl);


    //---------------------------------------------------------------------------------------------
    // \name Inquiry

    /** @see setSample */
    const std::string& getSamplePath() const { return samplePath; }

    /** Returns a const reference to our playback settings. */
    const std::vector<PlaybackSetting>& getSettings() const { return settings; }

    /** Tries to find a setting of the given type in our settings array and returns the index of 
    the last stored setting of given type (and with given index, if applicable - like for 
    controllers), if a setting of the given type was found or -1 if such a setting wasn't found. 
    Typically, we want to have only one entry of a particular type in the array, but if for some 
    reason (like a poorly written sfz file) there is more than one, it's the last one that counts 
    (it will overwrite anything that came before). That's why we return the last index, i.e. we do
    a linear search starting at the end of the array. */
    int findSetting(PlaybackSetting::Type type, int index = -1);

    // todo: float getSetting(PlaybackSetting::Type, int index) this should loop through the 
    // settings to see, if it finds it and if not call the same method on the parent or return the
    // default value, if the parent is nullptr.

    /** Returns a pointer to the parent level which encloses this level. In a Region, this would 
    point to its enclosing Group, in a Group to its enclosing Instrument and in an Instrument, it
    would remain nullptr (unless we introduce an even higher level such as an "Ensemble"). */
    const OrganizationLevel* getParent() const { return parent; }

    /** Returns the lowest key at which this region will be played. */
    uchar getLoKey() const { return loKey; }

    /** Returns the higest key at which this region will be played. */
    uchar getHiKey() const { return hiKey; }

    /** Returns the lowest velocity at which this region will be played. */
    uchar getLoVel() const { return loVel; }

    /** Returns the highest velocity at which this region will be played. */
    uchar getHiVel() const { return hiVel; }

    /** Returns the generic pointer for custom satellite data or objects that are associated with
    this region. This pointer is intended to be used for some sort of audio stream object that is 
    used for accessing the sample data. It has been made a generic void pointer to decouple 
    rsSamplerData from the AudioFileStream class that is used in rsSamplerEngine. The 
    sampler-engine assigns this pointer with appropriate stream object and when retriveing them, 
    does an appropriate type cast. ToDo: try to find a better design, maybe move up into 
    baseclass */
    const void* getCustomPointer() const { return custom; }
    // replaces getSampleStream


    bool shouldPlayForKey(uchar key) const { return key >= loKey && key <= hiKey; }
    // todo: later maybe also take loKey/hiKey parent into account, if not nulltpr


  protected:

    /** Sets the audio stream object that should be used for this region. */
    //void setSampleStream(const AudioFileStream* newStream) { sampleStream = newStream; }

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
    // Yes, this will introduce coupling but it gets rid of the redundant storage and the unelegant
    // (and potentially unsafe) typecast of the void pointer. Maybe trying to decouple it amounts 
    // to the "speculative generality" antipattern here - we'll see...
    // https://refactoring.guru/smells/speculative-generality

    const void* custom = nullptr;

    OrganizationLevel* parent = nullptr;



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

  };


  //-----------------------------------------------------------------------------------------------
  /** A region is the lowest organizational level ins sfz. It contains a sample along with 
  performance settings including information for which keys and velocities the sample should be 
  played and optionally other constraints for when the the sample should be played and also 
  settings for pitch, volume, pan, filter, envelopes, etc. */
  class Region : public OrganizationLevel
  {

  public:

    void copyDataFrom(const Region* src);


    /** Return a pointer to the group to which this region belongs. */
    const Group* getGroup() const { return (const Group*) getParent(); }
    // maybe rename to getParentGroup or getEnclosingGroup


    bool operator==(const Region& rhs) const;


  private:

    friend class Group;  // do we need this? if not, get rid.
    friend class rsSamplerData;
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
    const Instrument* getInstrument() const { return (const Instrument*) getParent(); }

    /** Returns a pointer to the region with the given index within the group. */
    Region* getRegion(int i) const;
    // maybe return a const pointer

    bool operator==(const Group& rhs) const;
    // the comparison is quite strict in the sense that the settings must occur in the same order

  private:


    std::vector<Region*> regions;
    /**< Pointers to the regions belonging to this group. */

    friend class rsSamplerData;
    //friend class rsSamplerEngine;  // try to get rid
  };

  //-----------------------------------------------------------------------------------------------
  /** The instrument is the highest organizational level in sfz. There is actually no section 
  header in the .sfz file format for the whole instrument that corresponds to this class. The whole
  content of the sfz file *is* the instrument. But for consistency, we represent it by a class as 
  well. Maybe later an additional "Ensemble" or "Orchestra" level can be added on top. */
  class Instrument : public OrganizationLevel
  {

  public:

    int getNumGroups() const { return (int) groups.size(); }


    /** Returns a pointer to the region with the given index within the group. */
    Group* getGroup(int i) { return groups[i]; }


    const std::vector<PlaybackSetting>& getGroupSettings(int groupIndex) const
    { return groups[groupIndex]->getSettings(); }


    /** Returns true, if the given index i refers to a valid group within this instrument. */
    bool isGroupIndexValid(int i) const { return i >= 0 && i < (int)groups.size(); }


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

  /** Adds a region to the group with the given groupIndex and returns the region index within the
  group of the newly added region or -1 if the grouIndex was invalid. */
  int addRegion(int groupIndex, uchar loKey = 0, uchar hiKey = 127);

  /** Removes the region with given index from the group with given index and returns true, if this
  was successful and false if it wasn't, which happens when the index pair is invalid. */
  bool removeRegion(int groupIndex, int regionIndex);


  void setRegionCustomPointer(int gi, int ri, void* ptr)
  { instrument.groups[gi]->regions[ri]->setCustomPointer(ptr); }

  void setRegionSample(int gi, int ri, const std::string& samplePath)
  { instrument.groups[gi]->regions[ri]->setSamplePath(samplePath); }

  rsReturnCode setRegionSetting(int gi, int ri, PlaybackSetting::Type type, float value);

  rsReturnCode setGroupSetting(int gi, PlaybackSetting::Type type, float value);

  rsReturnCode setInstrumentSetting(PlaybackSetting::Type type, float value);


  /** Clears the whole instrument definition. */
  void clearInstrument() { instrument.clearGroups(); }
  //{ instrument.groups.clear(); }
  // todo: may wrap into instrument.clear() - don't access the groups array directly


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  bool isIndexPairValid(int groupIndex, int regionIndex) const
  { 
    int gi = groupIndex, ri = regionIndex;
    return instrument.isGroupIndexValid(gi) && instrument.groups[gi]->isRegionIndexValid(ri);
  }

  bool isGroupIndexValid(int i) const { return instrument.isGroupIndexValid(i); }


  int getNumGroups() const { return instrument.getNumGroups(); }

  int getNumRegions(int i) const { return instrument.groups[i]->getNumRegions(); }

  //const Group& getGroupRef(size_t i) const { return instrument.groups[i]; }

  //const Group& getGroupRef(int i) const { return *instrument.groups[i]; }
  // replace by getGroupPtr ..maybe rename to just getGroup

  const Group* getGroup(int i) const { return instrument.groups[i]; }

  const Region* getRegion(int gi, int ri) const 
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


}
#endif