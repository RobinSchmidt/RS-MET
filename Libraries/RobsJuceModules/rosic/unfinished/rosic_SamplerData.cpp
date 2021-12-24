namespace rosic {
namespace Sampler {

//-------------------------------------------------------------------------------------------------
// The internal classes

float rsSamplerData::PlaybackSetting::getDefaultValue(Opcode type)
{
  // New:
  SfzOpcodeTranslator* t = SfzOpcodeTranslator::getInstance();
  return t->opcodeDefaultValue(type);


  /*
  using TP = Opcode;
  switch(type)
  {
  case TP::LoKey:          return 0.f;
  case TP::HiKey:          return 127.f;
  case TP::LoVel:          return 0.f;
  case TP::HiVel:          return 127.f;

  case TP::Volume:         return 0.f;
  case TP::Pan:            return 0.f;

  case TP::PitchKeyCenter: return 60.f;
  case TP::Transpose:      return  0.f;
  case TP::Tune:           return  0.f;

  case TP::Delay:          return 0.f;
  case TP::Offset:         return 0.f;


  // Filter: the spec says, if the cutoff is not specified, it should deactivate the filter...how
  // can we capture this here? maybe return -1 as code?

  // Extensions:
  case TP::DistShape:      return 0.f;  // maybe rename to DistortShape, ...
  case TP::DistDrive:      return 0.f;
  //case TP::DistGain:       return 0.f;

  }

  RAPT::rsError("Unknown type of PlaybackSetting, i.e. unknown sfz opcode.");
  return 0.f;  // maybe we should return NaN? but no - that would be evil!
  */
}

DspType rsSamplerData::PlaybackSetting::getTargetProcessorType(Opcode type)
{
  // New: 
  SfzOpcodeTranslator* t = SfzOpcodeTranslator::getInstance();
  return t->opcodeToProcessor(type);
  // try to get rid of this method entirely - let the caller do the stuff directly
  // ...does not yet work

  /*
  using TP = Opcode;
  using SP = DspType;
  switch(type)
  {
  case TP::LoKey:          return SP::SamplePlayer;
  case TP::HiKey:          return SP::SamplePlayer;
  case TP::LoVel:          return SP::SamplePlayer;
  case TP::HiVel:          return SP::SamplePlayer;
  case TP::Volume:         return SP::SamplePlayer;
  case TP::Pan:            return SP::SamplePlayer;
  case TP::PitchKeyCenter: return SP::SamplePlayer;
  case TP::Transpose:      return SP::SamplePlayer;
  case TP::Tune:           return SP::SamplePlayer;
  case TP::Delay:          return SP::SamplePlayer;
  case TP::Offset:         return SP::SamplePlayer;

  case TP::FilType:        return SP::Filter;
  case TP::Cutoff:         return SP::Filter;
  case TP::Resonance:      return SP::Filter;    // shorten to FilterReso or FiltReso

  case TP::DistShape:      return SP::WaveShaper;
  case TP::DistDrive:      return SP::WaveShaper;
  }

  // Maybe the the amount code can be reduced when we ensure that the TP types are sorted by their
  // SP in the enum, i.e. we could write:
  //   if(type < TP::PlayerOpcodesEnd) return SP::SamplePlayer;
  //   if(type < TP::FilterOpcodesEnd) return SP::Filter;
  // etc. But this requires discipline when adding new opcodes to respect the ordering.

  RAPT::rsError("Unknown type of PlaybackSetting, i.e. unknown sfz opcode.");
  return SP::Unknown;
  */
}
// This function should be dragged out of the class. It would be most convenient to have some 
// sort global or singleton object from which such information about the various relationships
// between opcode enum indices, their sfz strings, their associated signal processors, etc. can 
// be retrieved. Something like SfzOpcodeInfo with members like getProcessorForOpcode, 
// getStringForOpcode, getOpcodeForString, getOpcodesForProcessor. etc. It would act as a kind 
// of globally available database and could implement efficient mappings. Mapping opcode indices
// to anything should always be doable in O(1), Mapping opcode sfz-strings to indices may be 
// doable O(log(N)) where n is the number of opcodes, if we also store an alphabetically sorted
// list, etc. ...but these are implementation details. The class should document the complexity
// of each such mapping operation. Mapping the strings to indces could perhaps even be done by 
// a std::map to get O(1) on average - we'll see....

//-------------------------------------------------------------------------------------------------

void rsSamplerData::OrganizationLevel::ensureProcessorPresent(Opcode opcodeType)
{ 
  using namespace RAPT;
  using SPT = DspType;
  SPT dspType = PlaybackSetting::getTargetProcessorType(opcodeType);
  rsAssert( dspType != SPT::Unknown );
  if( dspType == SPT::SamplePlayer || dspType == SPT::Unknown )
    return;
    // The sample-player at the start of the processing chain doesn't really count as bona-fide DSP
    // processor. It's always there, there's always exactly one and it behaves quite differently 
    // from the rest. We need it among the types for consistency, though.
  if( !rsContains(signalProcessors, dspType) )
    signalProcessors.push_back(dspType);
}

void rsSamplerData::OrganizationLevel::setSetting(const PlaybackSetting& s)
{
  using TP = Opcode;
  TP t = s.getType();

  // Handle the lo/hi key/vel opcodes as special cases:
  if(t == TP::LoKey) { loKey = (uchar)s.getValue(); return; }
  if(t == TP::HiKey) { hiKey = (uchar)s.getValue(); return; }
  if(t == TP::LoVel) { loVel = (uchar)s.getValue(); return; }
  if(t == TP::HiVel) { hiVel = (uchar)s.getValue(); return; }
  // ToDo: maybe we should assert that the value is an integer in the range 0..127

  // All other settings are handled by either overwriting the last setting of that type in our 
  // array, if present or by appending the setting, if not present:
  int i = findSetting(t);
  if(i != -1)
    settings[i] = s;
  else
  {
    settings.push_back(s);
    ensureProcessorPresent(t);
    // The order in which the processors appear in the chain should reflect the order in which 
    // their opcodes appear in the sfz (or, if setup is done programmatically, the order in which
    // the opcodes were added). The first opcode applying to a particular kind of processor 
    // counts. For example, if the opcodes are added in the order FilterCutoff, DistDrive,
    // FilterResonance, the filter appears before the waveshaper in the DSP chain...maybe with 
    // some exceptions for opcodes that apply to processors that must be at fixed positions in the 
    // chain such as the SamplePlayer. ...hmm...but what, if we want two processors of the same 
    // kind? like filter1 -> waveshaper -> filter2. Maybe we could re-use the "index" member
    // PlaybackSetting (which was originally supposed to indicate a midi controller for 
    // control-change). But then, the sfz-parser would not only have to look for "cutoff" but also
    // cutoff1, cutoff2, etc. and translate that into an appropriate PlaybackSetting, i.e. one
    // with the index variable set. But that seems doable.
  }
}

bool rsSamplerData::OrganizationLevel::removeSetting(Opcode type)
{
  bool wasRemoved = false;
  for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
    if(settings[i].getType() == type) {
      RAPT::rsRemove(settings, i);
      wasRemoved = true;
    }
  }
  return wasRemoved;
  // We can't use size_t for i because the -1 would create an access violation when size() = 0
}

void rsSamplerData::OrganizationLevel::copyDataFrom(const OrganizationLevel* lvl)
{
  samplePath = lvl->samplePath;
  settings   = lvl->settings;

  // not sure, if the pointers should be copied - maybe not:
  //custom = lvl->custom;
  //parent = lvl->parent;
}

float rsSamplerData::OrganizationLevel::getSettingValue(
  Opcode type, int index, bool accumulate) const
{
  float val = PlaybackSetting::getDefaultValue(type);  // init to global fallback value
  int i = findSetting(type, index);

  if(i != -1)                        // setting was found 
    val = settings[i].getValue();    //   -> retrieve value
  else                               // setting was not found
    if(parent != nullptr)            //   -> try to fall back to parent's value   
      val = parent->getSettingValue(type, index, accumulate);

  if(accumulate && parent != nullptr)
    return val + parent->getSettingValue(type, index, accumulate);
  else
    return val;
}
// needs tests...maybe the accumulation is not always a straightforward addition. Think of
// Pan for example, which has a more complex accumulation behavior. 3 panners in series with 
// setting of +50 will not have the same effect as one with +150 (which is out of the legal range 
// anyway)...but maybe the accumulate feature is not even needed in this function - may get rid of
// it

int rsSamplerData::OrganizationLevel::findSetting(Opcode type, int index) const
{
  for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
    if(settings[i].getType() == type && settings[i].getIndex() == index)
      return (int)i;
  }
  return -1;
}

void rsSamplerData::Region::copyDataFrom(const Region* src)
{
  rsSamplerData::OrganizationLevel::copyDataFrom(src);
  loKey = src->loKey;
  hiKey = src->hiKey;
  loVel = src->loVel;
  hiVel = src->hiVel;
  int dummy = 0;
}

bool rsSamplerData::Region::operator==(const rsSamplerData::Region& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= loKey == rhs.loKey;
  equal &= hiKey == rhs.hiKey;
  equal &= loVel == rhs.loVel;
  equal &= hiVel == rhs.hiVel;
  equal &= samplePath == rhs.samplePath;
  return equal;
  // What about the customPointer? should we require that to be equal, too?
}

int rsSamplerData::Group::addRegion(uchar loKey, uchar hiKey)
{
  rsSamplerData::Region* r = new rsSamplerData::Region;
  r->setLoKey(loKey);
  r->setHiKey(hiKey);
  return addRegion(r);
}

int rsSamplerData::Group::addRegion(Region* r)
{
  r->setParent(this);
  regions.push_back(r);
  return ((int)regions.size()) - 1;
}

bool rsSamplerData::Group::removeRegion(int i)
{
  if(i < 0 || i >= (int)regions.size())
    return false;
  delete regions[i];
  RAPT::rsRemove(regions, i);
  return true;
}

void rsSamplerData::Group::copyDataFrom(const Group* src)
{
  rsSamplerData::OrganizationLevel::copyDataFrom(src);
  clearRegions();
  settings = src->getSettings();
  for(int i = 0; i < src->getNumRegions(); i++) {
    const rsSamplerData::Region* srcRegion = src->getRegion(i);
    rsSamplerData::Region* dstRegion = new rsSamplerData::Region;
    dstRegion->copyDataFrom(srcRegion);
    addRegion(dstRegion);
  }
}

void rsSamplerData::Group::clearRegions()
{
  for(size_t i = 0; i < regions.size(); i++)
    delete regions[i];
  regions.clear();
}

int rsSamplerData::Group::getRegionIndex(const rsSamplerData::Region* region) const
{
  for(size_t i = 0; i < regions.size(); i++)
    if(regions[i] == region)
      return (int)i;
  return -1;
}

rsSamplerData::Region* rsSamplerData::Group::getRegion(int i) const
{
  if(i < 0 || i >= (int)regions.size()) {
    RAPT::rsError("Invalid region index");
    return nullptr;
  }
  return regions[i];
}

bool rsSamplerData::Group::operator==(const rsSamplerData::Group& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= regions.size() == rhs.regions.size();
  if(!equal) return false;
  for(size_t i = 0; i < regions.size(); i++)
    equal &= *(regions[i]) == *(rhs.regions[i]);
  return equal;
}

int rsSamplerData::Instrument::addGroup()
{
  rsSamplerData::Group* g = new rsSamplerData::Group;
  return addGroup(g);
}

int rsSamplerData::Instrument::addGroup(rsSamplerData::Group* g)
{
  g->parent = this;
  groups.push_back(g);
  return ((int)groups.size()) - 1;
}

void rsSamplerData::Instrument::clearGroups()
{
  for(size_t i = 0; i < groups.size(); i++)
    delete groups[i];
  groups.clear();
}

bool rsSamplerData::Instrument::operator==(const rsSamplerData::Instrument& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= groups.size() == rhs.groups.size();
  if(!equal) return false;
  for(size_t i = 0; i < groups.size(); i++)
    equal &= *(groups[i]) == *(rhs.groups[i]);
  return equal;
}

//-------------------------------------------------------------------------------------------------
// The actual rsSamplerData class:

int rsSamplerData::addRegion(int gi, uchar loKey, uchar hiKey)
{
  if(gi < 0 || gi >= (int)instrument.groups.size()) {
    RAPT::rsError("Invalid group index");
    return -1;
  }
  int ri = instrument.groups[gi]->addRegion(loKey, hiKey);  // region index within its group
  return ri;
}

bool rsSamplerData::removeRegion(int gi, int ri)
{
  if(gi < 0 || gi >= (int)instrument.groups.size()) {
    RAPT::rsError("Invalid group index");
    return false;
  }
  return instrument.groups[gi]->removeRegion(ri);
}

rsReturnCode rsSamplerData::setRegionSetting(int gi, int ri, Opcode type, float value, int index)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex;
  }
  instrument.groups[gi]->regions[ri]->setSetting(PlaybackSetting(type, value, index));  // new

  //instrument.groups[gi]->regions[ri]->settings.push_back(PlaybackSetting(type, value)); // old
  // Preliminary. We need to figure out, if that setting already exists and if so, just change 
  // its value instead of pushing another value for the same parameter. Implement it in way so we 
  // can call it here as: settings->set(PlaybackSetting::Type type, float value)

  return rsReturnCode::success;
}

rsReturnCode rsSamplerData::removeRegionSetting(int gi, int ri, Opcode type)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex;
  }
  bool wasRemoved = instrument.groups[gi]->regions[ri]->removeSetting(type);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode rsSamplerData::setGroupSetting(int gi, Opcode type, float value, int index)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex;
  }

  instrument.groups[gi]->settings.push_back(PlaybackSetting(type, value, index));
  // Preliminary. We need to figure out, if that setting already exists and if so, just change 
  // its value instead of pushing another value for the same parameter

  return rsReturnCode::success;
}

rsReturnCode rsSamplerData::removeGroupSetting(int gi, Opcode type)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex;
  }
  bool wasRemoved = instrument.groups[gi]->removeSetting(type);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode rsSamplerData::setInstrumentSetting(Opcode type, float value, int index)
{
  instrument.settings.push_back(PlaybackSetting(type, value, index));
  // Preliminary. see above

  return rsReturnCode::success;
}

rsReturnCode rsSamplerData::removeInstrumentSetting(Opcode type)
{
  bool wasRemoved = instrument.removeSetting(type);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

void rsSamplerData::clearAllRegionSettings()
{
  for(size_t gi = 0; gi < instrument.groups.size(); gi++)
    for(size_t ri = 0; ri < instrument.groups[gi]->regions.size(); ri++)
      instrument.groups[gi]->regions[ri]->clearSettings();
}

void rsSamplerData::clearAllGroupSettings()
{
  for(size_t gi = 0; gi < instrument.groups.size(); gi++)
    instrument.groups[gi]->clearSettings();
}

void rsSamplerData::clearAllInstrumentSettings()
{
  instrument.clearSettings();
  //instrument.settings.clear();
}

void rsSamplerData::clearAllSettings()
{
  clearAllRegionSettings();
  clearAllGroupSettings();
  clearAllInstrumentSettings();
}
// needs test


std::string rsSamplerData::getAsSFZ() const
{
  auto writeSettingsToString = [](const OrganizationLevel* lvl, std::string& str)
  {
    auto toStr = [](const uchar c) { return std::to_string(c); }; // uchar to string
    const std::string& samplePath = lvl->getSamplePath();
    if(!samplePath.empty()) str += "sample=" + samplePath + '\n';
    if(lvl->getLoKey() !=   0) str += "lokey=" + toStr(lvl->getLoKey()) + '\n';
    if(lvl->getHiKey() != 127) str += "hikey=" + toStr(lvl->getHiKey()) + '\n';
    if(lvl->getLoVel() !=   0) str += "lovel=" + toStr(lvl->getLoVel()) + '\n';
    if(lvl->getHiVel() != 127) str += "hivel=" + toStr(lvl->getHiVel()) + '\n';
    using SettingsRef = const std::vector<PlaybackSetting> &;
    SettingsRef settings = lvl->getSettings();
    for(size_t i = 0; i < settings.size(); i++)
      writeSettingToString(settings[i], str);
  };

  std::string str;
  writeSettingsToString(&instrument, str);
  for(int gi = 0; gi < getNumGroups(); gi++) {
    str += "<group>\n";
    writeSettingsToString(getGroup(gi), str);
    for(int ri = 0; ri < getNumRegions(gi); ri++) {
      str += "<region>\n";
      writeSettingsToString(getRegion(gi, ri), str);
    }
  }
  return str;

  // ToDo: write lokey/hikey settings into the string, they are stored directly in the Region 
  // object and not also in the settings. Maybe they should be. That would simplify the 
  // serialization but that may complicate other things due to the introduced redundancy and 
  // therefore extra care to keep the data consistent. That raises the question, if groups and
  // instruments also can define lokey/hikey settings and how they are interpreted. If so, maybe
  // these lokey/hikey members should be moved into the OrganizationLevel baseclass. 
  // ...done ...verify and delete comment

  // ToDo: figure out how SFZPlayer behaves with respect to this maybe by defining those opcodes
  // at all 3 levels - i guess, it will use the most restrictive setting of all of them
}

void rsSamplerData::setFromSFZ(const std::string& str)
{
  clearInstrument();
  if(str.empty())
    return;
  size_t endOfFile = std::numeric_limits<size_t>::max();

  // Extracts the subtring starting at startIndex up to (and excluding) the next newline '\n' 
  // charcater. If there is no '\n', it will return the string from startIndex up to its end:
  std::string sep(" \n");  // allowed seperator characters
  auto getToken = [&](const std::string& str, size_t startIndex)
  {
    int start  = (int)startIndex;
    int length = -1;  // initial value should not matter
    rosic::rsFindToken(str, sep, &start, &length);
    return str.substr(start, length);
  };

  // Sets up one setting in lvl given in the format "opcode=value":
  auto setupSetting = [&](OrganizationLevel* lvl, const std::string& str)
  {
    size_t splitIndex = str.find('=', 0);
    std::string opcode = str.substr(0, splitIndex);
    std::string value  = str.substr(splitIndex+1, str.length() - splitIndex - 1);

    if(opcode == "sample") {     // needs to be treated in a special way
      lvl->setSamplePath(value);
      return;
    }

    PlaybackSetting ps = getSettingFromString(opcode, value);
    lvl->setSetting(ps);
    //lvl->addSetting(ps); 
    // The difference between setSetting and addSetting is that setSetting first scans the array
    // of settings to figure out, if such a setting is already present and if so overwrites it
    // whereas addSetting just appends the setting no matter what. Clearly the latter is faster but
    // the former is more fool-proof aginst badly written sfz files. I'm not yet sure, which 
    // version should be used...we'll see...
  };

  // Sets up the given level according to the given string which is supposed to contain one setting
  // per line in the format "opocde=value\n":
  auto setupLevel = [&](OrganizationLevel* lvl, const std::string& str)
  {
    size_t start = 0;
    while(true)
    {
      std::string token = getToken(str, start); // extract one token at at time
      if(token.length() == 0)
        break;
      setupSetting(lvl, token);                 // set a setting from this token
      start += token.length() + 1;
      //if(start >= str.length()) break;          // may be superfluous?
    }
  };


  //std::string group  = "<group>\n";   // not sure, whether we should include the \n
  //std::string region = "<region>\n";  // ..probably not
  std::string group  = "<group>";
  std::string region = "<region>";
  size_t Lg = group.length();
  size_t Lr = region.length();

  std::string tmp;                    // for extracted substrings (maybe use string_view)

  // Find start and end index in the string for the first group:
  size_t i0 = str.find(group, 0);
  size_t i1 = str.find(group, i0+1);

  // Set up instrument level settings:
  tmp = str.substr(0, i0);
  setupLevel(&instrument, tmp);

  // Loop over the the groups within the instrument definition:
  bool allGroupsDone = false;
  while(!allGroupsDone)
  {
    if(i1 == endOfFile) {
      allGroupsDone = true;
      i1 = str.length() - 1;
    }

// Extract substring with group definition and add a new group to the instrument:
    std::string groupDef = str.substr(i0, i1-i0); // group definition (ToDo: use string_view)
    int gi = instrument.addGroup();
    Group* g = instrument.getGroup(gi);
    g->parent = &instrument;

    // Find start and end index in the string for the first region within the current group:
    size_t j0 = str.find(region, i0);
    size_t j1 = str.find(region, i0+1);

    // Set up group level settings:
    tmp = str.substr(i0+Lg, j0-i0-Lg);
    setupLevel(g, tmp);

    // Loop over the the regions within the group definition:
    bool allRegionsDone = false;
    while(!allRegionsDone)
    {
      // Find start and end index of next region definition:
      j0 = groupDef.find(region, j1);
      RAPT::rsAssert(j0 != endOfFile);  // for debug - gets triggered when we have empty regions
      j1 = groupDef.find(region, j0+1);
      if(j1 == endOfFile) {
        allRegionsDone = true;
        j1 = groupDef.length() - 1;
      }

// Extract substring with region definition and add a new region to the group:
      std::string regionDef = groupDef.substr(j0, j1-j0); // region definition (ToDo: use string_view)
      int ri = g->addRegion();
      Region* r = g->getRegion(ri);
      r->setParent(g);

      // Set up region level settings:
      tmp = groupDef.substr(j0+Lr, j1-j0-Lr);
      setupLevel(r, tmp);
      int dummy = 0;
    }

    // Find start and end index of next group defintion:
    i0 = str.find(group, i1);      // start index of the group in the string
    i1 = str.find(group, i0+1);    // end index of the group in the string
    int dummy = 0;
  }

  // ToDo: 
  // -The general structure of the nested region is the similar to the enclosing group block 
  //  -> try to refactor to get rid of the duplication (maybe it can be implemented recursively)
  // -Maybe use string_view for the extracted substrings to avoid copying the data:
  //  https://en.cppreference.com/w/cpp/header/string_view
}

bool rsSamplerData::saveToSFZ(const char* path) const
{
  std::string sfz = getAsSFZ();
  return rsWriteStringToFile(path, sfz.c_str());
}
// this has no safeguards against overwriting an existing file!

bool rsSamplerData::loadFromSFZ(const char* path)
{
  // just for debug, to figure out, in which directory the mac expects the sfz file:
  //rsWriteStringToFile("TestFile.sfz", "blablabla");
  // that fails, too with an "Unable to open file" error. Could it have to do with permission?

  char* c_str = rsReadStringFromFile(path);
  if(c_str)
  {
    std::string sfz(c_str);
    setFromSFZ(sfz);
    free(c_str);
    return true;
    // ToDo:
    // Actually, setFromSFZ could also go wrong. This would indicate that the file loading 
    // succeeded but the content of the file could not be parsed (i.e. was malformed or we have a
    // bug in the parser). It could also mean that even though the sfz file itself is ok, we failed
    // to load one or more of the samples - maybe they are not found where they are supposed to be
    // Maybe we should return a return code which could be either of:
    // success, fileLoadError, sfzParseError
  }
  else
    return false;

  // This is clearly not elegant. Get rid of the intermediate c-string!
}

void rsSamplerData::writeSettingToString(const PlaybackSetting& setting, std::string& s)
{
  using PST = Opcode;
  PST  type = setting.getType();
  float val = setting.getValue();
  int index = setting.getIndex();

  // This makes the unit test fail - why?:
  //if(val = PlaybackSetting::getDefaultValue(type))
  //  return; // default values need not to be stored - todo: maybe optionally store them anyway

  // Helper function to add an opcode with value to the string str:
  auto add = [](std::string& str, const char* opcodeName, float value)
  {
    str += opcodeName + std::string("=") + to_string(value) + "\n";
  };
  // get rid - not needed anymore

  SfzOpcodeTranslator* t = SfzOpcodeTranslator::getInstance();
  add(s, t->opcodeToStringC(type), val);



  /*
  switch(type)
  {
  case PST::Volume:          { add(s, "volume", val);  } break;
  case PST::Pan:             { add(s, "pan", val);  } break;

  case PST::PitchKeyCenter:  { add(s, "pitch_keycenter", val);  } break;
  case PST::Transpose:       { add(s, "transpose", val);  } break;
  case PST::Tune:            { add(s, "tune", val);  } break;

  case PST::Delay:           { add(s, "delay", val);  } break;
  case PST::Offset:          { add(s, "offset", val);  } break;

  case PST::FilType:         { add(s, "fil_type", val);  } break;
  case PST::Cutoff:          { add(s, "cutoff", val);  } break;
  case PST::Resonance:       { add(s, "resonance", val);  } break;


  // Extensions:
  case PST::DistShape:       { add(s, "dist_shape", val);  } break;
  case PST::DistDrive:       { add(s, "dist_drive", val);  } break;


  default:                  { RAPT::rsError("Unknown Opcode"); }
  }
  */

  // ToDo: 
  // -Maybe use custom string conversion functions because the std::to_string just uses a 
  //  fixed number of 6 decimal digits after the point. Maybe that's suitable, but maybe not:
  //  https://www.cplusplus.com/reference/string/to_string/
  //  ...well, i think, it's not suitable for int params, but we may convert to int. I think, a 
  //  fixed number (maybe 8 or 9..whatever number ensures lossless roundtrips) of total decimal 
  //  digits is better
  // -Why are the lokey, hikey, lovel, hivel opcodes not handled here? I think, it's because they 
  //  are handled already by the caller because they require special treatment. Document this!
}

rsSamplerData::PlaybackSetting rsSamplerData::getSettingFromString(
  const std::string& opcode, const std::string& valStr)
{
  using PS  = PlaybackSetting;
  using PST = Opcode;
  float val = std::stof(valStr);  // maybe use cutom function later
  int   idx = -1;
  // todo: if applicable, exctract the index from the opcode and set it up in the setting by 
  // passing it as 3rd parameter to the constructor

  // Key range:
  if(opcode == "lokey")           return PS(PST::LoKey, val);
  if(opcode == "hikey")           return PS(PST::HiKey, val);
  if(opcode == "lovel")           return PS(PST::LoVel, val);
  if(opcode == "hivel")           return PS(PST::HiVel, val);

  // Amplitude:
  if(opcode == "volume")          return PS(PST::Volume, val);
  if(opcode == "pan")             return PS(PST::Pan,    val);

  // Pitch:
  if(opcode == "pitch_keycenter") return PS(PST::PitchKeyCenter, val);
  if(opcode == "transpose")       return PS(PST::Transpose,      val);
  if(opcode == "tune")            return PS(PST::Tune,           val);
  // todo:  pitch_keytrack, pitch_veltrack, bend_up, bend_down, bend_step 
  // maybe: pitch_random

  // Sample Player:
  if(opcode == "delay")           return PS(PST::Delay,  val);
  if(opcode == "offset")          return PS(PST::Offset, val);

  // todo:  offset, end, count, loop_mode, loop_start, loop_end
  // maybe: delay_random, delay_ccN, offset_random, offset_ccN, sync_beats, sync_offset

  // Filter:
  if(opcode == "fil_type")        return PS(PST::FilType,   val);
  if(opcode == "cutoff")          return PS(PST::Cutoff,    val);
  if(opcode == "resonance")       return PS(PST::Resonance, val);


  // Extensions:
  if(opcode == "dist_shape")      return PS(PST::DistShape, val);
  if(opcode == "dist_drive")      return PS(PST::DistDrive, val);


  // ...more to come...

  return PS(PST::Unknown, 0.f);  // fallback value
}

void rsSamplerData::copy(const rsSamplerData& src, rsSamplerData& dst)
{
  dst.clearInstrument();
  dst.instrument.copyDataFrom(&src.instrument);
  for(int i = 0; i < src.getNumGroups(); i++) {
    const Group* srcGroup = src.getGroup(i);
    Group* dstGroup = new Group;
    dstGroup->copyDataFrom(srcGroup);
    dst.addGroup(dstGroup);
  }
}

}} // namespaces

/*

Bugs:
-when the opcode value in the sfz file is given without decimal dot, e.g. volume=-6 instead of 
 volume=-6.0, it crashes

Notes:
-When writing .sfz files, make sure to use the forward slash "/" as seperator in the sample paths.
 Using the backslash "\" will lead to failure of file loading on mac.

-for the sfz-parsing, support separation between opcodes not only by newline but also by space 
 -any combination of newlines and spaces should be allowed 
 -maybe std::regex could be used for this or: 
 -use a simple state machine with states:
  A: within token, B: within separator, C: token finished
  -start in state A, index = starIndex
  -read symbols until C is reached
   if(state == A && isSeperator(symbol)) {
     state = A; index++; }
   if(state == A && !isSeperator(symbol) {
     state = B; index++ }
   if(state == B && isSeperator(symbol)) {
     state = B; index++ }
   if(state == B && !isSeperator(symbol)) {
     state = C; break; }
   return index;
   where isSeperator(c) returns true, if c is '\n' or ' '

   ...done, i think

ToDo:
-make a unit test that programmatically creates different .sfz files representing the same 
 instrument but with different formatting, render output and compare
-maybe have an inquiry function that takes an opcode and returns the standard in which this 
 opcode is defined (sfz, sfz2, aria, rs, ...)

-In the sfz-spec, it says that the pitch_keycenter opcode can also be specified as e.g. c#4
 -> support that syntax in the sfz parser!

Refactor:
-Maybe make a struct that contains the integer enum index for the opcode (e.g. PST::Pan), the 
 sfz opcode string (e.g. "pan"), the default value (e.g. 0.f). That avoids having to change 
 3 places (PlaybackSetting::getDefaultValue, writeSettingToString, getSettingFromString) when 
 adding a new opcode. It will most likely also reduce the overall amount of code. The code in these
 functions is very repetitive anyway. Then, somewhere in the code we need to have an array of such 
 structs. Disadvantage: The mapping between opcodes indices and their strings must be figured out 
 at runtime (i.e. linear search?)...or we use something like a std::map or some selfmade 
 "dictionary" datastructure...but it needs satellite data for the default value. Maybe a data 
 structure that allows to store an array of pointers to objects (the structs) and additionally 
 maintains an arbitrary number of arrays of indices into that array which are sorted according to
 different criteria, for example, using different fields as key. But where and when would we fill
 that data structure. Maybe it should be some sort of singleton object rsSamplerOpcodeList. Or 
 maybe an object of that class could be a static member of rsSamplerData. It could have member 
 functions:
 rsSamplerData::PlaybackSetting::Type getOpcodeIndex(const::std::string& opcodeString)
 const::std::string& getOpcodeString(...Type opcodeType)
 float getDefaultValue(...Type opcodeType)



*/
