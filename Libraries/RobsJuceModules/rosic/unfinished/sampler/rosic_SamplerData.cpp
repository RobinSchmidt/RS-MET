namespace rosic {
namespace Sampler {

//-------------------------------------------------------------------------------------------------
// The internal classes

float SfzInstrument::PlaybackSetting::getDefaultValue(Opcode type, int index)
{
  SfzCodeBook* t = SfzCodeBook::getInstance();
  return t->opcodeDefaultValue(type, index);
}
DspType SfzInstrument::PlaybackSetting::getTargetDspType(Opcode type)
{
  SfzCodeBook* t = SfzCodeBook::getInstance();
  return t->opcodeToProcessor(type);
}
// Maybe get rid of these two. The caller should use the translator himself. Or maybe it's more 
// convenient to keep them as convenience functions? We'll see....

//-------------------------------------------------------------------------------------------------

void SfzInstrument::HierarchyLevel::ensureDspsPresent(Opcode opcodeType, int howMany)
{ 
  using namespace RAPT;
  using SPT = DspType;
  SPT dspType = PlaybackSetting::getTargetDspType(opcodeType);
  rsAssert( dspType != SPT::Unknown );

  // use: SfzCodeBook::isDspSetting:
  if(dspType == SPT::SamplePlayer || dspType == SPT::Unknown)
    return;
    // The sample-player at the start of the processing chain doesn't really count as bona-fide DSP
    // processor. It's always there, there's always exactly one and it behaves quite differently 
    // from the rest. We need it among the types for consistency, though.

  // Insert however many processors of the required kind are missing to the end of the chain:
  int count   = (int)rsCount(dspTypes, dspType);   // rename to numPresent
  int missing = howMany - count;                   // rename howMany to numRequired
  for(int i = 0; i < missing; i++)
    dspTypes.push_back(dspType);
}

void SfzInstrument::HierarchyLevel::updateDspsArray()
{
  dspTypes.clear();
  for(size_t i = 0; i < settings.size(); i++)
  {
    Opcode op = settings[i].getOpcode();
    int idx   = settings[i].getIndex();
    if(idx != -1)  // opcodes that apply to the DSP array don't have index -1
      ensureDspsPresent(op, RAPT::rsMax(idx, 1));
  }
}

void SfzInstrument::HierarchyLevel::setSetting(const PlaybackSetting& s)
{
  using TP = Opcode;
  TP t = s.getOpcode();

  // Handle the lo/hi key/vel opcodes as special cases:
  if(t == TP::LoKey) { loKey = (uchar)s.getValue(); return; }
  if(t == TP::HiKey) { hiKey = (uchar)s.getValue(); return; }
  if(t == TP::LoVel) { loVel = (uchar)s.getValue(); return; }
  if(t == TP::HiVel) { hiVel = (uchar)s.getValue(); return; }
  // ToDo: maybe we should assert that the value is an integer in the range 0..127
  // we should also handle the "key" opcode which specifies lokey, hikey, 
  // pitch_keycenter simultaneously?

  // The "key" opcode specifies lokey, hikey and pitch_keycenter at the same time:
  if(t == TP::Key)
  {
    loKey = hiKey = (uchar)s.getValue();
    setSetting(PlaybackSetting(Opcode::PitchKeyCenter, s.getValue(), -1));
  }

  // All other settings are handled by either overwriting the last setting of that type in our 
  // array, if present or by appending the setting, if not present:
  int idx = s.getIndex();
  int foundAt = findSetting(t, idx);
  if(foundAt != -1)
    settings[foundAt] = s;
  else
  {
    settings.push_back(s);
    ensureDspsPresent(t, RAPT::rsMax(idx, 1));
    // The order in which the processors appear in the chain should reflect the order in which 
    // their opcodes appear in the sfz (or, if setup is done programmatically, the order in which
    // the opcodes were added). The first opcode applying to a particular kind of processor 
    // counts. For example, if the opcodes are added in the order cutoff, dist_drive, resonance, 
    // the filter appears before the waveshaper in the DSP chain...maybe with some exceptions for 
    // opcodes that apply to processors that must be at fixed positions in the chain such as the 
    // SamplePlayer. ...hmm...but what, if we want two processors of the same kind? like filter1 
    // -> waveshaper -> filter2. -> document behavior in case of indexed DSPs
  }
}

bool SfzInstrument::HierarchyLevel::removeSetting(Opcode type, int index)
{
  bool wasRemoved = false;
  if(index == -1) {
    for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
      if(settings[i].getOpcode() == type) {
        RAPT::rsRemove(settings, i);
        wasRemoved = true; }}}
  else {
    for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
      if(settings[i].getOpcode() == type && settings[i].getIndex() == index) {
        RAPT::rsRemove(settings, i);
        wasRemoved = true; }}}
  updateDspsArray();
  return wasRemoved;
  // We can't use size_t for i because the -1 would create an access violation when size() = 0
  // Maybe it should remove the DSP if it was the last setting that applied to it?
}

void SfzInstrument::HierarchyLevel::copyDataFrom(const HierarchyLevel* lvl)
{
  samplePath = lvl->samplePath;
  settings   = lvl->settings;
  dspTypes   = lvl->dspTypes;

  // not sure, if the pointers should be copied - maybe not:
  //custom = lvl->custom;  // this one may be, it's the pointer to the audio stream
  //parent = lvl->parent;  // no, this one definitely not (i think)
}

float SfzInstrument::HierarchyLevel::getSettingValue(
  Opcode type, int index, bool accumulate) const
{
  float val = PlaybackSetting::getDefaultValue(type, index);  // init to global fallback value
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

int SfzInstrument::HierarchyLevel::findSetting(Opcode type, int index) const
{
  for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
    if(settings[i].getOpcode() == type && settings[i].getIndex() == index)
      return (int)i;
  }
  return -1;
}

void SfzInstrument::Region::copyDataFrom(const Region* src)
{
  SfzInstrument::HierarchyLevel::copyDataFrom(src);
  loKey = src->loKey;
  hiKey = src->hiKey;
  loVel = src->loVel;
  hiVel = src->hiVel;
  int dummy = 0;
}

bool SfzInstrument::Region::operator==(const SfzInstrument::Region& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= dspTypes == rhs.dspTypes;
  equal &= loKey == rhs.loKey;
  equal &= hiKey == rhs.hiKey;
  equal &= loVel == rhs.loVel;
  equal &= hiVel == rhs.hiVel;
  equal &= samplePath == rhs.samplePath;
  return equal;
  // What about the customPointer? should we require that to be equal, too?
}

int SfzInstrument::Group::addRegion(uchar loKey, uchar hiKey)
{
  SfzInstrument::Region* r = new SfzInstrument::Region;
  r->setLoKey(loKey);
  r->setHiKey(hiKey);
  return addRegion(r);
}

int SfzInstrument::Group::addRegion(Region* r)
{
  r->setParent(this);
  regions.push_back(r);
  return ((int)regions.size()) - 1;
}

bool SfzInstrument::Group::removeRegion(int i)
{
  if(i < 0 || i >= (int)regions.size())
    return false;
  delete regions[i];
  RAPT::rsRemove(regions, i);
  return true;
}

void SfzInstrument::Group::copyDataFrom(const Group* src)
{
  SfzInstrument::HierarchyLevel::copyDataFrom(src);
  clearRegions();
  settings = src->getSettings();
  for(int i = 0; i < src->getNumRegions(); i++) {
    const SfzInstrument::Region* srcRegion = src->getRegion(i);
    SfzInstrument::Region* dstRegion = new SfzInstrument::Region;
    dstRegion->copyDataFrom(srcRegion);
    addRegion(dstRegion);
  }
}

void SfzInstrument::Group::clearRegions()
{
  for(size_t i = 0; i < regions.size(); i++)
    delete regions[i];
  regions.clear();
}

int SfzInstrument::Group::getRegionIndex(const SfzInstrument::Region* region) const
{
  for(size_t i = 0; i < regions.size(); i++)
    if(regions[i] == region)
      return (int)i;
  return -1;
}

SfzInstrument::Region* SfzInstrument::Group::getRegion(int i) const
{
  if(i < 0 || i >= (int)regions.size()) {
    RAPT::rsError("Invalid region index");
    return nullptr;
  }
  return regions[i];
}

bool SfzInstrument::Group::operator==(const SfzInstrument::Group& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= dspTypes == rhs.dspTypes;
  equal &= regions.size() == rhs.regions.size();
  if(!equal) return false;
  for(size_t i = 0; i < regions.size(); i++)
    equal &= *(regions[i]) == *(rhs.regions[i]);
  return equal;
}

int SfzInstrument::Global::addGroup()
{
  SfzInstrument::Group* g = new SfzInstrument::Group;
  return addGroup(g);
}

int SfzInstrument::Global::addGroup(SfzInstrument::Group* g)
{
  g->parent = this;
  groups.push_back(g);
  return ((int)groups.size()) - 1;
}

void SfzInstrument::Global::clearGroups()
{
  for(size_t i = 0; i < groups.size(); i++)
    delete groups[i];
  groups.clear();
}

bool SfzInstrument::Global::operator==(const SfzInstrument::Global& rhs) const
{
  bool equal = settings == rhs.settings;
  equal &= dspTypes == rhs.dspTypes;
  equal &= groups.size() == rhs.groups.size();
  if(!equal) return false;
  for(size_t i = 0; i < groups.size(); i++)
    equal &= *(groups[i]) == *(rhs.groups[i]);
  return equal;
}

//-------------------------------------------------------------------------------------------------
// The actual SfzInstrument class:

int SfzInstrument::addRegion(int gi, uchar loKey, uchar hiKey)
{
  if(gi < 0 || gi >= (int)global.groups.size()) {
    RAPT::rsError("Invalid group index");
    return -1; }
  int ri = global.groups[gi]->addRegion(loKey, hiKey);  // region index within its group
  return ri;
}

bool SfzInstrument::removeRegion(int gi, int ri)
{
  if(gi < 0 || gi >= (int)global.groups.size()) {
    RAPT::rsError("Invalid group index");
    return false; }
  return global.groups[gi]->removeRegion(ri);
}

rsReturnCode SfzInstrument::setRegionSetting(int gi, int ri, Opcode type, float value, int index)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->regions[ri]->setSetting(PlaybackSetting(type, value, index));
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::removeRegionSetting(int gi, int ri, Opcode type, int index)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  bool wasRemoved = global.groups[gi]->regions[ri]->removeSetting(type, index);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode SfzInstrument::clearRegionSettings(int gi, int ri)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->regions[ri]->clearSettings();
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::setGroupSetting(int gi, Opcode type, float value, int index)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->setSetting(PlaybackSetting(type, value, index));
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::removeGroupSetting(int gi, Opcode type, int index)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex; }
  bool wasRemoved = global.groups[gi]->removeSetting(type, index);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode SfzInstrument::setInstrumentSetting(Opcode type, float value, int index)
{
  global.setSetting(PlaybackSetting(type, value, index));
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::removeInstrumentSetting(Opcode type, int index)
{
  bool wasRemoved = global.removeSetting(type, index);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

void SfzInstrument::clearAllRegionSettings()
{
  for(size_t gi = 0; gi < global.groups.size(); gi++)
    for(size_t ri = 0; ri < global.groups[gi]->regions.size(); ri++)
      global.groups[gi]->regions[ri]->clearSettings();
}

void SfzInstrument::clearAllGroupSettings()
{
  for(size_t gi = 0; gi < global.groups.size(); gi++)
    global.groups[gi]->clearSettings();
}

void SfzInstrument::clearAllInstrumentSettings()
{
  global.clearSettings();
}

void SfzInstrument::clearAllSettings()
{
  clearAllRegionSettings();
  clearAllGroupSettings();
  clearAllInstrumentSettings();
}
// needs test

std::string SfzInstrument::getAsSFZ() const
{
  auto writeSettingsToString = [](const HierarchyLevel* lvl, std::string& str)
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
  writeSettingsToString(&global, str);
  for(int gi = 0; gi < getNumGroups(); gi++) {
    str += "<group>\n";
    writeSettingsToString(getGroup(gi), str);
    for(int ri = 0; ri < getNumRegions(gi); ri++) {
      str += "<region>\n";
      writeSettingsToString(getRegion(gi, ri), str); }}
  return str;

  // ToDo: write lokey/hikey settings into the string, they are stored directly in the Region 
  // object and not also in the settings. Maybe they should be. That would simplify the 
  // serialization but that may complicate other things due to the introduced redundancy and 
  // therefore extra care to keep the data consistent. That raises the question, if groups and
  // instruments also can define lokey/hikey settings and how they are interpreted. If so, maybe
  // these lokey/hikey members should be moved into the HierarchyLevel baseclass. 
  // ...done ...verify and delete comment

  // ToDo: figure out how SFZPlayer behaves with respect to this maybe by defining those opcodes
  // at all 3 levels - i guess, it will use the most restrictive setting of all of them
}

rsReturnCode SfzInstrument::setFromSFZ(const std::string& strIn) // rename to setFromSfz
{
  clearInstrument();
  if(strIn.empty())
    return rsReturnCode::failed;
  size_t endOfFile = std::numeric_limits<size_t>::max();

  // Pre-process the string to make parsing easier:
  std::string str = strIn;
  rsRemoveLineComments(str, '/');     // remove the comments
  rsReplaceCharacter(str, '\n', ' '); // replace newlines with whitespaces
  rsRemoveRepeats(str, ' ');          // replace sequences of whitespaces with single whitespace
  // -Factor out into a function preProcessSfz
  // -Include stripping away comments. A comment begins with a slash '/' and extends until the end
  //  of the line. It's really annyoing that we can't yet write any comments. Maybe hava a general
  //  function that removes all characters between a startTag and endTag and call it like
  //    str = removeBetween(str, "/", "\n", true, false)
  //  where the true/false flags indicate whether or not the startTag, endTag (here "/" and "\n") 
  //  characters themselves should also be removed. The slash itself shall be removed but the 
  //  newline should remain intact. Of course, it must be called before rsReplaceCharacter.



  // Extracts the subtring starting at startIndex up to (and excluding) the next separator ' ' 
  // charcater. If there is no ' ', it will return the string from startIndex up to its end:
  std::string sep(" ");
  auto getToken = [&](const std::string& str, size_t startIndex)
  {
    int start  = (int)startIndex;
    int length = -1;  // initial value should not matter
    rosic::rsFindToken(str, sep, &start, &length);
    return str.substr(start, length);
  };
  // todo: maybe use a simpler implementation that uses only one single seperator character and
  // assumes that it occurs only once - due to our new pre-processing stage for the string, we
  // can ensure this

  // Sets up one setting in lvl given in the format "opcode=value":
  auto setupSetting = [&](HierarchyLevel* lvl, const std::string& str)
  {
    size_t splitIndex = str.find('=', 0);
    std::string opcode = str.substr(0, splitIndex);
    std::string value  = str.substr(splitIndex+1, str.length() - splitIndex - 1);
    if(opcode == "sample") {     // needs to be treated in a special way
      lvl->setSamplePath(value);
      return;  }
    PlaybackSetting ps = getSettingFromString(opcode, value);
    lvl->setSetting(ps);
  };

  // Sets up the given level according to the given string which is supposed to contain one setting
  // per line in the format "opocde=value\n":
  auto setupLevel = [&](HierarchyLevel* lvl, const std::string& str)
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
  setupLevel(&global, tmp);

  // Loop over the the groups within the instrument definition:
  bool allGroupsDone = false;
  while(!allGroupsDone)
  {
    if(i1 == endOfFile) {
      allGroupsDone = true;
      i1 = str.length(); }

    // Extract substring with group definition and add a new group to the instrument:
    std::string groupDef = str.substr(i0, i1-i0); // group definition (ToDo: use string_view)
    int gi = global.addGroup();
    Group* g = global.getGroup(gi);
    g->parent = &global;

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
        j1 = groupDef.length(); }

      // Extract substring with region definition and add a new region to the group:
      //std::string regionDef = groupDef.substr(j0, j1-j0); // for debug?
      int ri = g->addRegion();
      Region* r = g->getRegion(ri);
      r->setParent(g);

      // Set up region level settings:
      tmp = groupDef.substr(j0+Lr, j1-j0-Lr);
      setupLevel(r, tmp);
    }

    // Find start and end index of next group defintion:
    i0 = str.find(group, i1);      // start index of the group in the string
    i1 = str.find(group, i0+1);    // end index of the group in the string
  }

  return rsReturnCode::success;

  // ToDo: 
  // -There are actually more points of failure which we may have to report. we need to report some
  //  sort of "parseError" or something - maybe we could be more specific about the kind of parse
  //  error
  // -The general structure of the nested region is the similar to the enclosing group block 
  //  -> try to refactor to get rid of the duplication (maybe it can be implemented recursively)
  // -Maybe use string_view for the extracted substrings to avoid copying the data:
  //  https://en.cppreference.com/w/cpp/header/string_view
}

bool SfzInstrument::saveToSFZ(const char* path) const
{
  std::string sfz = getAsSFZ();
  return rsWriteStringToFile(path, sfz.c_str());
}
// this has no safeguards against overwriting an existing file!

rsReturnCode SfzInstrument::loadFromSFZ(const char* path)
{
  // just for debug, to figure out, in which directory the mac expects the sfz file:
  //rsWriteStringToFile("TestFile.sfz", "blablabla");
  // that fails, too with an "Unable to open file" error. Could it have to do with permission?

  char* c_str = rsReadStringFromFile(path);
  if(c_str)
  {
    std::string sfz(c_str);
    rsReturnCode rc = setFromSFZ(sfz);
    free(c_str);
    return rc;

    //return true;
    // ToDo:
    // Actually, setFromSFZ could also go wrong. This would indicate that the file loading 
    // succeeded but the content of the file could not be parsed (i.e. was malformed or we have a
    // bug in the parser). It could also mean that even though the sfz file itself is ok, we failed
    // to load one or more of the samples - maybe they are not found where they are supposed to be
    // Maybe we should return a return code which could be either of:
    // success, fileLoadError, sfzParseError
  }
  else
    return rsReturnCode::fileLoadError;

  // This is clearly not elegant. Get rid of the intermediate c-string!
}

void SfzInstrument::writeSettingToString(const PlaybackSetting& setting, std::string& s)
{
  using PST = Opcode;
  PST    op = setting.getOpcode();
  float val = setting.getValue();
  int   idx = setting.getIndex();  // not yet used but shall be later

  // This makes the unit test fail - why?:
  //if(val = PlaybackSetting::getDefaultValue(type))
  //  return; // default values need not to be stored - todo: maybe optionally store them anyway

  SfzCodeBook* t = SfzCodeBook::getInstance();
  s += t->opcodeToString(op, idx) + std::string("=") + t->valueToString(op, val)  + "\n";

  // ToDo:
  // -Document why the lokey, hikey, lovel, hivel opcodes are not handled here. I think, it's 
  //  because they are handled already by the caller because they require special treatment.
}

SfzInstrument::PlaybackSetting SfzInstrument::getSettingFromString(
  const std::string& opStr, const std::string& valStr)
{
  using PS  = PlaybackSetting;
  using PST = Opcode;
  SfzCodeBook* t = SfzCodeBook::getInstance();
  int idx;
  PST   op  = t->stringToOpcode(opStr, &idx);
  float val = t->stringToValue(op, valStr);
  return PS(op, val, idx);
}

void SfzInstrument::copy(const SfzInstrument& src, SfzInstrument& dst)
{
  dst.clearInstrument();
  dst.global.copyDataFrom(&src.global);
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


*/
