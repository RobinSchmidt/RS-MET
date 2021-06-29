//-------------------------------------------------------------------------------------------------
// The internal classes

float rsSamplerData::PlaybackSetting::getDefaultValue(Type type)
{
  using TP = PlaybackSetting::Type;
  switch(type)
  {
  case TP::LoKey:          return 0.f;
  case TP::HiKey:          return 127.f;
  case TP::LoVel:          return 0.f;
  case TP::HiVel:          return 127.f;

  case TP::Volume:         return 0.f;
  case TP::Pan:            return 0.f;
  case TP::PitchKeyCenter: return 60.f;  // verify!
  }

  RAPT::rsError("Unknown type of PlaybackSetting, i.e. unknown sfz opcode.");
  return 0.f;  // maybe we should return NaN?
}

void rsSamplerData::OrganizationLevel::setSetting(const PlaybackSetting& s)
{
  using TP = PlaybackSetting::Type;
  TP t = s.getType();

  // Handle the lo/hi key/vel opcodes as special cases:
  if(t == TP::LoKey) { loKey = (uchar) s.getValue(); return; }
  if(t == TP::HiKey) { hiKey = (uchar) s.getValue(); return; }
  if(t == TP::LoVel) { loVel = (uchar) s.getValue(); return; }
  if(t == TP::HiVel) { hiVel = (uchar) s.getValue(); return; }
  // ToDo: maybe we should assert that the value is an integer in the range 0..127

  // All other settings are handled by either overwriting the last setting of that type in our 
  // array, if present or appending the setting, if not present:
  int i = findSetting(t);
  if(i != -1) 
    settings[i] = s;
  else
    settings.push_back(s);
}

void rsSamplerData::OrganizationLevel::copyDataFrom(const OrganizationLevel* lvl)
{
  samplePath = lvl->samplePath;
  settings   = lvl->settings;

  // not sure, if the pointers should be copied - maybe not:
  //custom = lvl->custom;
  //parent = lvl->parent;
}

int rsSamplerData::OrganizationLevel::findSetting(PlaybackSetting::Type type, int index)
{
  for(int i = ((int)settings.size()) - 1; i >= 0; i--) {
    if(settings[i].getType() == type && settings[i].getIndex() == index)
      return (int) i; }
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
  r->parent = this;
  regions.push_back(r);
  return ((int) regions.size()) - 1;
}

bool rsSamplerData::Group::removeRegion(int i)
{
  if(i < 0 || i >= (int) regions.size())
    return false;

  RAPT::rsError("Not yet implemented"); // something to do

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
    addRegion(dstRegion); }
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
      return (int) i;
  return -1;
}

rsSamplerData::Region* rsSamplerData::Group::getRegion(int i) const
{
  if(i < 0 || i >= (int)regions.size()) {
    RAPT::rsError("Invalid region index");
    return nullptr; }
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
  return ((int) groups.size()) - 1;
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
    return -1; }
  int ri = instrument.groups[gi]->addRegion(loKey, hiKey);  // region index within its group
  return ri;
}

bool rsSamplerData::removeRegion(int gi, int ri)
{
  if(gi < 0 || gi >= (int)instrument.groups.size()) {
    RAPT::rsError("Invalid group index");
    return false; }
  return instrument.groups[gi]->removeRegion(ri);
}

std::string rsSamplerData::getAsSFZ() const
{
  auto writeSettingsToString = [](const OrganizationLevel* lvl, std::string& str)
  {
    auto toStr = [](const uchar c) { return std::to_string(c); }; // uchar to string
    const std::string& samplePath = lvl->getSamplePath();
    if(!samplePath.empty()   ) str += "sample=" + samplePath + '\n';
    if(lvl->getLoKey() !=   0) str += "lokey=" + toStr(lvl->getLoKey()) + '\n';
    if(lvl->getHiKey() != 127) str += "hikey=" + toStr(lvl->getHiKey()) + '\n';
    if(lvl->getLoVel() !=   0) str += "lovel=" + toStr(lvl->getLoVel()) + '\n';
    if(lvl->getHiVel() != 127) str += "hivel=" + toStr(lvl->getHiVel()) + '\n';
    using SettingsRef = const std::vector<PlaybackSetting>&;
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
      writeSettingsToString(getRegion(gi, ri), str); }}
  return str;

  // ToDo: write lokey/hikey settings into the string, they are stored directly in the Region 
  // object and not also in the settings. Maybe they should be. That would simplify the 
  // serialization but that may complicate other things due to the introduced redundancy and 
  // therefore extra care to keep the data consistent. That raises the question, if groups and
  // instruments also can define lokey/hikey settings and how they are interpreted. If so, maybe
  // these lokey/hikey members should be moved into the OrganizationLevel baseclass. ...done 
  // ToDo: figure out how SFZPlayer behaves with respect to this maybe by defining those opcodes
  // at all 3 levels - i guess, it will use the most restrictive setting of all of them
}

void rsSamplerData::setFromSFZ(const std::string& str)
{
  clearInstrument();
  size_t endOfFile = std::numeric_limits<size_t>::max();

  // Extracts the subtring starting at startIndex up to (and excluding) the next newline '\n' 
  // charcater. If there is no '\n', it will return the string from startIndex up to its end:
  std::string sep(" \n");  // allowed seperator characters
  auto getToken = [&](const std::string& str, size_t startIndex)
  {
    int start  = (int) startIndex;
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
      return;  }

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
      i1 = str.length() - 1; }

    // Extract substring with group definition and add a new group to the instrument:
    std::string groupDef = str.substr(i0, i1-i0); // group definition (todo: use string_view)
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
      j1 = groupDef.find(region, j0+1);
      if(j1 == endOfFile) {
        allRegionsDone = true;
        j1 = groupDef.length() - 1; }

      // Extract substring with region definition and add a new region to the group:
      std::string regionDef = groupDef.substr(j0, j1-j0); // region definition (todo: use string_view)
      int ri = g->addRegion();
      Region* r = g->getRegion(ri);
      r->parent = g;

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
  char* c_str = rsReadStringFromFile(path);
  if(c_str)
  {
    std::string sfz(c_str);
    setFromSFZ(sfz);
    free(c_str);
    return true; 
    // Actually, setFromSFZ could also go wrong. This would indicate that the file loading 
    // succeeded but the content of the file could not be parsed (i.e. was malformed or we have a
    // bug in the parser). Maybe we should return a return code which could be either of:
    // success, fileLoadError, sfzParseError
  }
  else
    return false;

  // This is clearly not elegant. Get rid of the intermediate c-string!
}

void rsSamplerData::writeSettingToString(const PlaybackSetting& setting, std::string& s)
{
  using PST = PlaybackSetting::Type;
  PST  type = setting.getType();
  float val = setting.getValue();
  int index = setting.getIndex();

  // This makes the unit test fail - why?:
  //if(val = PlaybackSetting::getDefaultValue(type))
  //  return; // default values need not to be stored - todo: maybe optionally store them anyway

  switch(type)
  {
  case PST::Volume:         { s += "volume="          + to_string(val) + "\n";  } break;
  case PST::Pan:            { s += "pan="             + to_string(val) + "\n";  } break;
  case PST::PitchKeyCenter: { s += "pitch_keycenter=" + to_string(val) + "\n";  } break;
    // more to come....
  }

  // todo: 
  // -Maybe use custom string conversion functions because the std::to_string just uses a 
  //  fixed number of 6 decimal digits after the point. Maybe that's suitable, but maybe not:
  //  https://www.cplusplus.com/reference/string/to_string/
  //  ...well, i think, it's not suitable for int params, but we may convert to int. I think, a 
  //  fixed number (maybe 8 or 9..whatever number ensures lossless roundtrips) of total decimal 
  //  digits is better
}

rsSamplerData::PlaybackSetting rsSamplerData::getSettingFromString(
  const std::string& opcode, const std::string& valStr)
{
  using PS  = PlaybackSetting;
  using PST = PS::Type;
  float val = std::stof(valStr);  // maybe use cutom function later
  int   idx = -1;
  // todo: if applicable, exctract the index from the opcode and set it up in the setting by 
  // passing it as 3rd parameter to the constructor

  if(opcode == "lokey")          return PS(PST::LoKey,         val);
  if(opcode == "hikey")          return PS(PST::HiKey,         val);
  if(opcode == "lovel")          return PS(PST::LoVel,         val);
  if(opcode == "hivel")          return PS(PST::HiVel,         val);

  if(opcode == "volume")          return PS(PST::Volume,         val);
  if(opcode == "pan")             return PS(PST::Pan,            val);

  if(opcode == "pitch_keycenter") return PS(PST::PitchKeyCenter, val);
  // ...more to come...

  return PS(PST::Unknown, 0.f);  // fallback value
}
// todo: implement writeToSFZ, loadSFZ (taking filenames as parameters)

void rsSamplerData::copy(const rsSamplerData& src, rsSamplerData& dst)
{
  dst.clearInstrument();
  for(int i = 0; i < src.getNumGroups(); i++) {
    const Group* srcGroup = src.getGroup(i);
    Group* dstGroup = new Group;
    dstGroup->copyDataFrom(srcGroup);
    dst.addGroup(dstGroup); }
}

/*

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



*/