using namespace RAPT;  // try to get rid

int rsSamplerData::addRegion(int gi, uchar loKey, uchar hiKey)
{
  if(gi < 0 || gi >= (int)instrument.groups.size()) {
    rsError("Invalid group index");
    return -1; }
  int ri = instrument.groups[gi]->addRegion(loKey, hiKey);  // region index within its group
  return ri;
}

void rsSamplerData::OrganizationLevel::copyDataFrom(const OrganizationLevel* lvl)
{
  samplePath = lvl->samplePath;
  settings   = lvl->settings;

  // not sure, if the pointers should be copied - maybe not:
  //custom = lvl->custom;
  //parent = lvl->parent;
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

void rsSamplerData::Group::copyDataFrom(const Group* src)
{
  rsSamplerData::OrganizationLevel::copyDataFrom(src);
  clearRegions();
  settings = src->getSettings();
  for(size_t i = 0; i < src->getNumRegions(); i++) {
    const rsSamplerData::Region* srcRegion = src->getRegion((int)i);
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
    rsError("Invalid region index");
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


std::string rsSamplerData::getAsSFZ() const
{
  std::string str;

  using SettingsRef = const std::vector<PlaybackSetting>&;
  auto writeSettingsToString = [](SettingsRef settings, std::string& str)
  {
    for(size_t i = 0; i < settings.size(); i++)
      writeSettingToString(settings[i], str);

    // todo:
    // -If the settings define a sample opcode, write that into the string also - it hink, the 
    //  sample opcode needs special handling because its value is a string, not corresponidng to 
    //  any enum value...what, if we later introduce more such string-valued opcodes, such as the
    //  sample-directory? ...hmm...currently, the sample opcode is treated specially below - maybe 
    //  this should be done here by means of not passing the settings but a reference or pointer
    //  to the OrganizationLevel
    // -Maybe factor that function out into a static member function
    // -Maybe factor out a function settingToString...maybe that should be a member of 
    //  PlaybackSetting. toString/fromString
  };


  writeSettingsToString(instrument.getSettings(), str);
  for(size_t gi = 0; gi < getNumGroups(); gi++)
  {
    str += "<group>\n";
    const Group& g = getGroupRef(gi);
    writeSettingsToString(g.getSettings(), str);
    for(size_t ri = 0; ri < g.getNumRegions(); ri++)
    {
      str += "<region>\n";
      const Region* r = getRegionPtr(gi, ri);
      const std::string& samplePath = r->getSamplePath();
      if(!samplePath.empty())
        str += "sample=" + samplePath + '\n';
      writeSettingsToString(r->getSettings(), str);
    }
  }
  return str;

  // ToDo: write lokey/hikey settings into the string, they are stored directly in the Region 
  // object and not also in the settings. Maybe they should be. That would simplify the 
  // serialization but that may complicate other things due to the introduced redundancy and 
  // therefore extra care to keep the data consistent. That raises the question, if groups and
  // instruments also can define lokey/hikey settings and how they are interpreted. If so, maybe
  // these lokey/hikey members should be moved into the OrganizationLevel baseclass. ToDo: figure
  // out how SFZPlayer behaves with respect to this maybe by defining those opcoded at all 3 
  // levels - i guess, it will use the most restrictive setting of all of them
}

void rsSamplerData::setFromSFZ(const std::string& str)
{
  clearInstrument();
  size_t endOfFile = std::numeric_limits<size_t>::max();

  // Extracts the subtring starting at startIndex up to (and excluding) the next newline '\n' 
  // charcater. If there is no '\n', it will return the string from startIndex up to its end:
  auto getLine = [&](const std::string& str, size_t startIndex)
  {
    size_t endIndex = str.find('\n', startIndex);
    if(endIndex >= str.length())
      return str.substr(startIndex, str.length()-startIndex);
    else
      return str.substr(startIndex, endIndex-startIndex);
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
    lvl->addSetting(ps);  // todo: use setSetting (add or overwrite)
  };

  // Sets up the given level according to the given string which is supposed to contain one setting
  // per line in the format "opocde=value\n":
  auto setupLevel = [&](OrganizationLevel* lvl, const std::string& str)
  {
    size_t start = 0;
    while(true)
    {
      std::string line = getLine(str, start);   // extract one line at at time
      if(line.length() == 0) break;
      setupSetting(lvl, line);                  // set a setting from this line
      start += line.length() + 1;
      if(start >= str.length()) break;
    }
  };


  std::string group  = "<group>\n";   // not sure, whether we should include the \n
  std::string region = "<region>\n";
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



// Low-level helper functions for file I/O:
bool rosic::rsWriteStringToFile(const char* path, const char* str)
{
  // ToDo: provide optional argument for the length, defaults to 0 in which case we use strlen

  FILE* f = fopen(path, "w");
  if(f != NULL){
    fwrite(str, 1, strlen(str), f);
    fclose(f); 
    return true; }
  else {
    rsError("Unable to open file");
    return false; }

  // https://www.tutorialspoint.com/cprogramming/c_file_io.htm
}

char* rosic::rsReadStringFromFile(const char *filename)
{
  char *buffer = NULL;
  //int string_size, read_size;  // original
  size_t string_size, read_size;
  FILE *handler = fopen(filename, "r");
  if(handler)
  {
    fseek(handler, 0, SEEK_END);  // seek the last byte of the file
    string_size = ftell(handler); // offset from the first to the last byte, filesize
    rewind(handler);              // go back to the start of the file
    buffer = (char*) malloc(sizeof(char) * (string_size + 1) );    // allocate a string
    read_size = fread(buffer, sizeof(char), string_size, handler); // read it all in one go
    //buffer[string_size] = '\0';   // this was the original code
    buffer[read_size] = '\0';   // put a \0 in the last position



    //if(string_size != read_size)
    //{
    //  // Something went wrong, free memory and set the buffer to NULL
    //  free(buffer);
    //  buffer = NULL;
    //}
    // We actually run into this branch, apparently due to CR/LF line-ending stuff. In the unit 
    // test, the file has 7 lines and 116 characters in total with CR/LF line endings, but the 
    // function only reads 109 characters. The extra 7 characters apparently are the additional
    // line ending symbols


    fclose(handler);    // close the file
  }
  return buffer;

  // Code adapted from here:
  // https://stackoverflow.com/questions/3463426/in-c-how-should-i-read-a-text-file-and-print-all-strings
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
    // actually, setFromSFZ could also go wrong - we should do int rc = setFromSFZ and return rc
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
  for(size_t i = 0; i < src.getNumGroups(); i++) {
    const Group* srcGroup = src.getGroupPtr(i);
    Group* dstGroup = new Group;
    dstGroup->copyDataFrom(srcGroup);
    dst.addGroup(dstGroup); }
}