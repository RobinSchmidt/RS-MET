//=================================================================================================
// Function definitions for the helper classes:

template<class T>
bool AudioFileStreamPreloaded<T>::setData(
  T** newData, int numFrames, int numDataChannels, T sampleRate, int numStreamChannels, 
  const std::string& path)
{
  // Deallocate old and allocate new memory:
  clear();
  int numChannelsMin = rsMin(numDataChannels, numStreamChannels);
  flatData = new T[numChannelsMin*numFrames];
  channelPointers = new T*[numStreamChannels];
  if(flatData == nullptr || channelPointers == nullptr) {
    clear(); return false; }  // memory allocation failed

  // Copy the new data into the freshly allocated memory:
  //for(int c = 0; c < numChannels; c++) {
  //  channelPointers[c] = &flatData[c*numFrames];
  //  for(int n = 0; n < numFrames; n++)
  //    channelPointers[c][n] = newData[c][n]; }
  // Maybe we should have a version of this function which does not need to copy data but instead
  // just takes over ownership of the passed array. But this would need a parameter for the flat 
  // data array, too. We'll see, how this meshes with the wavefile loading functions...


  for(int c = 0; c < numChannelsMin; c++) 
    for(int n = 0; n < numFrames; n++)
      flatData[c*numFrames + n] = newData[c][n];
  for(int c = 0; c < numStreamChannels; c++)
    channelPointers[c] = &flatData[c % numChannelsMin];

  // Update metadata members and report success:
  this->numChannels = numStreamChannels;
  this->numFrames   = numFrames;
  this->sampleRate  = sampleRate;
  this->path        = path;
  return true; // success
}

/*
void AudioFileStreamPreloaded::setNumOutputChannels(int newNumChannels)
{
  int dummy = 0;
}
*/

template<class T>
void AudioFileStreamPreloaded<T>::clear()
{
  numChannels = 0;
  numFrames   = 0; 

  delete[] channelPointers;
  channelPointers = nullptr; 

  delete[] flatData;
  flatData = nullptr;
}

//-------------------------------------------------------------------------------------------------

template<class T>
int SamplePool<T>::findSample(const std::string& path) const
{
  for(size_t i = 0; i < samples.size(); i++) {
    if(samples[i]->getPath() == path)
      return (int) i; }
  return -1;
}

template<class T>
void SamplePool<T>::clear()
{
  for(size_t i = 0; i < samples.size(); i++)
    delete samples[i];
  samples.clear();
}

//-------------------------------------------------------------------------------------------------
// rsSamplerData:

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
bool rsWriteStringToFile(const char* path, const char* str)
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

char* rsReadStringFromFile(const char *filename)
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

//=================================================================================================
// rsSamplerEngine

rsSamplerEngine::rsSamplerEngine(int maxNumLayers)
{
  // factor out into setMaxNumLayers
  int L = maxNumLayers;
  playerPool.resize(L);
  idlePlayers.resize(L);
  activePlayers.reserve(L);
  for(int i = 0; i < L; i++)
    idlePlayers[i] = &playerPool[i];
}

rsSamplerEngine::~rsSamplerEngine()
{

}

//-------------------------------------------------------------------------------------------------
// Setup:

int rsSamplerEngine::addSampleToPool(
  float** data, int numFrames, int numChannels, float sampleRate, const std::string& path)
{
  // todo: 
  // -check, if a sample with the same path already exists - if so, we have nothing to 
  //  do and may return early with an appropriate code
  //   if(isSampleInPool(..)) return ReturnCode::nothingToDo;
  // -currently, we do such a check in addSamplesUsedIn - maybe it should be done here instead 
  //  and/or in loadSampleToPool

  AudioFileStreamPreloaded<float>* stream = new AudioFileStreamPreloaded<float>;
  bool allocOK = stream->setData(data, numFrames, numChannels, sampleRate, 2, path);
  if(allocOK == false)
    return ReturnCode::memAllocFail;
  return samplePool.addSample(stream);
}

int rsSamplerEngine::loadSampleToPool(const std::string& path)
{
  if(isSampleInPool(path))
    return ReturnCode::nothingToDo;

  int numFrames;
  int numChannels;
  int sampleRate;
  float** data = rosic::readFloatFromWaveFile(path.c_str(), numChannels, numFrames, sampleRate);
  if(data == nullptr)
    return ReturnCode::fileLoadError;
    // This could mean that the file was not found or the memory allocation failed. ToDo: be more 
    // specific which of the two conditions happened

  int rc = addSampleToPool(data, numFrames, numChannels, (float) sampleRate, path);

  // Clean up data (todo: wrap into utility function - there is actually already one):
  for(int c = 0; c < numChannels; c++)
    delete[] data[c];
  delete[] data;

  return rc;
}

int rsSamplerEngine::unUseSample(int i)
{
  if(i < 0 || i >= samplePool.getNumSamples()){
    rsError("Invalid sample index");
    return ReturnCode::invalidIndex; }
  int numRegions = 0;
  using StreamPtr = const AudioFileStream<float>*;
  StreamPtr stream = samplePool.getSampleStream(i);
  for(int gi = 0; gi < getNumGroups(); gi++) {
    for(int ri = 0; ri < getNumRegions(gi); ri++) {
      Region* r = getRegion(gi, ri);
      StreamPtr regionStream = (StreamPtr) r->getCustomPointer();
      if(regionStream == stream) {
        r->setCustomPointer(nullptr);
        r->setSamplePath("");
        numRegions++;  }}}
  return numRegions;

  // todo: maybe try to implement a sort of "forAllRegions" macro that implements the boildplate to
  // loop over all regions once and for all
}

int rsSamplerEngine::unUseSample(const std::string& samplePath)
{
  int i = findSampleIndexInPool(samplePath);
  if(i == -1) return 0;
  return unUseSample(i);
}

// todo: implement removeSample - should call unUseSample before actually removing the sample from 
// the pool to avoid dangling pointers

int rsSamplerEngine::addRegion(int gi, uchar loKey, uchar hiKey)
{
  int ri = sfz.addRegion(gi, loKey, hiKey);
  if(ri == -1)
    return ReturnCode::invalidIndex;    // gi was an invalid group index
  Region* r = getRegion(gi, ri);
  for(uchar k = loKey; k <= hiKey; k++)
    addRegionForKey(k, r);
  return ri;
}

int rsSamplerEngine::setRegionSample(int gi, int ri, int si)
{
  if(!isIndexPairValid(gi, ri)) {
    rsError("Invalid group- and/or region index");
    return ReturnCode::invalidIndex; }
  if(!isSampleIndexValid(si)) {
    rsError("Invalid sample index");
    return ReturnCode::invalidIndex; }

  const AudioFileStream<float>* s = samplePool.getSampleStream(si);
  sfz.setRegionCustomPointer(gi, ri, (void*) s);
  sfz.setRegionSample(gi, ri, s->getPath());
  return ReturnCode::success;
}

int rsSamplerEngine::setRegionSetting(int gi, int ri, PlaybackSetting::Type type, float value)
{
  if(!isIndexPairValid(gi, ri)) {
    rsError("Invalid group- and/or region index");
    return ReturnCode::invalidIndex; }

  sfz.setRegionSetting(gi, ri, type, value);
  return ReturnCode::success;
}

int rsSamplerEngine::setupFromSFZ(const rsSamplerData& newSfz)
{
  removeSamplesNotUsedIn(newSfz);     // remove samples that are not needed anymore from memory
  int rc1 = addSamplesUsedIn(newSfz); // load samples that are needed but not yet loaded
  sfz = newSfz;                       // replace old sfz instrument definition member with new
  int rc2 = setupAudioStreams();      // connect regions in new sfz with appropriate stream objects
  setupRegionsForKey();               // updates regionsForKey array
  if(rc1 >= 0 && rc2 == ReturnCode::success)
    return ReturnCode::success;
  else
    return ReturnCode::fileLoadError;  
    // ToDo: be more specific about the error condition

  // Maybe have state variables numSamplesAdded, numSamplesRemoved, numSamplesNotFound that can be 
  // inquired from client
  // code. this is useful mainly for unit tests but maybe it makes sense to display such info on 
  // the GUI, too, so the user has some feedback about what patches have a lot of samples in common


  // This function needs to be mutexed with anything that accesses the stream-pointers in the sfz
  // objects
}

int rsSamplerEngine::loadFromSFZ(const char* path)
{
  rsSamplerData newSfz;
  bool wasLoaded = newSfz.loadFromSFZ(path);
  if(!wasLoaded)
    return ReturnCode::fileLoadError;
  return setupFromSFZ(newSfz);
}


//-------------------------------------------------------------------------------------------------
// Inquiry:

rsSamplerEngine::Region* rsSamplerEngine::getRegion(int groupIndex, int regionIndex)
{
  int gi = groupIndex, ri = regionIndex;
  if(gi < 0 || gi >= (int)sfz.instrument.groups.size()) {
    rsError("Invalid group index");
    return nullptr; }
  return sfz.instrument.groups[gi]->getRegion(ri);
}

int rsSamplerEngine::getNumRegionsUsing(int i) const
{
  if(i < 0 || i >= samplePool.getNumSamples()){
    rsError("Invalid sample index");
    return ReturnCode::invalidIndex; }
  int numRegions = 0;
  using StreamPtr = const AudioFileStream<float>*;
  StreamPtr stream = samplePool.getSampleStream(i);
  for(int gi = 0; gi < getNumGroups(); gi++){
    for(int ri = 0; ri < getNumRegions(gi); ri++) {
      const Region* r = getRegionConst(gi, ri);
      StreamPtr regionStream = (StreamPtr) r->getCustomPointer();
      if(regionStream == stream)
        numRegions++; }}
  return numRegions;

  // todo: maybe try to implement a sort of "forAllRegions" macro that implements the boildplate to
  // loop over all regions once and for all
}

int rsSamplerEngine::getNumRegionsUsing(const std::string& samplePath) const
{
  int i = findSampleIndexInPool(samplePath);
  if(i == -1) 
    return 0;
  return getNumRegionsUsing(i);
}

int rsSamplerEngine::findSampleIndexInPool(const std::string& sample) const
{
  return samplePool.findSample(sample);

  /*
  for(int i = 0; i < samplePool.getNumSamples(); i++)
  {
    //const std::string& poolSample = samplePool.getSamplePath(i);
    std::string poolSample = samplePool.getSamplePath(i);
    if(sample == poolSample)
      return i;
  }
  return -1;
  */
  // ToDo: 
  // -factor out the implementation into the SamplePool class, so we can just do:
  //  return samplePool.findSample(sample) 
  // -Maybe keep the samplePool sorted, so we can use binary search here. That requires to
  //  reorder it when new samples are added, but addition of new samples is costly anyway due to 
  //  disk access, so that probably doesn't really matter.
}

//-------------------------------------------------------------------------------------------------
// Processing:

void rsSamplerEngine::processFrame(double* left, double* right)
{
  rsFloat64x2 out = 0.0;
  for(size_t i = 0; i < activePlayers.size(); i++)
  {
    out += activePlayers[i]->getFrame();
    if(activePlayers[i]->hasFinished()) {
      deactivateRegionPlayer(i);
      i--; }
  }
  *left  = out[0];
  *right = out[1];

  // ToDo: Test, if it's more efficient to loop through the activePlayers array backwards. But then
  // we would need a signed int as loop index to make it work also when size() == 0 because we 
  // would need to start at size()-1, which would be around 2^64 for size_t. Then, the i-- could be
  // removed, but that's not the main point. The main point is that the deactivation/removal would 
  // need less data copying.
}

void rsSamplerEngine::processFrame(float* left, float* right)
{
  double L, R;
  processFrame(&L, &R);
  *left  = (float) L;
  *right = (float) R;
}

void rsSamplerEngine::processBlock(float** block, int numFrames)
{

}

void rsSamplerEngine::handleMusicalEvent(const rsMusicalEvent<float>& ev)
{
  using Type = rsMusicalEvent<float>::Type;
  Type  type = ev.getType();
  float val1 = ev.getValue1();
  float val2 = ev.getValue2();
  switch(type)
  {
  case Type::noteOn:  handleNoteOn( (uchar)val1, (uchar)val2); break;
  case Type::noteOff: handleNoteOff((uchar)val1, (uchar)val2); break;
  // Maybe we should pass the floats directly to noteOn/Off. This would allow for microtuning. 
  // However, it complicates identifying to which notes a noteOff should apply. Or rather, it's 
  // not really complicated but requires a design decision: should we require an exact match or 
  // allow the range key +-0.5? For the moment, this here is good enough, though.

  }
}

int rsSamplerEngine::stopAllPlayers()
{
  int numPlayers = (int) activePlayers.size();
  for(int i = numPlayers-1; i >= 0; i--)        // just seems nicer to stop them in reverse order
    idlePlayers.push_back(activePlayers[i]);
  activePlayers.clear();
  return numPlayers;
}


//-------------------------------------------------------------------------------------------------
// Internal:

bool rsSamplerEngine::shouldRegionPlay(const Region* r, uchar key, uchar vel)
{
  // Check, if key/vel are in the required ranges for the region. When this function was called 
  // from handleNoteOn, the first check is actually redundant because if key is not within the 
  // key-range of r, r is not supposed to be in regionsForKey[key]. However, we do the check 
  // anyway, anticipating that the function might get called from somwhere else, too):
  if(key < r->loKey || key > r->hiKey) return false;
  if(vel < r->loVel || vel > r->hiVel) return false;


  // Check, if all other playback constraints for the given region as defined in r->settings are 
  // satisfied. If any of them isn't, return false:
  // ...more to do...

  // -to do this efficiently, we should have the settings sorted in some way - the most often
  //  accessed settings should come first. or last. so we can have two sets of most often accessed
  //  settings - maybe those which are accessed often from the audio-rendering should come first 
  //  and those accessed often from the event-processing last. settings that are accessed rarely
  //  go into the middle section. ..but maybe we shoud avoid accessing the settings from the 
  //  rendering anyway


  // The region has passed all constraint filters and should indeed be played:
  return true;
}

void rsSamplerEngine::addRegionForKey(uchar k, const Region* r)
{
  if(r->shouldPlayForKey(k))
    regionsForKey[k].addRegion(r);
  // What, if the region is already there? We should check that before. It's probably not supposed
  // to happen, but anyway. Well...maybe it is, when the user tweaks loKey/hiKey settings on a GUI.
  // On the other hand, it may actually be useful to be able to duplicate regions, especially on
  // a GUI: duplicate-and-edit could be a common workflow
}

void rsSamplerEngine::findRegion(const rsSamplerEngine::Region* r, int* gi, int* ri)
{
  *gi = -1;
  *ri = -1;
  for(size_t i = 0; i < sfz.instrument.groups.size(); i++) {
    int j = sfz.instrument.groups[i]->getRegionIndex(r);
    if(j != -1) {
      *gi = (int) i;
      *ri = j;
      return; }}
  rsError("Region not found");
  // A region should always be found in one (and only one group). If we don't find it, it means
  // the caller has passed a pointer to a region object that is not part of this instrument. If 
  // this happens, it indicates a bug at the call site.
}

rsSamplerEngine::RegionPlayer* rsSamplerEngine::getRegionPlayerFor(
  const Region* r, uchar key, uchar vel)
{
  // The behavior for this in situations where the region r is already being played by some 
  // voice/player should depend on the playback-mode of the region. If it's in one-shot mode, a new
  // player should be handed and the old one should just continue to play along with the new. 
  // Otherwise, the old player should be re-used. But this may cause clicks. Maybe a new player 
  // should be handed and the old one should be flagged for a quick fade-out (a few milliseconds). 
  // Maybe that feature should be optional, controlled by a retriggerFadeOutTime parameter that is 
  // zero by default (indicating hard, clicking retriggers). Or maybe this fade-out could also be 
  // set per region instead of globally? Check, if sfz may actually have an opcode for that.
  if(idlePlayers.empty())
    return nullptr; // Maybe we should implement more elaborate voice stealing?
  RegionPlayer* rp = rsGetAndRemoveLast(idlePlayers);
  rp->setKey(key);
  rp->setRegionToPlay(r, sampleRate);
  activePlayers.push_back(rp);
  return rp;
}

bool rsSamplerEngine::isSampleUsedIn(
  const AudioFileStream<float>* stream, const rsSamplerData& sfz)
{
  const std::string& streamPath = stream->getPath();
  for(size_t gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group& g = sfz.getGroupRef(gi);
    for(size_t ri = 0; ri < g.getNumRegions(); ri++) {
      Region* r = g.getRegion((int)ri);        // the conversion is unelegant - try to get rid
      const std::string& regionPath = r->getSamplePath();
      if(regionPath == streamPath)
        return true; }}
  return false;
}

int rsSamplerEngine::deactivateRegionPlayer(size_t i)
{
  if(i >= activePlayers.size()) {
    rsError("Invalid player index");
    return ReturnCode::invalidIndex; }
  RegionPlayer* p = activePlayers[i];
  rsRemove(activePlayers, i);
  idlePlayers.push_back(p);
  return ReturnCode::success;
}

const AudioFileStream<float>* rsSamplerEngine::getSampleStreamFor(const Region* r)
{
  //return r->getSampleStream();
  return (const AudioFileStream<float>*) r->getCustomPointer();
  // todo: 
  // -some sanity checks may be appropriate here
  // -maybe, if the region itself has stored a nullptr, we should check, if the enclosing group 
  //  defines a stream and if it also doesn't, check the instrument, i.e. walk up the hierarchy 
  //  ladder for fallback streams
  // -maybe that should be done in getCustomPointer already
}

int rsSamplerEngine::handleNoteOn(uchar key, uchar vel)
{
  if(vel == 0) { return handleNoteOff(key, vel); }

  int numRegions = 0;  // number of regions that were triggered by this noteOn
  for(size_t i = 0; i < regionsForKey[key].getNumRegions(); i++) 
  {
    const Region* r  = regionsForKey[key].getRegion(i);
    if(!shouldRegionPlay(r, key, vel))
      continue;
    RegionPlayer* rp = getRegionPlayerFor(r, key, vel);
    if(rp == nullptr) {
      // Roll back the addition of all players so far to the activePlayers and move them back into
      // the idlePlayers again. We don't really want notes to play with an incomplete set of 
      // samples. It's all or nothing - either all samples for the given key get triggered or none
      // of them:
      for(int j = 0; i < numRegions; j++) {
        rp = rsGetAndRemoveLast(activePlayers);
        idlePlayers.push_back(rp); }
      return ReturnCode::voiceOverload;
    }
    else
    {
      //rp->setKey(key);
      numRegions++;
    }
  }
  return ReturnCode::success;
  // Another possibility for the return value would have been to return the number of layers that
  // have been triggered, but we don't do that because then it would be not quite clear what we
  // should return from noteOff to make the functions somewhat consistent. In noteOff, we could
  // either return the number of released regions or the number of triggered release samples. Both
  // would make just as much sense. So, for unambiguous consistency, we let both just return a 
  // success or failure report. 
}

int rsSamplerEngine::handleNoteOff(uchar key, uchar vel)
{
  // Note-off events may also trigger the playback of regions (note-off samples are a thing)...
  int numRegions = 0;
  for(int i = int(activePlayers.size()) - 1; i >= 0; i--)
  {
    if(activePlayers[i]->getKey() == key)
      deactivateRegionPlayer(size_t(i));
      // ToDo: refine this later: we may not want to immediately stop the player but rather 
      // trigger the release phase and mark for quick fade-out
  }
  return ReturnCode::success; // preliminary

  // ToDo:
  // -Mark them for going into release state, if they have an amp-env or stop them immediately, if
  //  they don't. Maybe as a later refinement, we could also apply a quick fade-out (over a couple
  //  of millisecond) instead of a hard stop. That could be useful also in other contexts, such as
  //  different retrigger modes
  // -later: trigger all regions that should play a release sample for the given key/vel combo
  // -check, if looping forward or backward (implying using int or size_t for i) is more
  //  efficient
}

int rsSamplerEngine::removeSamplesNotUsedIn(const rsSamplerData& sfz)
{
  numSamplesRemoved = 0;
  for(int i = samplePool.getNumSamples()-1; i >= 0; i--) {
    const AudioFileStream<float>* sample = samplePool.getSampleStream(i);
    if(!isSampleUsedIn(sample, sfz)) {
      samplePool.removeSample(i);
      numSamplesRemoved++;
      i--; }}
  return numSamplesRemoved;
}

int rsSamplerEngine::addSamplesUsedIn(const rsSamplerData& sfz)
{
  auto isValidSamplePath = [](const std::string& path) { return !path.empty(); };
  // -maybe make this a member function
  // -maybe add more sophisticated tests (but not too costly stuff to not slow down loading)

  numSamplesLoaded = 0;
  numSamplesFailed = 0;
  bool allOK = true;
  for(size_t gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group& g = sfz.getGroupRef(gi);
    for(size_t ri = 0; ri < g.getNumRegions(); ri++) {
      Region* r = g.getRegion((int)ri);     // the conversion is unelegant - try to get rid
      const std::string& path = r->getSamplePath();
      if(isValidSamplePath(path) && !isSampleInPool(path)) {
        int rc = loadSampleToPool(path);
        if(rc >= 0)
          numSamplesLoaded++;
        else
          numSamplesFailed++; }}}
  if(numSamplesFailed == 0)
    return numSamplesLoaded;
  else
    return ReturnCode::fileLoadError;
}

int rsSamplerEngine::setupAudioStreams()
{
  // Function to connect a Region object with one of our stream objects:
  auto setupStream = [this](Region* r, const std::string& path)
  {
    if(path.empty()) {
      r->setCustomPointer(nullptr);
      return true; }
    int si = findSampleIndexInPool(path);
    if(si == -1)
      return false;
    const AudioFileStream<float>* stream = samplePool.getSampleStream(si);
    r->setCustomPointer(stream);
    return true;
  };
  // ToDo: Maybe change Region* to rsSamplerData::OrganizationLevel*, so we can assign streams to 
  // groups and instruments also.

  bool allOK = true;
  for(size_t gi = 0; gi < sfz.getNumGroups(); gi++) 
  {
    const Group& g = sfz.getGroupRef(gi);
    for(size_t ri = 0; ri < g.getNumRegions(); ri++) 
    {
      Region* r = g.getRegion((int)ri);     // the conversion is unelegant - try to get rid
      const std::string& path = r->getSamplePath();
      allOK &= setupStream(r, path);
    }
  }
  // ToDo: Maybe, if getSamplePath returns the empty string, meaning that a region does not define
  // the sample opcode, assign the stream from the outlying group. If that also doesn't define a
  // sample, use the stream from the instrument.

  if(allOK)
    return ReturnCode::success;
  else
    return ReturnCode::notFound;  // stream was not found for one or more samples
}

void rsSamplerEngine::setupRegionsForKey()
{
  for(uchar k = 0; k < numKeys; k++)
    regionsForKey[k].clear();
  for(size_t gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group& g = sfz.getGroupRef(gi);
    for(size_t ri = 0; ri < g.getNumRegions(); ri++) {
      Region* r = g.getRegion((int)ri);  // the conversion is unelegant - try to get rid
      for(uchar k = 0; k < numKeys; k++)
        addRegionForKey(k, r); }}
  int dummy = 0;
}

//-------------------------------------------------------------------------------------------------
// rsSamplerEngine::RegionPlayer

void rsSamplerEngine::RegionPlayer::setRegionToPlay(
  const rsSamplerEngine::Region* regionToPlay, double fs)
{
  region = regionToPlay;
  stream = getSampleStreamFor(region);
  prepareToPlay(fs);
}

rsFloat64x2 rsSamplerEngine::RegionPlayer::getFrame()
{
  if(sampleTime < 0.0) {             // Negatively initialized sampleTime implements delay.
    sampleTime += 1.0;               // We just increment the time and return 0,0. Actual output
    return rsFloat64x2(0.0, 0.0); }  // will be produced as soon as sampleTime reaches zero.  

  float L, R;                        // left and right output

  //int intTime = round(sampleTime);
  int intTime = (int) floor(sampleTime);
  stream->getFrameStereo(intTime, &L, &R);
  //stream->getFrameStereo(sampleTime, &L, &R);
  // This does implicit conversion of double to int by truncation which amounts to the worst 
  // possible quality resampling scheme. ToDo:
  // -implement linear interpolation and later also better methods (sinc, elephant, ...)
  // -keep sample time as combination of int and float to avoid computation of the fractional part
  //  at each sample and to avoid losing precision for the fractional part when the integer part
  //  is large (thereby eating up significant digits).
  // -the interpolation should probably be handled by the AudioStream class to make it re-usable
  //  also for resampled audio playback in other contexts. Maybe here, we should just call
  //  stream->getFrameStereo(sampleTimeInt, sampleTimeFrac, &L, &R);
  // -it should probably also receive the increment in order to make a decision for time-scaling
  //  the sinc, if necessary for anti-aliasing

  // more stuff to do:
  // -apply pitch envelope and lfo - these should affect (scale?) the effective increment that we 
  //  add to sampleTime - but our increment *member* should not be modified, instead, do something
  //  like sampleTime += increment * incScaler; where incScaler is computed from the pitch 
  //  modifiers. Maybe create a subclass of SignalProcessor called PitchShifter that has just a
  //  dummy callback that just stores the desired shift value and we read it out here
  // -apply the DSP processes

  sampleTime += increment;
  return this->amp * rsFloat64x2(L, R); 
}

void rsSamplerEngine::RegionPlayer::processBlock(rsFloat64x2* y, int N)
{
  for(int n = 0; n < N; n++)
    y[n] = getFrame();
  // preliminary - todo: run the different DSP processes, one after another, over the whole block,
  // using in-place processing, the steps are (in that order)
  // -fill y with pitch envelope (including pitch LFO)
  // -fill y with interpolated raw sample values (or: maybe compute pitch envelope on the fly)
  // -apply filter (maybe the filter envelope can be computed on the fly)
  // -apply amp-envelope
}

bool rsSamplerEngine::RegionPlayer::isPlayable(const Region* region)
{
  bool ok = true;
  ok &= region != nullptr;
  ok &= region->getGroup() != nullptr;
  ok &= region->getGroup()->getInstrument() != nullptr;
  ok &= region->getCustomPointer() != nullptr;           // should point to a stream object
  return ok;
}

void rsSamplerEngine::RegionPlayer::prepareToPlay(double fs)
{
  rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  rsAssert(stream != nullptr);   // Ditto.

  bool ok = buildProcessingChain();
  if(!ok)
  {
    //return rsSamplerEngine::ReturnCode::voiceOverload;  // rename to layerOverload
    // This should actually not happen in therory (as by the sfz spec, and unlimited number of 
    // layers is available), but in practice, it may happen in extreme situations like triggering a
    // whole lot of layers at once or in very short succession while already being close to the 
    // limit, such that we don't have enough pre-allocated players and/or dsp objects available and 
    // the required additional allocation of more is not fast enough.
  }

  // To set up the settings, we call setupDspSettings 3 times to:
  // (1) set up the general instrument-wide settings
  // (2) set up group specific settings (this may override instrument settings)
  // (3) set up region specific settings (this may override group and/or instrument settings)
  resetDspState();        // Needs to be done after building the chain
  resetDspSettings();     // Reset all DSP settings to default values
  setupDspSettings(region->getGroup()->getInstrument()->getSettings(), fs);
  setupDspSettings(region->getGroup()->getSettings(), fs);
  setupDspSettings(region->getSettings(), fs);

  // return rsSamplerEngine::ReturnCode::success;
}

bool rsSamplerEngine::RegionPlayer::hasFinished()
{
  //int numFrames = stream->getNumFrames();
  //int tmp = stream->getNumFrames() - 1;
  if( sampleTime >= stream->getNumFrames() )
    return true;

  // todo: 
  // -check also, if the amplitude envelope has reached its end
  // -hmm - maybe, if we allow the frequency envelope to go negative, we could also move 
  //  backward through the sample, so having reached the end of the stream may not actually be an
  //  appropriate condition. Or maybe, we should allow more general time-warping envelopes. 
  //  We'll see

  return false;
}

void rsSamplerEngine::RegionPlayer::resetDspState()
{
  for(size_t i = 0; i < dspChain.size(); i++)
    dspChain[i]->resetState();
  for(size_t i = 0; i < modulators.size(); i++)
    modulators[i]->resetState();
}

bool rsSamplerEngine::RegionPlayer::buildProcessingChain()
{
  // ToDo: build the chain of DSP processors and the set of modulators and wire everything up, as
  // defined by the region settings...

  return true;  // preliminary
}

void rsSamplerEngine::RegionPlayer::resetDspSettings()
{
  // Initialize all values and DSP objects to default values (maybe factor out):
  amp = 1.0;
  sampleTime = 0.0;
  increment  = 1.0;
  // ...more to do... 
  // reset all processors and modulators
}

void rsSamplerEngine::RegionPlayer::setupDspSettings(
  const std::vector<PlaybackSetting>& settings, double fs)
{
  // Loop through the settings of the region and for each setting that is present, change the 
  // value from its default to the stored value:
  
  using PS = PlaybackSetting;
  using TP = PS::Type;


  double tmp = stream->getSampleRate();
  increment  = tmp/fs;  // or fs/tmp?

  double rootKey = 69.0;
  double amp = 1.0;
  double pan = 0.0;
  int panRule = PlaybackSetting::PanRule::linear;

  for(size_t i = 0; i < settings.size(); i++)
  {

    PlaybackSetting setting = settings[i];
    TP type = setting.getType();
    double val = (double) setting.getValue();
    switch(type)
    {
    // Amp settings:
    case TP::Volume:  { amp      = rsDbToAmp(val); } break;
    case TP::Pan:     { pan      = val;            } break;
    case TP::PanRule: { panRule  = (int)val;       } break;

    // Pitch settings:
    case TP::PitchKeyCenter: { rootKey = val; } break;

    // Filter settings:
      //case TP::FilterCutoff: { flt.setCutoff(val);  } break;

    // Equalizer settings:
    // .....

    }
  }

  double pitchOffset = double(key) - rootKey;
  increment *= pow(2.0, pitchOffset / 12.0);
  //increment *= rsPitchOffsetToFreqFactor(pitchOffset); // faster, less precise - but probably
  // precise enough, when we switch to int+float for representing increments




  // From the computed local amp/pan/panRule variables, compute the amp member (which is
  // a rsFloat64x2)
  double t1, t2;  // temporaries
  switch(panRule)
  {
  case PS::PanRule::linear:
  {
    t1 = (pan/200.0) + 0.5; // -100..+100 -> 0..1
    t2 = 1.0 - t1;
    this->amp = 2.0 * amp * rsFloat64x2(t2, t1);
  } break;
  case PS::PanRule::sinCos:
  {
    rsError("not yet implemented");
  } break;
  }
  // ToDo: support the "width" opcode. Implement it by a M/S matrix - we need two rsFloat64
  // members to represent both rows of it. Then, in the realtime processing call, we let the stream
  // produce the stereo pair, compute the product of that pair with both rows by element-wise 
  // multiply and summing the whole vector and then assigning the first output to first sum and
  // the 2nd to the 2nd. But should width be applied before or after the pan? -> test with 
  // sfzplayer. I would say, width-before-pan makes more sense from a usability perspective. If 
  // sfz thinks otherwise, maybe provide both options, switched by an additional opcode. sfz also 
  // has the position opcode. maybe that's post-width and apn pan is pre-width?

  // ToDo:
  // -Maybe within the switch statement set up some flags that indicate, if a particular setting is
  //  used. If the flag is false, we may skip the associated DSP process in getFrame/processBlock. 
  //  We may need inquiry functions such as hasFilter, hasAmpEnv, hasPitchEnv, hasFilterEnv, 
  //  hasPitchLFO. But this makes things more complicated, so maybe it's not really a good idea.
}

//=================================================================================================
/*
void rsSamplerEngineLoaderSFZ::setFromSFZ(rsSamplerEngine* se, const std::string& str)
{


}

std::string rsSamplerEngineLoaderSFZ::getAsSFZ(const rsSamplerData& sfz)
{
  std::string str;

  // todo: 
  // -retrieve settings of the instrument
  // -write them into the string (using a function writeSettingsToString(settings, str)
  // -loop over the groups
  //  -write a group header
  //  -retrieve group settings
  //  -write them into the string
  //  -loop over the regions
  //   -write a region header
  //   -retrieve region settings
  //   -write them into the string

  return str;
}
*/
// implement writeToSFZ, loadSFZ (taking filenames as parameters)




//=================================================================================================

/*

Goals: 
-Implement (a subset of) the feature set of the sfz specification, perhaps with some extensions 
 that are specifically necessary for the drum sampler. The general architecture should be such 
 that it will possible (and hopefully easy) to implement the full feature set of sfz (and sfz2?) 
 later.
-It should be able to parse sfz files and set itself up accordingly. But maybe that should go into
 a separate class. It should also be able to export its settings to sfz. Maybe make a class
 rsSamplerEngineSFZ. Should also warn when features from an imported sfz file are not supported and
 when the setup is not representable by an sfz file when exporting to sfz.
-It should support different ways of streaming the audio - at least: 
   (1) all samples are preloaded into RAM
   (2) direct-from-disk streaming (DFD)
 To implement this, it should use some sort of suitable abstraction of an audio-streamer that can 
 be plopped in. At first, we implement just the (much simpler) preloading version but it should be
 straightforward to add DFD later.
-In th sfz spec, the performance parameters defined for the whole instrument and for groups work as
 fallback values that can be overriden by respective values on a lower level of the hierarchy. The
 drum-sampler needs them to work in an accumulative fashion. For example a cutoff defined for a 
 region should be used as is and the cutoff for the group should be for a 2nd filter through which 
 the signal of the whole group is passed. That's a different semantic. In the original sfz-spec, 
 the group cutoff would just have been overriden by the region setting (i think -> verify). Maybe 
 we should have switches like: group/instrumentSettingsAccumulate and/or maybe that should be done
 in a subclass
-Loop-points and start/end points shall be floating point numbers. When parsing an sfz file, we 
 just split the number string into pre-dot and post-dot part, parse them separately and store the 
 post-dot part in a double (or TPar) and the pre-dot part in an integer. That way, we won't suffer 
 precision loss for numbers bigger pre-dot part. This should be downward compatible with sfz spec 
 (which has only integer loop points (i think -> verify))

ToDo:
-The GUI should show a warning message, when the maximum number of voices is exceeded. SFZ 
 specifies a practically infinite number of voices, so in order to be compliant to the spec, we 
 should always have enough voices available. The original SFZ.exe by rgcaudio has 256 voices, 
 where the term voice refers to a single region here, so in this terminology, a signle key can 
 already play multiple voice due to key- and velocity crossfades, for example.
-Allow the user to select between storing samples as float or short int (16 Bit) in memory. Allows
 different trade-off between memory and cpu usage: float needs twice as much memory but doesn't 
 need on-the-fly type conversion.

Notes:

obsolete:
Maybe numChannels should be a member variable instead of being passed to process. If the loaded
sample has a different number of channels than what we need to produce as output, we need sensible 
rules to deal with that situation, such as: If output is stereo and sample is mono: both outpus 
receive the same signal. If the situation is reversed, just the left channel goes into the mono 
output. Maybe something like: outChannel[i] = sampleChannel[i % numSampleChannels]. I think, the
most common situations are output: stereo, sample: mono or stereo. But it would be nice, if we 
could handle more general situations in a sensible way. For stereo-to-mono, it could be argued 
that a mix of both channels would be more sensible - but that rule doesn't seem to generalize 
well. But maybe we should have an (optional) exceptional rule for mono outputs to use a mixdown of
all channels.

-maybe make a nested namespace Sampler(Engine) (should later become part of rosic) - the nested
 classes are getting a bit unwieldy
-maybe rapt should be organized using nested namespaces - maybe look at the doxygen-generated
 API documentation, how this looks like

maybe rename to rsSampler, rsSoundFontPlayer, rsSamplerSFZ

If client code wants to modify regions and groups, it needs to do this by calling appropriate
functions on the rsSamplerEngine object with a pointer to the region or group to be modified
as parameter. For example:

  rsSamplerEngine sampler;
  // ...more stuff...
  using PST = rsSamplerEngine::PlaybackSetting::Type;
  sampler.setRegionSetting(region, PST::PitchKeyCenter, 69.f);

This is realized by having Region and Group define only private setters and letting rsSamplerEnigne
be a friend class, so it may access them. The goal is to prevent client code to modify regions and
groups behind the back of the sampler, because the sampler may need to take additional actions on 
such modifications. So, the sampler always acts as "man-in-the-middle" for any such changes to 
regions and groups. This might actually be a design pattern (-> figure out, if it's a known one).
This pattern should be used only for very closely coupled (ideally: nested) classes such is the 
case here.

Ideas:
-Could it make sense to define a level above the instrument - maybe an ensemble? Different 
 instruments in an ensemble could respond to different midi-channels. This would resemble the
 "multi-timbral" feature commonly seen in hardware romplers. But maybe that should be done in a
 class that contains a bunch (16) objects of type rsSamplerEngine. Maybe it should be called
 rsSamplerEnsemble or something. Maybe the samplePool object should than be shared among the
 embedded engines. Maybe it should be shared anyway to allow embedding the sampler as plugin in 
 a DAW and share the imported sample-content with it.
-Maybe at some point, we may want to provide more advanced envelope-generators such as the ones
 seen in Straightliner.
-The playback restriction based on last received controller values can be used to do a 
 waldorf-style wavetable synthesis: Create "WaveTable128" samples that contain 128 single cycles, 
 assign them to 128 regions and each such region is played only when the last received controller
 matches, like cycle 50 is played when the controller c that controls the cycle-number is in 
 49.5 <= c < 50.5. Maybe we can also use smaller wavetables (like 32) that use crossfading based
 on the controller...actually we should probably always crossfade between 2 cycles. At 49.5, we 
 would actually hear a 50/50 mix between cycle 49 and cycle 50
-To enable that feature, we should probably store the most recently received values of all 
 controllers
-Write a general "instrumentify" algorithm that takes as input a multisample of an instrument 
 (sampled at multiple keys and maybe also velocities) and produces a (extended) sfz instrument from
 it. The algorithm should include:
 -modeling (anti)formants via the sfz eq bands 
  -should be done on the instrument/masetr level - just one EQ for the sum
  -may need more than 3 bands
  -or, use a general pole/zero model (tweakable as sfz filter via its spectral centroid)
 -splitting of transient and body, use a sort of LA synthesis
 -maybe try to model the body with a delayline with a pole/zero model filter in its feedback path
  to model the freq-dependent decay (via its magnitudes) and the inharmonicity (via its phases)
 -perhaps split body further into harmonic/inharmonic/noisy parts
 -create an ambience sample that can be mixed in by
  -randomizing FFT phases (perhaps with a smeared magnitude spectrum)
  -convolving the sample with exponentially enveloped gaussian noise
  -these samples should get a filter evelope in the sfz, ideally, of a slope filter, i.e. a filter
   which doesn't have its cutoff changing over time, but its slope - this models faster decay of 
   high frequencies


Problem:
-SFZ actually allows for an unlimited number of regions and it seems that each region needs its
 own chain of DSP objects. Worse, when a region is retriggered and it is in "one-shot" mode, it is 
 supposed to be played twice (i.e. overlap with itself), etc. so it seems, we can't really 
 reasonably allocate "enough" DSP objects to be able to deal with any situation.
-Ideas: we could have a pool for any kind of supported DSP object (filter, eq, env-gen, lfo, etc.) 
 and whena noteOn is received, we build up the chain of DSP objects as needed by the region by 
 grabbing objects from the pools. Problem: the pools may run out of available objects.
-If we dynamically resize the pools, objects that are currently in use would be deallocated, so 
 that doesn't work. What we would need would be a sort of dynamically growing array that never
 deallocates - when it need to grow, it keeps the allocated memory allocated as is and allocates
 new memory somewher else - it wouldn't be contiguous anymore, but we would be safe from 
 deallocation. make classes rsGrowingArray and rsThreadedGrowingArray
-Whenever half of the objects are used up, we would allocate a new chunk of memory equal to the
 current size, so it would grow exponentially like dynamic arrays typically do.
-We should probably delegate the allocation of more memory to a worker thread to to its 
 nondeterministic runtime.
->figure out, how other sfz/sampler engines deal with this problem
-maybe implement a simple RegionPlayer class without any DSP (just pure sample playback) and 
 subclasses with various DSP objects
-We have 4 biquad objects in series. If the complete filter could be re-expressed as parallel
 structure, we could make more use of SIMD and could cut down the amount of state by factor 2
 by using single precision. Maybe we can express the filters as parallel connection of an allpass
 with a direct path (see DAFX). Maybe this can be done for envelopes and LFOs, too - there, it's
 even easier bacause they can be computed in parallel anyway. Make classes rsBiquadFloat64x4, etc.
-try to minimize the state of the RegionPlayer. maybe for the envelopes, it's enough to maintain
 the sample-time and compute the envelope without keeping state for each - but these are 
 optimization concerns to be addressed later

SFZ - Resources:
https://en.wikipedia.org/wiki/SFZ_(file_format)
https://github.com/sfz/tests/   test sfz files demonstrating various features
https://sfzformat.com/legacy/   opcode reference
https://sfzformat.com/headers/  reference for section headers in sfz files
http://www.drealm.info/sfz/plj-sfz.xhtml  description of the sfz format
https://www.kvraudio.com/forum/viewtopic.php?f=42&t=508861  kvr forum thread with documentation
https://sfzinstruments.github.io/  collection of sfz instruments
http://ariaengine.com/overview/sfz-format/
https://www.linuxsampler.org/sfz/    has convenient list of opcodes, also for sfz v2
http://doc.linuxsampler.org/sfz/

https://noisesculpture.com/cakewalk-synthesizers/
https://noisesculpture.com/cakewalk-synthesizers-downloads/


https://sfzformat.com/software/players/  players (also open source)
https://plugins4free.com/plugin/217/   sfz by rgcaudio

open source sfz players:
https://github.com/swesterfeld/liquidsfz/
https://sfz.tools/sfizz/downloads
https://github.com/altalogix/SFZero/
https://github.com/s-oram/Grace/

sfz compatibel samplers
https://github.com/christophhart/HISE/

deeper into the codebases:

https://github.com/swesterfeld/liquidsfz/tree/master/lib
https://github.com/swesterfeld/liquidsfz/blob/master/lib/synth.hh
This seems to do it the simple way: it has a fixed number of voices and if they are used up, no
more can be added - if i understand it correctly (see alloc_voice, line 230)



about float vs double:
https://randomascii.wordpress.com/2012/03/21/intermediate-floating-point-precision/

Ideas for new opcodes:
sample_dir=factory  (other options: user, here, E:/Samples/MySamples, ../../Samples/Piano, 
                     default: here)

maybe define a subregion header. Idea use the same sample with the mostly same settings but one or
a few settings differently for different keys...but no - this places too much burden on the 
playback engine - it would have to scan each region for subregions - no good idea!

*/