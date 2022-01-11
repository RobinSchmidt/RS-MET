namespace rosic {
namespace Sampler {

int rsSamplerEngine::instanceCounter = 0;

//-------------------------------------------------------------------------------------------------
// Lifetime:

rsSamplerEngine::rsSamplerEngine(int maxNumLayers)
{
  setMaxNumLayers(maxNumLayers);
  instanceCounter++;
  if(instanceCounter == 1)
    SfzCodeBook::createInstance();
    // This is a singleton object which is used for various translation tasks from various places.
    // We manage its lifetime here in this rather high-level class and lower level objects who live
    // only inside this one, such as rsSamplerData, will access the instance for various 
    // operations. We will also take care of the cleanup here.
}

rsSamplerEngine::~rsSamplerEngine()
{
  instanceCounter--;
  RAPT::rsAssert(instanceCounter >= 0);
  if(instanceCounter == 0)
    SfzCodeBook::deleteInstance();
}

//-------------------------------------------------------------------------------------------------
// Setup:

void rsSamplerEngine::setMaxNumLayers(int newMax)
{
  int L = newMax;
  playerPool.resize(L);
  idlePlayers.resize(L);
  activePlayers.reserve(L);
  for(int i = 0; i < L; i++)
  {
    playerPool[i].setDspResourcePool(&dspPool);
    playerPool[i].allocateMemory();
    idlePlayers[i] = &playerPool[i];
  }
}

void rsSamplerEngine::clearInstrument() 
{ 
  //RAPT::rsAssert(isMainThread()); 
  // todo: implement isMainThread and isAudioThread ...not yet sure how - maybe we should have a 
  // singleton in rosic storing two pointers to threads which are set up on construction maybe
  // rosic::rsThreadInfo::isMainThread() or something like that

  reset();                // Don't clear while still playing - that would be a crash!
  sfz.clearInstrument(); 
  samplePool.clear();
  setupRegionsForKey();   // clears regionsForKey array
}

int rsSamplerEngine::addSampleToPool(
  float** data, int numFrames, int numChannels, float sampleRate, const std::string& path)
{
  // todo: 
  // -check, if a sample with the same path already exists - if so, we have nothing to 
  //  do and may return early with an appropriate code
  //   if(isSampleInPool(..)) return rsReturnCode::nothingToDo;
  // -currently, we do such a check in addSamplesUsedIn - maybe it should be done here instead 
  //  and/or in loadSampleToPool

  AudioFileStreamPreloaded<float>* stream = new AudioFileStreamPreloaded<float>;
  bool allocOK = stream->setData(data, numFrames, numChannels, sampleRate, 2, path);
  if(allocOK == false)
    return rsReturnCode::memAllocFail;
  return samplePool.addSample(stream);
}

int rsSamplerEngine::loadSampleToPool(const std::string& path)
{
  if(isSampleInPool(path))
    return rsReturnCode::nothingToDo;

  int numFrames;
  int numChannels;
  int sampleRate;
  float** data = rosic::readFloatFromWaveFile(path.c_str(), numChannels, numFrames, sampleRate);
  if(data == nullptr)
    return rsReturnCode::fileLoadError;
    // This could mean that the file was not found or the memory allocation failed. ToDo: be more 
    // specific which of the two conditions happened

  int rc = addSampleToPool(data, numFrames, numChannels, (float) sampleRate, path);

  // Clean up data (todo: wrap into utility function - there is actually already one):
  //for(int c = 0; c < numChannels; c++)
  //  delete[] data[c];
  delete[] data[0];
  delete[] data;
  return rc;
}

int rsSamplerEngine::unUseSample(int i)
{
  if(i < 0 || i >= samplePool.getNumSamples()){
    RAPT::rsError("Invalid sample index");
    return rsReturnCode::invalidIndex; }
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
    return rsReturnCode::invalidIndex;    // gi was an invalid group index
  Region* r = getRegion(gi, ri);
  for(uchar k = loKey; k <= hiKey; k++)
    addRegionForKey(k, r);
  return ri;
}

int rsSamplerEngine::removeRegion(int gi, int ri)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }


  // Delete all our pointers to the region that will be deleted (maybe factor out into 
  // deleteAllPointersTo(r)):


  Region* r = getRegion(gi, ri);
  RAPT::rsAssert(r != nullptr);

  // We must stop all region players that make use of region r:
  for(int i = (int)activePlayers.size() - 1; i >= 0; i--) 
  {
    if(activePlayers[i]->getRegionToPlay() == r) 
    {
      rsReturnCode rc = stopRegionPlayer(i);
      RAPT::rsAssert(rc == rsReturnCode::success); 
    }
  }

  /*
  for(size_t i = 0; i < activePlayers.size(); i++)
  {
    if(activePlayers[i]->getRegionToPlay() == r)
    {
      activePlayers[i]->setRegionToPlay(nullptr, 0.0, groupSettingsOverride, regionSettingsOverride);
      // the idlePlayers are supposed to have a nullptr anyway
    }
  }
  */


  // We must remove all pointers to the region r from our regionsForKey array:
  for(int k = 0; k < numKeys; k++)
    regionsForKey[k].removeRegion(r);

  // Remove region from the rsSamplerData object
  bool success = sfz.removeRegion(gi, ri);
  if(success) 
    return rsReturnCode::success;
  else        
    return rsReturnCode::notFound;
}

int rsSamplerEngine::setRegionSample(int gi, int ri, int si)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  if(!isSampleIndexValid(si)) {
    RAPT::rsError("Invalid sample index");
    return rsReturnCode::invalidIndex; }

  const AudioFileStream<float>* s = samplePool.getSampleStream(si);
  sfz.setRegionCustomPointer(gi, ri, (void*) s);
  sfz.setRegionSample(gi, ri, s->getPath());
  return rsReturnCode::success;
}

int rsSamplerEngine::setRegionSetting(int gi, int ri, Opcode type, float value, int index)
{
  return sfz.setRegionSetting(gi, ri, type, value, index);
}

int rsSamplerEngine::setGroupSetting(int i, Opcode type, float value, int index)
{
  return sfz.setGroupSetting(i, type, value, index);
}

int rsSamplerEngine::setInstrumentSetting(Opcode type, float value, int index)
{
  return sfz.setInstrumentSetting(type, value, index);
}

int rsSamplerEngine::setupFromSFZ(const rsSamplerData& newSfz)
{
  removeSamplesNotUsedIn(newSfz);     // remove samples that are not needed anymore from memory
  int rc1 = addSamplesUsedIn(newSfz); // load samples that are needed but not yet loaded
  sfz = newSfz;                       // replace old sfz instrument definition member with new
  int rc2 = setupAudioStreams();      // connect regions in new sfz with appropriate stream objects
  setupRegionsForKey();               // updates regionsForKey array
  preAllocateDspMemory();             // allocates memory for the DSP objects, etc.
  if(rc1 >= 0 && rc2 == rsReturnCode::success)
    return rsReturnCode::success;
  else
    return rsReturnCode::fileLoadError;  
    // ToDo: be more specific about the error condition

  // Maybe have state variables numSamplesAdded, numSamplesRemoved, numSamplesNotFound that can be 
  // inquired from client
  // code. this is useful mainly for unit tests but maybe it makes sense to display such info on 
  // the GUI, too, so the user has some feedback about what patches have a lot of samples in common


  // This function needs to be mutexed with anything that accesses the stream-pointers in the sfz
  // objects
}

bool rsSamplerEngine::setSfzRootDir(const char* path) 
{ 
  sfzDir = path;
  return true;

  // preliminary - todo:

  //if(!rsFile::isValidPath(path))
  //  return false;
  //else {
  //  sfzDir = path;
  //  return true; }

  // ToDo: 
  // -Check, if the last character is the seperator character, i.e. "/" or "\". This is currently 
  //  assumed to be the case. If it's not, either append it here to fix the situation or raise an
  //  error...
}

bool rsSamplerEngine::saveToSFZ(const char* path, bool pathIsAbsolute) const
{
  std::string absPath = getAbsolutePath(path, pathIsAbsolute);
  return sfz.saveToSFZ(path);
}

int rsSamplerEngine::loadFromSFZ(const char* path, bool pathIsAbsolute)
{
  std::string absPath = getAbsolutePath(path, pathIsAbsolute);
  rsSamplerData newSfz;
  rsReturnCode rc = newSfz.loadFromSFZ(absPath.c_str());
  if(rc != rsReturnCode::success) {
    clearInstrument();
    return rc; }
  return setupFromSFZ(newSfz);

  // Maybe this function should load the file and then call setFromSFZ instead of using 
  // newSfz.loadFromSFZ etc. There's a lot of code duplication between this function and 
  // setFromSFZ. Doing it the other way may require use to duplicate some code from 
  // newSfz.loadFromSFZ but maybe that's the lesser evil.
}

int rsSamplerEngine::setFromSFZ(const std::string& sfzFileContents)
{
  rsSamplerData newSfz;
  rsReturnCode rc = newSfz.setFromSFZ(sfzFileContents);
  if(rc != rsReturnCode::success) {
    clearInstrument();
    return rc; }
  return setupFromSFZ(newSfz);
}

//-------------------------------------------------------------------------------------------------
// Inquiry:

rsSamplerEngine::Region* rsSamplerEngine::getRegion(int groupIndex, int regionIndex)
{
  int gi = groupIndex, ri = regionIndex;
  if(gi < 0 || gi >= (int)sfz.instrument.groups.size()) {
    RAPT::rsError("Invalid group index");
    return nullptr; }
  return sfz.instrument.groups[gi]->getRegion(ri);
}

int rsSamplerEngine::getNumRegionsUsing(int i) const
{
  if(i < 0 || i >= samplePool.getNumSamples()){
    RAPT::rsError("Invalid sample index");
    return rsReturnCode::invalidIndex; }
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

std::string rsSamplerEngine::getAbsolutePath(const char* path, bool pathIsAbsolute) const
{
  std::string absPath;
  if(pathIsAbsolute)
    absPath = std::string(path);
  else
    absPath = sfzDir + std::string(path); // or do we need to insert a separator "/" or "\"?
  return absPath;
}

//-------------------------------------------------------------------------------------------------
// Processing:

void rsSamplerEngine::processFrame(double* left, double* right)
{
  rsFloat64x2 out = 0.0;

  for(int i = 0; i < (int)activePlayers.size(); i++) {
    out += activePlayers[i]->getFrame();
    if(activePlayers[i]->hasFinished()) {
      stopRegionPlayer(i);
      i--;  }}
    // ToDo: Test, if it's more efficient to loop through the activePlayers array backwards. Then, 
    // the i-- in the loop body could be removed, but that's not the main point. The main point is 
    // that the deactivation/removal would need less data copying.


  // probably obsolete - that stuff is handled in the subclass:
  //if(groupSettingsOnTop)
  //{
  //  for(int i = 0; i < (int)activeGroupPlayers.size(); i++)
  //    out += activeGroupPlayers[i]->getFrame();
  //}
  //else
  //{
  //  for(int i = 0; i < (int)activePlayers.size(); i++) {
  //    out += activePlayers[i]->getFrame();
  //    if(activePlayers[i]->hasFinished()) {
  //      deactivateRegionPlayer(i);
  //      i--;  }}
  //  // ToDo: Test, if it's more efficient to loop through the activePlayers array backwards. Then, 
  //  // the i-- in the loop body could be removed, but that's not the main point. The main point is 
  //  // that the deactivation/removal would need less data copying.
  //}


  *left  = out[0];
  *right = out[1];

  // ToDo:
  // -Try to get rid of the conditional by always using the groupPlayers array. If we are in
  //  non-on-top mode, just use always just a single active GroupPlayer (i.e. the 
  //  activeGroupPlayers just has a length of 1). But then, the (then empty) GroupPlayer's 
  //  dspChain wil always be applied...but maybe with block-processing, the cost will be negligible
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

rsSamplerEngine::PlayStatusChange rsSamplerEngine::handleMusicalEvent(
  const rsMusicalEvent<float>& ev)
{
  using Type = rsMusicalEvent<float>::Type;
  Type  type = ev.getType();
  float val1 = ev.getValue1();
  float val2 = ev.getValue2();
  switch(type)
  {
  case Type::noteOn:  return handleNoteOn( (uchar)val1, (uchar)val2); break;
  case Type::noteOff: return handleNoteOff((uchar)val1, (uchar)val2); break;
  // Maybe we should pass the floats directly to noteOn/Off. This would allow for microtuning. 
  // However, it complicates identifying to which notes a noteOff should apply. Or rather, it's 
  // not really complicated but requires a design decision: should we require an exact match or 
  // allow the range key +-0.5? For the moment, this here is good enough, though.

  default: return PlayStatusChange();
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
  if(key < r->getLoKey() || key > r->getHiKey()) return false;
  if(vel < r->getLoVel() || vel > r->getHiVel()) return false;
  if(r->getCustomPointer() == nullptr) return false;


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
  RAPT::rsError("Region not found");
  // A region should always be found in one (and only one group). If we don't find it, it means
  // the caller has passed a pointer to a region object that is not part of this instrument. If 
  // this happens, it indicates a bug at the call site.
}

RegionPlayer* rsSamplerEngine::getRegionPlayerFor(
  const Region* r, uchar key, uchar vel)
{
  RAPT::rsAssert(r->getCustomPointer() != nullptr); // No stream connected

  // The behavior for this in situations where the region r is already being played by some 
  // voice/player should depend on the playback-mode of the region. If it's in one-shot mode, a new
  // player should be handed and the old one should just continue to play along with the new. 
  // Otherwise, the old player should be re-used. But this may cause clicks. Maybe a new player 
  // should be handed and the old one should be flagged for a quick fade-out (a few milliseconds). 
  // Maybe that feature should be optional, controlled by a retriggerFadeOutTime parameter that is 
  // zero by default (indicating hard, clicking retriggers). Or maybe this fade-out could also be 
  // set per region instead of globally? Check, if sfz may actually have an opcode for that.
  if(idlePlayers.empty())
    return nullptr;  // Maybe we should implement more elaborate voice stealing?
  RegionPlayer* rp = RAPT::rsGetAndRemoveLast(idlePlayers);
  rp->setKey(key);  // why is it not enough to do it inside the "if"
  rsReturnCode rc = rp->setRegionToPlay(r, getSampleStreamFor(r), sampleRate, busMode);
  if(rc == rsReturnCode::success)
  {
    //rp->setKey(key);
    activePlayers.push_back(rp);
    return rp;
  }
  else
  {
    idlePlayers.push_back(rp);
    return nullptr;
  }
  //activePlayers.push_back(rp);
  //return rp;

  // ToDo: Let setRegionToPlay return a bool to indicate, if it's possible to play the region. It
  // may fail due to lack of ressources, such as DSP processors and/or modulators. If it returns 
  // false, roll it back, i.e. move the player back to the idlePlayers array (maybe reset any 
  // options that have been set such as setKey) and return a nullptr
}

bool rsSamplerEngine::isSampleUsedIn(
  const AudioFileStream<float>* stream, const rsSamplerData& sfz)
{
  const std::string& streamPath = stream->getPath();
  for(int gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group* g = sfz.getGroup(gi);
    for(int ri = 0; ri < g->getNumRegions(); ri++) {
      Region* r = g->getRegion(ri);
      const std::string& regionPath = r->getSamplePath();
      if(regionPath == streamPath)
        return true; }}
  return false;
}

rsReturnCode rsSamplerEngine::stopRegionPlayer(int i)
{
  if(i < 0 || i >= activePlayers.size()) {
    RAPT::rsError("Invalid player index");
    return rsReturnCode::invalidIndex; }
  RegionPlayer* p = activePlayers[i];
  RAPT::rsRemove(activePlayers, i);
  p->releaseResources();
  idlePlayers.push_back(p);
  return rsReturnCode::success;
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

rsSamplerEngine::PlayStatusChange rsSamplerEngine::handleNoteOn(uchar key, uchar vel)
{
  if(vel == 0) { return handleNoteOff(key, vel); }

  PlayStatusChange psc;
  for(int i = 0; i < regionsForKey[key].getNumRegions(); i++) 
  {
    const Region* r  = regionsForKey[key].getRegion(i);
    if(!shouldRegionPlay(r, key, vel))
      continue;
    RegionPlayer* rp = getRegionPlayerFor(r, key, vel);
    if(rp == nullptr) {
      // Roll back the addition of all players so far to the activePlayers and move them back into
      // the idlePlayers again. We don't really want notes to play with an incomplete set of 
      // samples. It's all or nothing - either all regions for the given key get triggered or none
      // of them:
      for(int j = 0; i < psc.numLayersStarted; j++) {
        rp = RAPT::rsGetAndRemoveLast(activePlayers);
        idlePlayers.push_back(rp); }
      psc.numLayersStarted = 0;
    }
    else
    {
      //rp->setKey(key);
      psc.numLayersStarted++;
    }
  }
  return psc;


  //return rsReturnCode::success;
  // Another possibility for the return value would have been to return the number of layers that
  // have been triggered, but we don't do that because then it would be not quite clear what we
  // should return from noteOff to make the functions somewhat consistent. In noteOff, we could
  // either return the number of released regions or the number of triggered release samples. Both
  // would make just as much sense. So, for unambiguous consistency, we let both just return a 
  // success or failure report.

  // ToDo: 
  // -maybe introduce a sort of PlayStatusChange struct, containing fields for: numLayersStarted,
  //  numLayersStopped
}

rsSamplerEngine::PlayStatusChange rsSamplerEngine::handleNoteOff(uchar key, uchar vel)
{
  // Note-off events may also trigger the playback of regions (note-off samples are a thing)...
  //int numRegions = 0;
  PlayStatusChange psc;
  for(int i = int(activePlayers.size()) - 1; i >= 0; i--)
  {
    if(activePlayers[i]->getKey() == key)
    {
      stopRegionPlayer(size_t(i));
      psc.numLayersStopped++;
      // ToDo: refine this later: we may not want to immediately stop the player but rather 
      // trigger the release phase and mark for quick fade-out
    }
  }
  return psc;


  //return rsReturnCode::success; // preliminary

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
      numSamplesRemoved++; }}
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
  for(int gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group* g = sfz.getGroup(gi);
    for(int ri = 0; ri < g->getNumRegions(); ri++) {
      Region* r = g->getRegion(ri); 
      const std::string& relativePath = r->getSamplePath();
      std::string path = getAbsolutePath(relativePath.c_str());  // preliminary
      if(isValidSamplePath(path) && !isSampleInPool(path)) {
        int rc = loadSampleToPool(path);
        if(rc >= 0)
          numSamplesLoaded++;
        else
          numSamplesFailed++; }}}
  if(numSamplesFailed == 0)
    return numSamplesLoaded;
  else
    return rsReturnCode::fileLoadError;

  // ToDo:
  // -the preliminary path = getAbsolutePath(relativePath.c_str()); should probably not use the
  //  getAbsolutePath function which is actually meant for the sfz files. Instead, we should 
  //  perhaps interpret the r->getSamplePath() as being relative to the .sfz file. At the moment,
  //  we interpret it as relative to the sfz root folder. If the sfz file actually resides in that
  //  sfz root folder, all is well. But in general, we should not assume this to be the case.
  // -do this chaneg also in setupAudioStreams below
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
  for(int gi = 0; gi < sfz.getNumGroups(); gi++) 
  {
    const Group* g = sfz.getGroup(gi);
    for(int ri = 0; ri < g->getNumRegions(); ri++) 
    {
      Region* r = g->getRegion(ri);
      const std::string& relativePath = r->getSamplePath();
      std::string path = getAbsolutePath(relativePath.c_str());  // preliminary
      allOK &= setupStream(r, path);
    }
  }
  // ToDo: 
  // -See comment in addSamplesUsedIn - the same applies here as well
  // -Maybe, if getSamplePath returns the empty string, meaning that a region does not define
  //  the sample opcode, assign the stream from the outlying group. If that also doesn't define a
  //  sample, use the stream from the instrument.

  if(allOK)
    return rsReturnCode::success;
  else
    return rsReturnCode::notFound;  // stream was not found for one or more samples
}

void rsSamplerEngine::setupRegionsForKey()
{
  for(uchar k = 0; k < numKeys; k++)
    regionsForKey[k].clear();
  for(int gi = 0; gi < sfz.getNumGroups(); gi++) {
    const Group* g = sfz.getGroup(gi);
    for(int ri = 0; ri < g->getNumRegions(); ri++) {
      Region* r = g->getRegion(ri);
      for(uchar k = 0; k < numKeys; k++)
        addRegionForKey(k, r); }}
}

void rsSamplerEngine::preAllocateDspMemory()
{
  //RAPT::rsWarning("rsSamplerEngine::preAllocateDspMemory still under construction");

  // ToDo:
  // -allocate enough RegionPlayers in our playerPool - how much is enough should depend on how 
  //  many regions could potentially paly simultaneously. To figure out that number, we need to
  //  count for each key the number of regions it may play...but it gets more complicated when
  //  release-samples and overlapping one-shot samples are considered, so maybe that's 
  //  impractical and we should just use some reasonably large-ish number...like 512 maybe?
  // -allocate "enough" DSP objects in our dspPool - here, the number will be even more complicated
  //  to calculate...maybe allocate for each region player 2 filters and 5 equalizers...maybe let
  //  the user control the amount of "allocation-headroom" and provide feedback about used-up 
  //  resources. trying to make that bulletproof seems impossible: you can't have "potentially 
  //  infinite" resources while at the same time guaranteeing no drop-outs to happen ever. We 
  //  either have to re-allocate if needed or put a fixed upper limit on the resources. But maybe
  //  we can let the user choose which of the two it should be
  // -allocate enough memory for the RegionPlayer's arrays: 
  //    std::vector<Modulator*> modulators;
  //    std::vector<ModulationConnection*> modMatrix;
  //    SignalProcessorChain dspChain;
  // -to do this, we must make a judgement how much is "enough" - if it later turns out to be not 
  //  enough, we may either not play back the voice or do the forbidden thing: allocate more memory
  //  in the audio thread or find some other way to handle the situation. Maybe some strategy that
  //  detects that we may soon run out of pre-allocated memory and firing up another thread to 
  //  allocate us some more could work? But that wouldn't be totally foolproof either. We'll see...
  //  ...see what hise, sfizz and linuxsampler do
}

//=================================================================================================

rsSamplerEngine2::rsSamplerEngine2(int maxNumLayers)
{
  setMaxNumLayers(maxNumLayers);
}

void rsSamplerEngine2::setMaxNumLayers(int newMax)
{
  rsSamplerEngine::setMaxNumLayers(newMax);
  int L = newMax;
  groupPlayerPool.resize(L);
  idleGroupPlayers.resize(L);
  activeGroupPlayers.reserve(L);
  for(int i = 0; i < L; i++) {
    idleGroupPlayers[i] = &groupPlayerPool[i];
    //idleGroupPlayers[i]->engine = this; // backlink to owner (us) currently not needed there
  }
}

rsSamplerEngine::PlayStatusChange rsSamplerEngine2::handleNoteOn(uchar key, uchar vel)
{
  PlayStatusChange psc = rsSamplerEngine::handleNoteOn(key, vel);
  if(!canFallBackToBaseclass()) 
  { 
    // Don't update the group players, if they are not used anyway
    updateGroupPlayers(psc);

  }
  return psc;
}

rsSamplerEngine::PlayStatusChange rsSamplerEngine2::handleNoteOff(uchar key, uchar vel)
{
  PlayStatusChange psc = rsSamplerEngine::handleNoteOff(key, vel);
  if(!canFallBackToBaseclass())
  {
    updateGroupPlayers(psc);
  }
  return psc;
}

rsReturnCode rsSamplerEngine2::stopRegionPlayer(int activeIndex)
{
  // ToDo: add code that verifies that rp is in exactly one of the active groupPlayers - maybe have
  // rsAssert(getNumContainingActiveGroupPlayers(rp) == 1) or something like that. that maybe 
  // useful to call in other places, too

  // Figure out, to which GroupPlayer the RegionPlayer with given activeIndex belongs and remove it
  // from the GroupPlayer. If this removal causes the GroupPlayer to become empty, also remove the 
  // GroupPlayer itself from the active ones:
  RegionPlayer* rp = activePlayers[activeIndex];
  for(int i = 0; i < (int)activeGroupPlayers.size(); i++) {
    if(activeGroupPlayers[i]->contains(rp)) {
      GroupPlayer* gp = activeGroupPlayers[i];
      gp->removeRegionPlayer(rp);
      if(gp->hasNoRegionPlayers()) {
        stopGroupPlayer(i);          // todo: take return code into account
        break; }}}                   // there can be only one and we have found it - we are done!
  return rsSamplerEngine::stopRegionPlayer(activeIndex);
}

void rsSamplerEngine2::processFrame(double* left, double* right)
{
  // Fall back to more efficient (because simpler) baseclass code, if possible:
  if(canFallBackToBaseclass()) { 
    rsSamplerEngine::processFrame(left, right); return; } 

  // Accumulate output from the group players:
  rsFloat64x2 out = 0.0;
  for(int i = 0; i < (int)activeGroupPlayers.size(); i++)
    out += activeGroupPlayers[i]->getFrame();

  // Stop region players that have finished playing:
  for(int i = 0; i < (int)activePlayers.size(); i++) {
    if(activePlayers[i]->hasFinished()) {
      stopRegionPlayer(i);
      i--; }}

  // Assign outputs:
  *left  = out[0];
  *right = out[1];
}

int rsSamplerEngine2::stopAllPlayers()
{
  int numRegionPlayersStopped = rsSamplerEngine::stopAllPlayers();
  int numGroupPlayers = (int) activeGroupPlayers.size();
  for(int i = numGroupPlayers-1; i >= 0; i--) {
    activeGroupPlayers[i]->releaseResources();
    idleGroupPlayers.push_back(activeGroupPlayers[i]); }
  activeGroupPlayers.clear();
  return numRegionPlayersStopped;
}
// ToDo: 
// -in addition to move the players back into their idle pool, we also need to move the dsp
//  objects back


void rsSamplerEngine2::updateGroupPlayers(PlayStatusChange psc)
{
  // Figure out, how many regions were triggered and add pointers to the freshly triggered 
  // RegionPlayers to an appropriate GroupPlayer - either to one that is already playing or grab a 
  // new one from the idle ones and add it to the active ones:
  int numLayersNow    = getNumActiveLayers();
  int numLayersBefore = numLayersNow - psc.numLayersStarted;
  for(int i = numLayersBefore; i < numLayersNow; i++) {
    RegionPlayer* rp = activePlayers[i];
    const rsSamplerData::Group* grp = rp->getRegionToPlay()->getGroup();
    int gpi = getActiveGroupPlayerIndexFor(grp);
    if(gpi != -1)
      activeGroupPlayers[gpi]->addRegionPlayer(rp);
    else
      startGroupPlayerFor(rp); }

  // If nothing was playing before, we have to start the instrumPlayer, too:
  if(numLayersBefore == 0)
    startInstrumPlayer();
}

int rsSamplerEngine2::getActiveGroupPlayerIndexFor(const rsSamplerData::Group* group)
{
  for(size_t i = 0; i < activeGroupPlayers.size(); i++)
    if(activeGroupPlayers[i]->group == group)
      return (int) i;
  return -1;
}

void rsSamplerEngine2::startGroupPlayerFor(RegionPlayer* rp)
{
  GroupPlayer* gp = RAPT::rsGetAndRemoveLast(idleGroupPlayers);
  const rsSamplerData::Group* grp = rp->getRegionToPlay()->getGroup();
  gp->addRegionPlayer(rp);
  gp->group = grp;

  gp->assembleDspChain(busMode); // ToDo: use boolean return value: if false (unlikely), roll back the RegionPlayer rp
  //gp->setupDspSettings(); // to do

  activeGroupPlayers.push_back(gp);
}

int rsSamplerEngine2::stopGroupPlayer(int i)
{
  if(i < 0 || i >= activeGroupPlayers.size()) {
    RAPT::rsError("Invalid player index");
    return rsReturnCode::invalidIndex; }
  GroupPlayer* p = activeGroupPlayers[i];
  p->releaseResources();
  RAPT::rsRemove(activeGroupPlayers, i);
  idleGroupPlayers.push_back(p);

  if(activeGroupPlayers.size() == 0)  // If nothing is playing anymore...
    stopInstrumPlayer();              // ...stop the whole instrumPlayer, too
  // hmmm ...maybe we should allow for a tail ringout?

  return rsReturnCode::success;
}

void rsSamplerEngine2::startInstrumPlayer()
{

  return;
}

void rsSamplerEngine2::stopInstrumPlayer()
{


  return;
}



}} // namespaces

//=================================================================================================

/*

Bugs:
-when loading a new instrument while a region is playing, it crashes
-group volume seems to work now, but instrument volume seems to be ignored...or it even get muted 
 when having a top-level volume. might be an issue with the parser -> check sfz spec


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
-Optionally resample the samples on loading to the output sample rate using a very high quality 
 (much better than realtime) resampling algorithm (say Sinc_512)

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


maybe rename to rsSampler, rsSamplerSFZ

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
-For the filter, make a special class rsSamplerFilter or rsMultiModeFilter which may interpret its 
 states and coeffs differently, depending on the selected mode (can implement ladder, svf (2 in 
 series or parallel), biquad, etc.).
-Maybe implement the accumulating mode already in the baseclass but there, do indeed use 3 filters
 per region player. The subclass is responsible for the optimization of using only 1 filter per 
 group and only one for the whole instrument.
-Maybe for the signals, use rsFloat32x4. fits better with usage of float for the parameters. Maybe
 we can use the 2 extra channels to our advantage as well. Maybe they can store the Hilbert trafo
 of the signals in the actual channels. Or maybe it can carry along some metadata that can be 
 used in subsequent DSP algos. Maybe instantaneous amp/phase? Maybe mid/side signals? We'll find a 
 way of using them. ...but it should be some data that can be easily generated in realtime because
 we don't want to store it in the RAM. Or maybe we could process 2 sample frames at a time? But if
 this turns out to be doable, there's no reason to stop at 2. If we have enough simd width, we 
 could also process 8 samples at a time. This could work for waveshapers and FIR filters...but IIR
 filters? maybe not....
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

-maybe have a multi-output version (e.g. 16 stereo outs) and allow groups to be routed to different
 output busses. 16 seems nice because of the possible correspondence to the 16 available midi 
 channels. an/or maybe introduce another "ensemble" level above the instrument and allows whole 
 instruments to go to different busses - should go into rsSamplerEngine2

Questions:
-How are the group-DSP processes supposed to be modulated? For this to make any sense, we 
 would indeed need per-group modulators. But that does not really fit with the idea that 
 RegionPlayers need to duplicate the modulators when there's group modulation on top. Take an 
 envelope for filter cutoff..the per-region envelope and per-group envelope would get added to
 obtain the *region's* cutoff, if modulators accuumulate. ..but an additonal per-group filter would
 just remain static?


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





about float vs double:
https://randomascii.wordpress.com/2012/03/21/intermediate-floating-point-precision/

Ideas for new opcodes:
sample_dir=factory  (other options: user, here, E:/Samples/MySamples, ../../Samples/Piano, 
                     default: here)

maybe define a subregion header. Idea use the same sample with the mostly same settings but one or
a few settings differently for different keys...but no - this places too much burden on the 
playback engine - it would have to scan each region for subregions - no good idea!

*/