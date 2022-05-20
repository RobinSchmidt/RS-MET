namespace rosic {
namespace Sampler {

//-------------------------------------------------------------------------------------------------
// The internal classes

float PlaybackSetting::getDefaultValue(Opcode type, int index)
{
  SfzCodeBook* t = SfzCodeBook::getInstance();
  return t->opcodeDefaultValue(type, index);
}
OpcodeType PlaybackSetting::getTargetOpcodeType(Opcode type)
{
  SfzCodeBook* t = SfzCodeBook::getInstance();
  return t->getOpcodeType(type);
}
// Maybe get rid of these two. The caller should use the translator himself. Or maybe it's more 
// convenient to keep them as convenience functions? We'll see....

//-------------------------------------------------------------------------------------------------

void SfzInstrument::HierarchyLevel::ensureDspsPresent(Opcode opcodeType, int howMany)
{ 
  using namespace RAPT;
  using SPT = OpcodeType;
  SPT dspType = PlaybackSetting::getTargetOpcodeType(opcodeType);
  rsAssert(dspType != SPT::Unknown, "Opcode refers to unknown kind of DSP object");

  // use: SfzCodeBook::isEffectSetting:
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

void SfzInstrument::HierarchyLevel::setFilterEnvDepth(float depthInCents)
{
  using OT = OpcodeType;
  int numFilters = (int)RAPT::rsCount(dspTypes, OT::Filter);
  for(int i = 0; i < numFilters; i++)
    setModulation(OT::FilterEnv, 1, Opcode::cutoffN, i+1, depthInCents, ModMode::cents);
}

void SfzInstrument::HierarchyLevel::setFilterLfoDepth(float depthInCents)
{
  using OT = OpcodeType;
  int numFilters = (int)RAPT::rsCount(dspTypes, OT::Filter);
  for(int i = 0; i < numFilters; i++)
    setModulation(OT::FilterLfo, 1, Opcode::cutoffN, i+1, depthInCents, ModMode::cents);
}

bool isEffect(OpcodeType ot)
{
  using OT = OpcodeType;
  return ot > OT::_TagEffectsStart && ot < OT::_TagEffectsEnd;
}
// maybe move as static member function into SfzCodeBook..there actually already is such a function
// -> use that and remove this!

void SfzInstrument::HierarchyLevel::setAmpEnvDepth(float depthInPercent)
{
  // If the last effect DSP in the chain is not an Amplifier with amplitudeN==0, we need to insert 
  // another Amplifier at the end.

  using OT = OpcodeType;
  using OC = Opcode;

  // Figure out the amplitudeN opcode/parameter of the last Amplifier in the chain, assume 1.f if 
  // there are no Amplifiers:
  int numAmps = (int)RAPT::rsCount(dspTypes, OT::Amplifier);
  float lastAmpParam = 1.f;
  if(numAmps > 0)
    lastAmpParam = getSettingValue(OC::amplitudeN, numAmps);

  // Append another amplifier with amplitudeN=0 if necessary:
  bool newAmpNeeded = !isLastEffectAmplifier() || lastAmpParam != 0.f;
  if(newAmpNeeded) {
    numAmps++; setSetting(PlaybackSetting(OC::amplitudeN, 0.f, numAmps)); }

  // OK, now we have ensured that the last effect in the chain is indeed an Amplifier with a 
  // setting of amplitudeN=0 and its index all other Amplifiers (i.e. the N in the amplitudeN 
  // opcode - *not* the index in the dspTypes array). This is the amplifier to which we route the
  // amp-env:
  setModulation(OT::AmpEnv, 1, OC::amplitudeN, numAmps, depthInPercent, ModMode::absolute);
}
// needs unit tests under various circumstances. The logic is quite complicated...

void SfzInstrument::HierarchyLevel::setAmpLfoDepth(float depth)
{
  // The logic to figure out whether or not we need to append an Amplifier is similar to 
  // setAmpEnvDepth except that we don't need the lastAmpParam == 0.f constraint.

  int numAmps = (int)RAPT::rsCount(dspTypes, OpcodeType::Amplifier);
  if(!isLastEffectAmplifier()) {
    numAmps++; setSetting(PlaybackSetting(Opcode::volumeN, 0.f, numAmps)); }
  setModulation(OpcodeType::AmpLfo, 1, Opcode::volumeN, numAmps, depth, ModMode::absolute);
}

void SfzInstrument::HierarchyLevel::setSetting(const PlaybackSetting& s)
{
  using OT = OpcodeType;
  using OC = Opcode;
  OC op = s.getOpcode();

  // Handle the lo/hi key/vel opcodes as special cases:
  if(op == OC::LoKey) { loKey = (uchar)s.getValue(); return; }
  if(op == OC::HiKey) { hiKey = (uchar)s.getValue(); return; }
  if(op == OC::LoVel) { loVel = (uchar)s.getValue(); return; }
  if(op == OC::HiVel) { hiVel = (uchar)s.getValue(); return; }
  // ToDo: maybe we should assert that the value is an integer in the range 0..127
  // we should also handle the "key" opcode which specifies lokey, hikey, 
  // pitch_keycenter simultaneously?

  // The "key" opcode specifies lokey, hikey and pitch_keycenter at the same time:
  if(op == OC::Key)
  {
    loKey = hiKey = (uchar)s.getValue();
    setSetting(PlaybackSetting(Opcode::PitchKeyCenter, s.getValue(), -1));
    // Shouldn't we return here? If not, document why not.
  }

  // Some modulation routing opcodes need to be also handled as special cases:
  if(op == OC::fileg_depth)  { setFilterEnvDepth(s.getValue()); return; }
  if(op == OC::fillfo_depth) { setFilterLfoDepth(s.getValue()); return; }
  if(op == OC::ampeg_depth)  { setAmpEnvDepth(   s.getValue()); return; }
  if(op == OC::amplfo_depth) { setAmpLfoDepth(   s.getValue()); return; }
  // ToDo: pitcheg_depth, pitchlfo_depth

  // All other settings are handled by either overwriting the last setting of that type in our 
  // array, if present or by appending the setting, if not present:
  int idx = s.getIndex();
  int foundAt = findSetting(op, idx);
  if(foundAt != -1)
    settings[foundAt] = s;
  else
  {
    settings.push_back(s);
    ensureDspsPresent(op, RAPT::rsMax(idx, 1));
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

  // We should probably also remove all modulation routings that have this setting as target?
  // ...or maybe not?

  // ToDo:
  // -Try to clean this up by getting rid of the index == -1 branch: Maybe don't allow the index 
  //  to be -1. If no index applies to a given Opcode, use index 1 as default. But in some places
  //  we rely on the encoding that -1 means "not applicable"...this would then need to be changed.
  // -Maybe implement a general library function removeAllMatches taking a vector and a predicate
  //  and returning the number of removed items. The predicate can be passed as lambda. We can then
  //  just do (assuming 1st branch was removed before): 
  //    int numRemoved = RAPT::rsRemoveAllMatches(settings, 
  //      [&](PlaybackSetting& s){ return s.getOpcode() == type && s.getIndex() == index; } );
  //    updateDspsArray();
  //    return numRemoved > 0;
  //  see also implementation of removeModulation.
}

void SfzInstrument::HierarchyLevel::setModulation(OpcodeType modSrcType, int modSrcIndex,
  Opcode modTarget, int modTargetIndex, float modDepth, ModMode modMode)
{
  // ToDo: verify that modTarget is a modulatable opcode, i.e. an opcode for setting a continuous 
  // parameter - something like cutoffN but not something like loop_mode - but maybe that should be
  // checked when assembling the modMatrix...nope - let's do it here - that's more efficient 
  // because it's called less often

  if(!SfzCodeBook::isModSourceSetting(modSrcType)) {
    RAPT::rsError("modSrcType must be some sort of modulation source");
    return; }
  ModulationRouting newRouting(
    modSrcType, modSrcIndex, modTarget, modTargetIndex, modDepth, modMode);

  // new:
  setModulation(newRouting);

  /*
  // old:
  for(size_t i = 0; i < modRoutings.size(); i++) {
    if( modRoutings[i].hasMatchingEndpoints(newRouting) ) {
      modRoutings[i].setDepth(modDepth);  // Overwrite existing routing
      modRoutings[i].setMode(modMode);
      return;  }}                         // ...and return early
  modRoutings.push_back(newRouting);      // or append a new routing
  */
}
// Maybe return the index where the setting was modified or the last index if it was pushed.
// Maybe we should have a sanity check that source and target exist...but maybe it's allowable to 
// set up the modulation connections before these exist. But in this case, we should perhaps make
// a sanity check for all connections after the region was fully set up. Maybe a member function
// checkModRoutingSanity ..or hasDanglingRoutings or something

void SfzInstrument::HierarchyLevel::setModulation(const ModulationRouting& r)
{
  for(size_t i = 0; i < modRoutings.size(); i++) {
    if( modRoutings[i].hasMatchingEndpoints(r) ) {
      modRoutings[i].setDepth(r.getDepth());  // Overwrite existing routing
      modRoutings[i].setMode(r.getMode());
      return;  }}                             // ...and return early
  modRoutings.push_back(r);                   // or append a new routing
}

template<class T, class P>
size_t rsRemoveIf(std::vector<T>& vec, P pred)
{
  if(vec.empty())
    return 0;
  size_t numRemoved = 0;
  size_t i = vec.size()-1;
  while(true) {
    if(pred(vec[i])) {
      RAPT::rsRemove(vec, i);
      numRemoved++; }
    if(i == 0)
      break;
    --i; }
  return numRemoved;
}
// move to RAPT, implement unit tests
// This is actually similar std::remove_if but the std function returns an iterator instead of a 
// number of removed items...which is not so useful, in my opinion. Maybe this could be implemented
// more efficiently using an int for i - we could just have a while(i >= 0) loop...but maybe we can
// also do it with size_t by using while(i != -1) using the fact that -1 corresponds to 
// max(size_t)...but that's ugly

bool SfzInstrument::HierarchyLevel::removeModulation(OpcodeType modSrcType, int modSrcIndex, 
  Opcode modTargetOpcode, int modTargetIndex)
{
  size_t numRemoved = rsRemoveIf(modRoutings, [&](ModulationRouting& s){ 
    return s.getSourceType() == modSrcType && s.getSourceIndex() == modSrcIndex && 
      s.getTargetOpcode() == modTargetOpcode && s.getTargetIndex() == modTargetIndex; });
  return numRemoved > 0;
}
// I think, numRemoved should always be either 0 or 1. Any other number would mean that we had two
// or more connections between the same pair of pins, which should not happen. When setModulation 
// is called for a pair of pins between which there already is a connection, the existing 
// connection should be updated rather than adding a second connection. In principle, the 
// architecture of the modulation system should allow mutliple connections between the same pair of
// pins but that doesn't seem to be a useful feature, so we don't do that.

void SfzInstrument::HierarchyLevel::connectAmpEnv()
{
  // New:
  int numAmpEnvs = (int)RAPT::rsCount(dspTypes, OpcodeType::AmpEnv);
  if(numAmpEnvs > 0) {
    RAPT::rsAssert(numAmpEnvs == 1, "Number of AmpEns should be zero or one");
    setAmpEnvDepth(100.f); }

  int dummy = 0;
}

void SfzInstrument::HierarchyLevel::copyDataFrom(const HierarchyLevel* lvl)
{
  samplePath  = lvl->samplePath;
  settings    = lvl->settings;
  modRoutings = lvl->modRoutings;
  dspTypes    = lvl->dspTypes;

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

int SfzInstrument::HierarchyLevel::getNumEffects() const
{
  int count = 0;
  for(size_t i = 0; i < dspTypes.size(); ++i)
    if(isEffect(dspTypes[i]))
      count++;
  return count;
}

int SfzInstrument::HierarchyLevel::getNumProcessorsOfType(OpcodeType type) const
{
  int count = 0;
  for(size_t i = 0; i < dspTypes.size(); ++i) {
    if(dspTypes[i] == type)
      count++; }
  return count;
}

bool SfzInstrument::HierarchyLevel::isLastEffectAmplifier() const
{
  // Figure out the index of the last Amplifier in the dspTypes array, -1 if none:
  int lastAmp = -1;        
  for(int i = (int)dspTypes.size()-1; i >= 0; --i) {
    if(dspTypes[i] == OpcodeType::Amplifier) {
      lastAmp = i; break; }}

  // Figure out the index in the dspTypes array of the last effect that is not Amplifier. We need 
  // the isEffect() test because our dspTypes array stores also the modulators and they shouldn't 
  // count here:
  int lastNonAmpEff = -1; 
  for(int i = (int)dspTypes.size()-1; i >= 0; --i) {
    if(isEffect(dspTypes[i]) && dspTypes[i] != OpcodeType::Amplifier) {
      lastNonAmpEff = i; break; }}

  // Return true, iff the last effect in the chain is an Amplifier. We use < rather than <= to 
  // catch the case when lastAmp == lastNonAmpEff == -1 which happens when there are no effects in 
  // the dspTypes array (it may be empty or there are only other kinds of devices such as 
  // modulators):
  return lastNonAmpEff < lastAmp;
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
  equal &= modRoutings == rhs.modRoutings;
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
  equal &= modRoutings == rhs.modRoutings;
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
  equal &= modRoutings == rhs.modRoutings;
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

rsReturnCode SfzInstrument::setGroupSetting(int gi, Opcode type, float value, int index)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->setSetting(PlaybackSetting(type, value, index));
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::setInstrumentSetting(Opcode type, float value, int index)
{
  global.setSetting(PlaybackSetting(type, value, index));
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::setRegionModulation(int gi, int ri, OpcodeType srcType, int srcIndex, 
  Opcode tgtParam, int tgtIndex, float depth, ModMode modMode)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->regions[ri]->setModulation(
    srcType, srcIndex, tgtParam, tgtIndex, depth, modMode);
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::setGroupModulation(int gi, OpcodeType srcType, int srcIndex, 
  Opcode tgtParam, int tgtIndex, float depth, ModMode modMode)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex; }
  global.groups[gi]->setModulation(srcType, srcIndex, tgtParam, tgtIndex, depth, modMode);
  return rsReturnCode::success;
}

rsReturnCode SfzInstrument::setInstrumentModulation(OpcodeType srcType, int srcIndex, 
  Opcode tgtParam, int tgtIndex, float depth, ModMode modMode)
{
  global.setModulation(srcType, srcIndex, tgtParam, tgtIndex, depth, modMode);
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

rsReturnCode SfzInstrument::removeGroupSetting(int gi, Opcode type, int index)
{
  if(!isGroupIndexValid(gi)) {
    RAPT::rsError("Invalid group index");
    return rsReturnCode::invalidIndex; }
  bool wasRemoved = global.groups[gi]->removeSetting(type, index);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode SfzInstrument::removeInstrumentSetting(Opcode type, int index)
{
  bool wasRemoved = global.removeSetting(type, index);
  if(wasRemoved) return rsReturnCode::success;
  else           return rsReturnCode::nothingToDo;
}

rsReturnCode SfzInstrument::removeRegionModulation(int gi, int ri, OpcodeType modSrcType,
  int modSrcIndex, Opcode modTarget, int modTargetIndex)
{
  if(!isIndexPairValid(gi, ri)) {
    RAPT::rsError("Invalid group- and/or region index");
    return rsReturnCode::invalidIndex; }
  bool wasRemoved = global.groups[gi]->regions[ri]->removeModulation(
    modSrcType, modSrcIndex, modTarget, modTargetIndex);
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

    const std::vector<PlaybackSetting>& settings = lvl->getSettings();
    for(size_t i = 0; i < settings.size(); i++)
      writeSettingToString(settings[i], str);

    const std::vector<ModulationRouting>& routings = lvl->getModRoutings();
    for(size_t i = 0; i < routings.size(); i++)
      writeModRoutingToString(routings[i], str);
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
  // -prepend a <group> tag if the file doesn't start with one, i.e. starts with a naked region
  // -Include stripping away comments (done?). A comment begins with a slash '/' and extends until 
  //  the end of the line. It's really annyoing that we can't yet write any comments. Maybe hava a
  //  general function that removes all characters between a startTag and endTag and call it like
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

  // A vector of opcodes that could not be handled in the first pass. Typcially, these are the 
  // modulation rotuing settings:
  std::vector<std::pair<std::string, std::string>> unhandledOpcodes;

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
    if(ps.getOpcode() != Opcode::Unknown) // It's a regular opcode setting. Handle it immediately.
      lvl->setSetting(ps);
    else                              // It may be a modulation routing opcode. Postpone handling.
      unhandledOpcodes.push_back(std::pair<std::string, std::string>(opcode, value));
  };

  // Sets up the given level according to the given string which is supposed to contain opcode 
  // settings separated by single whitepaces in the format "opocde=value ":
  auto setupLevel = [&](HierarchyLevel* lvl, const std::string& str)
  {
    // In a first pass, we tokenize the string into the "opcode=value" tokens and handle one such
    // token at a time. Some modulation-routing related opcodes cannot be handled in this first 
    // pass. These will be stored in the "unhandledOpcodes" array and may be handled in a second
    // pass:
    unhandledOpcodes.clear();
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

    // In the second pass, we loop through the unhandledOpcodes array to handle those opcodes that 
    // could not be handled in the first pass. The modulation routing opcodes must be handled in a 
    // second pass because some of them require our "settings" and/or "dspTypes" members to be 
    // already fully established which happens in the first pass:
    for(size_t i = 0; i < unhandledOpcodes.size(); i++) {
      ModulationRouting mr = getModRoutingFromString(
        unhandledOpcodes[i].first, unhandledOpcodes[i].second);
      if(mr.isInvalid())
        RAPT::rsError("String could not be parsed as modulation routing");
      else
        lvl->setModulation(mr); }

    // If some ampeg_ opcodes exist, we need to make sure that the AmpEnvGen is routed to an 
    // Amplifier at the end of the effect chain. If no such amplifier exists, i.e. if the last 
    // effect unit is not an Amplifier, an additional amplifier will be inserted. Otherwise, the 
    // existing one will be used:
    lvl->connectAmpEnv();
    //
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

      RAPT::rsAssert(j0 != endOfFile);  
      // for debug - gets triggered when we have empty regions ...but also in other case, i 
      // think

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
  //  error. Maybe instead of using a return code, return a std::string which is empty, if parsing 
  //  works and otherwise contains an error report that we can show the user in the GUI.
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

  SfzCodeBook* cb = SfzCodeBook::getInstance();
  s += cb->opcodeToString(op, idx) + std::string("=") + cb->valueToString(op, val)  + "\n";

  // ToDo:
  // -Document why the lokey, hikey, lovel, hivel opcodes are not handled here. I think, it's 
  //  because they are handled already by the caller because they require special treatment.
}

void SfzInstrument::writeModRoutingToString(const ModulationRouting& r, std::string& s)
{
  // Retrieve data from routing object:
  using OT = OpcodeType;
  using OC = Opcode;
  using MM = ModMode;
  OT  srcTp  = r.getSourceType();
  OT  tgtTp  = r.getTargetType();
  OC  tgtOp  = r.getTargetOpcode();
  int srcIdx = r.getSourceIndex();
  int tgtIdx = r.getTargetIndex();
  MM  mode   = r.getMode();

  // Handle special cases:
  if(srcTp == OT::AmpLfo && tgtOp == OC::volumeN && mode == MM::absolute)
  {
    RAPT::rsAssert(srcIdx == 1, "There should be at most one AmpLfo");
    s += "amplfo_depth=" + rsFloatToString(r.getDepth()) + "\n";
    return;
  }
  if(srcTp == OT::AmpEnv && tgtOp == OC::amplitudeN && mode == MM::absolute)
  {
    //return;
    // Uncommenting makes unit test fail. I think, we need to connect the ampeg at the end of 
    // setFromSFZ or something like that

    // The ampeg is always implicitly routed to the last Amplifier in the chain.
    // Maybe an additional constraint should be that tgtIdx is the index of the last Amplifier such
    // that routings of the ampeg_ to other amplifiers will be stored. Currently we may throw away
    // too many ampeg routings
  }
  if(srcTp == OT::FilterEnv && tgtOp == OC::cutoffN && mode == MM::cents)
  {
    RAPT::rsAssert(srcIdx == 1, "There should be at most one FilterEnv");
    // ...because in the sfz spec, the fileg_ opcodes are not indexed, i.e. there's no such thing
    // as fileg2_cutoff3 or anything like that. There is just fileg_depth which we interpret as
    // fileg1_cutoff1 at the moment. Maybe later, we should interpret it as fileg1_cutoffAll but
    // for that, we will first generally have to support such an "All" syntax.

    RAPT::rsAssert(tgtIdx == 1);
    // Hmm...we currently interpret the fileg_depth opcode as applying to the first cutoff, i.e. 
    // the cutoff of the first filter in the chain. But maybe it should apply to all filters? But 
    // then we should probably only store it once in the sfz string...hmm...

    s += "fileg_depth=" + rsFloatToString(r.getDepth()) + "\n"; 
    return;
  }
  // ToDo: Add other special cases here one by one. For example, fillfo, pitcheg, amplfo, etc.
  // Maybe factor out a function for handling the special cases, returning true, if the case was 
  // handled, false if not and call it here...or maybe not...

  // Handle the general case:
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  std::string tmp;
  tmp  = cb->modSourceToString(srcTp, srcIdx) + "_";
  tmp += cb->modTargetToString(tgtTp, tgtIdx, tgtOp) + "=";
  tmp += cb->modDepthToString(r.getDepth(), mode, tgtOp) + "\n";
  s += tmp;
}

PlaybackSetting SfzInstrument::getSettingFromString(
  const std::string& opStr, const std::string& valStr)
{
  using PS  = PlaybackSetting;
  using PST = Opcode;
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  int idx;
  PST   op  = cb->stringToOpcode(opStr, &idx);
  float val = cb->stringToValue(op, valStr);
  return PS(op, val, idx);
}
// maybe move into a static "fromString" member function of PlaybackSetting

ModulationRouting SfzInstrument::getModRoutingFromString(
  const std::string& opStr, const std::string& valStr)
{
  SfzCodeBook* cb = SfzCodeBook::getInstance();

  // Split off the substrings for the modulation source and target by using the underscore as 
  // seperator:
  std::string sep("_");
  int start  =  0;
  int length = -1;  // initial value should not matter
  rosic::rsFindToken(opStr, sep, &start, &length);
  std::string srcStr = opStr.substr(start, length);
  std::string tgtStr = opStr.substr(length+1, opStr.length() - length - 1);

  // Figure out type and index of the modulation source:
  int srcIndex;
  OpcodeType srcType = cb->stringToModSource(srcStr, &srcIndex);

  // Figure out Opcode and index of the modulation target:
  int tgtIndex;
  Opcode tgtOpcode = cb->stringToOpcode(tgtStr, &tgtIndex);
  //OpcodeType targetType = cb->getOpcodeType(targetOpcode);  // not needed

  // Figure out modulation depth and mode:
  ModMode mode;
  float depth = cb->stringToModDepth(valStr, &mode, tgtOpcode);

  return ModulationRouting(srcType, srcIndex, tgtOpcode, tgtIndex, depth, mode);

  // ToDo:
  // -Maybe this whole code should go into SfzCodeBook. It seems to fit in there better than here.
  //  But the codebook doesn't know about class ModulationRouting, so how should we return the 
  //  data? Maybe a function that takes all the members of that class as pointer arguments? But maybe
  //  it's ok to keep the code here

  //return ModulationRouting();  // standard constructor will create an invalid object
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
-when the sfz file immediately starts with a <region> without a <group> before, the parsing 
 fails. We should probably prepend a group statement if there is none in the pre-processing stage.
-patches starting with a comment don't load - the comment is not stripped off 
 -> bug in pre-processor?
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
