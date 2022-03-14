namespace rosic {
namespace Sampler {

//=================================================================================================
// rsSamplerEngine::SignalProcessorChain

void EffectChain::processFrame(float* L, float* R)
{
  for(size_t i = 0; i < processors.size(); i++)
    processors[i]->processFrame(L, R);
}

void EffectChain::processBlock(float* L, float* R, int N)
{
  for(int n = 0; n < N; n++)
    processFrame(L, R);
}

size_t EffectChain::getNumEffects(OpcodeType type) const
{
  size_t count = 0;
  for(size_t i = 0; i < processors.size(); i++) {
    if(processors[i]->getType() == type)
      count++;
  }
  return count;
}

Effect* EffectChain::getEffect(OpcodeType type, int index)
{
  RAPT::rsAssert(index >= 1 || index == -1);
  index = RAPT::rsMax(index-1, 0);
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int)processors.size(); i++) {
    Effect* dsp = getEffect(i);
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++;
    }
  }
  return nullptr;
}

//=================================================================================================
// SamplePlayer

/** Counts the number of occurences of elem in array a of length N. */
template<class T>
int rsCount(const T* a, int N, T elem)
{
  int c = 0;  // counter
  for(int i = 0; i < N; i++)
    if(a[i] == elem)
      c++;
  return c;
}
// move into rsArrayTools


// make (static) member of SamplePlayer ...maybe it should take a vector of Processor* and then we
// can use it also instead of effectChain.getNumEffects(opType) to match both branches more 
// closely. Then, we need to rename it
// rename to getNumProcessorsOfType(const std::vector<Processor*>& processors, OpcodeType type):
size_t getNumModulators(const std::vector<Processor*>& modSources, OpcodeType type)
{
  size_t count = 0;
  for(size_t i = 0; i < modSources.size(); i++) {
    if(modSources[i]->getType() == type)
      count++;
  }
  return count;
}

int findProcessorIndex(Processor* processors, int numProcessors, OpcodeType type, int index)
{
  // ToDo: 
  // -Check, if index < 0 and if so, modify it processors.size() + abs(index) to use indices 
  //  -1,-2,-3 for the hardwired modulators for amp, cutoff, pitch (maybe -4 for cutoff2)

  int count = 0;
  for(int i = 0; i < numProcessors; ++i) {
    if(processors[i].getType() == type) {
      count++;
      if(count == index)
        return i; }}
  return -1;
}


Processor* findProcessor(Processor* processors, int numProcessors, OpcodeType type, int index)
{
  int i = findProcessorIndex(processors, numProcessors, type, index);
  if(i == -1)
    return nullptr;
  else
    return &processors[i]; 

  /*
  // old:
  int count = 0;
  for(int i = 0; i < numProcessors; ++i) {
    if(processors[i].getType() == type) {
      count++;
      if(count == index)
        return &processors[i]; }}
  return nullptr;
  */
}






bool SamplePlayer::augmentOrCleanProcessors(const std::vector<OpcodeType>& dspTypeChain)
{
  //SfzCodeBook* cb = SfzCodeBook::getInstance();

  for(int i = 0; i < (int)dspTypeChain.size(); i++)
  {
    // Figure out the type of the effect or modulator that may need to be added to the effectChain 
    // or modSources array:
    OpcodeType opType = dspTypeChain[i];

    // We either add an effect to the effectChain or a modulator to the modSources array:
    if(SfzCodeBook::isEffectSetting(opType))
    {
      // Figure out, if we actually need to add another effect to the chain. If not, there's 
      // nothing more to do in this iteration:
      int sfzIndex = rsCount(&dspTypeChain[0], i, opType) + 1;
      if(effectChain.getNumEffects(opType) >= sfzIndex)
        continue;

      // OK - now we actually need to grab another effect of given type from the pool:
      Effect* eff = getEffect(opType);
      if(eff)
      {
        eff->setParametersToDefaults(sfzIndex);
        effectChain.addEffect(eff);
      }
      else 
      {
        disassembleProcessors();
        return false;
        // Not enough effects of desired type are available in the pool so we roll back any partially 
        // built chain and report failure. 
      }
    }
    else if(SfzCodeBook::isModSourceSetting(opType))
    {
      // The logic for adding modulation sources is the same as for adding effect processors:
      int sfzIndex = rsCount(&dspTypeChain[0], i, opType) + 1; 
      if(getNumModulators(modSources, opType) >= sfzIndex)
        continue;
      Processor* mod = getModulator(opType);  // maybe use a general getProcessor function
      if(mod) {
        mod->setParametersToDefaults(sfzIndex);
        modSources.push_back(mod);     }
      else {
        disassembleProcessors();
        return false;  }
      // Try to get rid of the duplication by making the branches even more similar and then try to
      // refactor....
    }
  }
  return true;
  // When we arrive here, we have successfully finished the loop which means that we either could
  // add enough effects/modulators of the desired types to the chain/array or they were already 
  // present before. In both cases, the effectChain or modSources is now in the required state so 
  // we can report success.
  //
  // Maybe remove the calls to setParametersToDefaults - instead we should perhaps rest them to
  // defaults when the object is reposited - but: when the default value depends on the index as
  // in the eqN_freq opcode, we don't really know the index when the object is just sitting in the
  // pool, so maybe it's indeed more appropriate to do it here. Maybe do both? The rationale is to
  // have all objects in the pool in a well defined state. That's not actually important though but
  // it somehow seems cleaner.
}

bool SamplePlayer::assembleModulations(const std::vector<ModulationSetting>& modSettings)
{
  RAPT::rsAssert(dspPool);
  RAPT::rsAssert(modMatrix.empty(), "Someone has not cleaned up the modMatrix");

  // Check, if enough connectors are available:
  if(dspPool->getNumIdleConnectors() < (int)modSettings.size())
    return false;

  // Do the assembly:
  for(size_t i = 0; i < modSettings.size(); i++)
  {
    ModulationSetting    ms = modSettings[i];
    ModulationConnector* mc = dspPool->grabConnector();
    RAPT::rsAssert(mc); // we already verified that enough are available


    // Determine pointer to modulation source and its index in our modSources array and set it up 
    // in the connector:
    /*
    //old:
    Processor* prc = findProcessor(
      modSources[0], (int) modSources.size(), ms.getSourceType(), ms.getSourceIndex());
    Modulator* src = dynamic_cast<Modulator*> (prc);
    RAPT::rsAssert(src);
    mc->setSource(src);
    // ToDo: instead of storing the pointer, store the index into the modSources array 
    // -> findProcessor should return an int
    */
    // new:
    int j = findProcessorIndex(modSources[0], (int) modSources.size(), ms.getSourceType(), 
      ms.getSourceIndex());
    RAPT::rsAssert(j >= 0);
    Processor* src = modSources[j];
    RAPT::rsAssert(src);
    mc->setSourceIndex(j);
    mc->setSource(src);


    // Determine pointer to modulation target (Processor and Parameter) and set it up in the
    // connector:
    Processor* prc;
    if(SfzCodeBook::isModSourceSetting(ms.getTargetType())) {
      prc = findProcessor(modSources[0], (int) modSources.size(), // Receiver is another modulator
        ms.getTargetType(), ms.getTargetIndex());  }
    else {
      prc = findProcessor(effectChain.processors[0], (int) effectChain.processors.size(), 
        ms.getTargetType(), ms.getTargetIndex());  }              // Receiver is an effect
    RAPT::rsAssert(prc);
    Parameter* param = prc->getParameter(ms.getTargetOpcode());
    RAPT::rsAssert(param);
    mc->setTarget(prc, param);
    RAPT::rsAppendIfNotAlreadyThere(modTargetProcessors, prc);
    RAPT::rsAppend(modTargetParams, param);
    // ToDo: instead of storing the pointer, store the index -> findProcessor should return an int
    // -> append should return the index...maybe not...i think, we need that index only for the 
    // modSource

    // Set up modulation depth and mode and add the connection to the modMatrix:
    mc->setDepth(ms.getDepth());
    mc->setMode(ms.getMode());
    modMatrix.push_back(mc);
  }

  return true;

  // Assembling the mod-connections may also fail if not enough are available - in this 
  // case we also need to roll back and return false. But: maybe we should do it in a way that 
  // can't fail by not grabbing pre-allocated connection objects from the pool but rather using a
  // std::vector<ModulationConnection> instead of std::vector<ModulationConnection*> Or maybe the
  // Region/Group etc. object could maintain such an array itself such we do not need to assemble 
  // it at all...or only need to re-connect pins, i.e. update the source/target pointers...but no -
  // I'm confusing again Regions with Layers here - we may have several layers playing the same
  // region and they will need different pointers. But nevertheless, it may make sense to use
  // a vector of direct ModulationConnection objects rather than using pre-allocated ones from the
  // pool. This may simplify the code but it may make the memory occupied by the modMatrix larger 
  // because now it stores objects instead of pointers...but this may actually help with caching. 
  // We would need less pointer chasing during playback.
}

bool SamplePlayer::assembleProcessors(
  const std::vector<OpcodeType>& dspTypes, const std::vector<ModulationSetting>& modSettings) 
{
  if(!areProcessorsEmpty()) {    // Sanity check
    RAPT::rsError("Someone has not cleaned up after finishing playback!");
    disassembleProcessors();  }  // ...so we do it here. But this should be fixed elsewhere!

  if(!augmentOrCleanProcessors(dspTypes)) 
    return false;

  if(!assembleModulations(modSettings)) {
    disassembleProcessors();
    return false; }

  return true;
}

void SamplePlayer::disassembleProcessors()
{
  for(int i = 0; i < effectChain.getNumEffects(); i++)
    dspPool->repositEffect(effectChain.getEffect(i));
  effectChain.clear();

  for(size_t i = 0; i < modSources.size(); ++i)
    dspPool->repositModulator(modSources[i]);
  modSources.clear();

  for(size_t i = 0; i < modMatrix.size(); ++i) {
    modMatrix[i]->reset();
    dspPool->repositConnector(modMatrix[i]); }
  modMatrix.clear();


  modTargetProcessors.clear();

  // ToDo: 
  // -let effectChain just be std::vector<Effect> to handle it uniformly with the modulators.
  //  Member functions should become free functions (maybe wrapped into a class...maybe this class)
  // -benchmark whether its faster to traverse the array from the back
}

void SamplePlayer::setupProcessorSetting(const PlaybackSetting& s)
{
  Effect* dsp = effectChain.getEffect(s.getTargetOpcodeType(), s.getIndex());
  if(dsp != nullptr)
    dsp->setParameter(s.getOpcode(), s.getValue());
  else
    RAPT::rsError("No processor available for DSP opcode");
    // We could not find a suitable processor in our dspChain to which the given setting could be
    // applied. If this happens, something went wrong (i.e. we have a bug) in assembleDspChain or 
    // dspChain.getProcessor.
}

void SamplePlayer::setupDspSettings(const std::vector<PlaybackSetting>& settings,
  double sampleRate, RegionPlayer* rp, bool busMode, PlayStatus* iv)
{
  SfzCodeBook* codebook = SfzCodeBook::getInstance();
  for(size_t i = 0; i < settings.size(); i++)
  {
    PlaybackSetting s = settings[i];
    Opcode op = s.getOpcode();
    if(     codebook->isEffectSetting(op))    { setupProcessorSetting(s);              }
    else if(codebook->isPlayerSetting(op)) { setupPlayerSetting(s, sampleRate, rp, iv); }
  }
}

//=================================================================================================
// RegionPlayer

rsReturnCode RegionPlayer::setRegionToPlay(const Region* regionToPlay,
  const AudioFileStream<float>* sampleStream, uchar key, uchar vel, double fs, bool busMode, 
  PlayStatus* iv)
{
  releaseResources(); // actually, it should not hold any at this point - or should it?
  region = regionToPlay;
  if(region == nullptr)
    return rsReturnCode::nothingToDo;
  stream = sampleStream;
  return prepareToPlay(key, vel, fs, busMode, iv);
}

void RegionPlayer::processFrame(float* L, float* R)
{
  // Negatively initialized sampleTime implements delay. If we are still "waiting", we just 
  // increment the time and return 0,0. Actual output will be produced as soon as sampleTime 
  // reaches zero. It actually sucks to have to do this per sample just to support offset. Try to
  // implement offset differently...but how?
  if(sampleTime < 0.0) {
    sampleTime += 1.0;
    if(sampleTime >= 0.0)
      sampleTime += offset;
    *L = *R = 0.f;
    return;  }
  if(sampleTime == 0.0)
    sampleTime = offset;

  stream->getFrameStereo((float)sampleTime, L, R);  
  // try to avoid the conversion to float - use a 2nd template parameter for the time



  //---------------------------------------------
  // Under construction - handle modulations:

  int numSources = (int) modSources.size();
  std::vector<float> modBuffer;
  modBuffer.resize(2*numSources);
  // preliminary - should become a member of PlayStatus, the size should be pre-allocated according
  // to the maximum number of modulators that a region/Group/Instrument has


  // Update our modulators:
  for(size_t i = 0; i < modSources.size(); ++i)
  {
    //modSources[i]->updateModValue();
    modSources[i]->processFrame(&modBuffer[2*i], &modBuffer[2*i+1]);
    // ToDo: Try to use processFrame instead. But then we need to store the output frames of all 
    // modulators in a local buffer here (maybe use a member modBuffer) and the 
    // ModulationConnection must somehow maintain an index into that buffer. Maybe the 
    // ModualtionConnection could store array indices into our modSources, 
    // modTargetProcessors arrays instead of pointers to Modulator, Processor. The index into the
    // modSources array could then double as index into the buffer of stored modulation signals.
    // ...done: we currently hold the index and the pointer - that's redundant because the pointer
    // could be obtained from the index - but it may save one indirection in per-sample processing
    // to store it redundantly...or does it...we'll see....if not, keep only the index.
    //
    // Actually, if we assume a single audio thread, the modBuffer could and should be shared among
    // all the RegionPlayers (if multi-threading is added later, it could be declared 
    // thread_local). Maybe the PlayStatus object could be an appropriate place to hold such a 
    // buffer. It could also generally hold the signal buffers needed for block processing.
  }

  // Initialize modulated parameters to non-modulated values:
  for(size_t i = 0; i < modTargetParams.size(); ++i)
    modTargetParams[i]->initModulatedValue();

  // Apply the modulations to the affected parameters:
  for(size_t i = 0; i < modMatrix.size(); ++i)
  {
    ModulationConnector* con = modMatrix[i];
    Parameter* par = con->getTargetParam();
    int   si = con->getSourceIndex();
    float u  = par->getValue();                 // unmodulated value
    //float m  = modSources[si]->modValue;        // modulator output
    float m  = modBuffer[2*i];                  // modulator output
    float c  = con->getContribution(m, u);      // compute modulation contribution
    par->applyModulation(c);                    // accumulate the contribution
    int dummy = 0;
    // We could avoid one dereferencing by avoiding the si variable and instead directly use a
    // con->getSource() function directly returning the pointer...but we want to refactor this 
    // later in a way such that the output of the modSource is not stored in the modSource object 
    // but rather in a modBuffer - for this we will need the index si because the same index will
    // be used for this buffer...
    // Here, we should perhaps also keep track of whether or not the modulatedValue in the target
    // parameter actually changes and if so, set a dirty flag - and avoid subsequent update 
    // compuations, if it didn't change (because all modulators have produced the same output as 
    // in the previous sample - which is common in envelope sustain and for midi based modulators).
    // But maybe this dirtification should be handled by Parameter itself.
  }

  // Let the affected Processors update their algo-parameters according to the new user parameters
  // now containing the modulations:
  for(size_t i = 0; i < modTargetProcessors.size(); ++i)
  {
    //if(modTargetProcessors[i].isDirty())
    //  modTargetProcessors[i].updateCoefficients();
    // Maybe updateCoefficients should itself check if coeffs are dirty..but maybe not because this
    // is a boilerplate function and we want to keep the amount of boilerplate low. 
    // updateCoefficients should be a purely virtual member function and do the stuff that we 
    // currently do in prepareToPlay - we can call it from there and maybe make prepareToPlay
    // non-virtual
  }

  // End of modulation handling. Maybe factor this out into a member function of the baseclass 
  // because we may need it in GroupPlayer etc. too in order to apply group modulations to the
  // group effects.
  //---------------------------------------------




  // Update our sampleTime counter:
  sampleTime += increment;
  if(loopMode == LoopMode::loop_continuous)
  {
    if(sampleTime >= loopEnd)
      sampleTime -= (loopEnd - loopStart);
  }

  // Apply the effect chain:
  effectChain.processFrame(L, R);


  // ToDo:
  // -implement better interpolation methods (sinc, elephant, ...)
  // -keep sample time as combination of int and float to avoid computation of the fractional part
  //  at each sample and to avoid losing precision for the fractional part when the integer part
  //  is large (thereby eating up significant digits).
  // -the interpolation should probably be handled by the AudioStream class to make it re-usable
  //  also for resampled audio playback in other contexts. Maybe here, we should just call
  //  stream->getFrameStereo(sampleTimeInt, sampleTimeFrac, &L, &R);
  // -it should probably also receive the increment in order to make a decision for time-scaling
  //  the sinc, if necessary for anti-aliasing. Maybe it should implement various anti-aliasing
  //  algorithms: sinc-interpolation, mip-mapping, integrate -> interpolate -> differentiate, 
  //  oversample etc. ... and maybe combinations of them...maybe have an opcode resample=sinc
  //  etc....but maybe handle that in the xml
  // -apply pitch envelope and lfo - these should affect (scale?) the effective increment that we 
  //  add to sampleTime - but our increment *member* should not be modified, instead, do something
  //  like sampleTime += increment * incScaler; where incScaler is computed from the pitch 
  //  modifiers. Maybe create a subclass of SignalProcessor called PitchShifter that has just a
  //  dummy callback that just stores the desired shift value and we read it out here
  // -Implement loop_sustain and all the other modes. The one_shot mode requires us to hanlde
  //  noteOffs differently (the handling is preliminary anyway)
}

void RegionPlayer::processBlock(float* L, float* R, int N)
{
  for(int n = 0; n < N; n++)
    processFrame(&L[n], &R[n]);
    // Preliminary - ToDo: implement actual proper block processing using signal buffers. These
    // buffers should be shared among all the player objects. Perhaps the best place for them is
    // in the PlayStatus object.
}

bool RegionPlayer::isPlayable(const Region* region)
{
  bool ok = true;
  ok &= region != nullptr;
  ok &= region->getGroup() != nullptr;
  ok &= region->getGroup()->getInstrument() != nullptr;
  ok &= region->getCustomPointer() != nullptr;           // should point to a stream object
  return ok;
}

void RegionPlayer::releaseResources()
{
  RAPT::rsAssert(dspPool, "This pointer should be assigned soon after creation");
  disassembleProcessors();  // Move effects, modulators and connectors back into the pool
  stream = nullptr;         // Invalidate our pointers
  region = nullptr;
  //key = 0;                // shouldn't we do this?
}

void RegionPlayer::allocateMemory()
{
  modSources.reserve(8);
  //modTargets.reserve(8);
  modMatrix.reserve(32);
  effectChain.reserve(8);
  // These values are ad-hoc arbitrarily chosen. Maybe give the function some parameters to choose
  // how many of each should be pre-allocated. It should be enough to avoid having to allocate more
  // later
}



rsReturnCode RegionPlayer::prepareToPlay(uchar key, uchar vel, double fs, bool busMode, 
  PlayStatus* iv)
{
  RAPT::rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  RAPT::rsAssert(stream != nullptr);   // Ditto.

  this->key = key;

  if(!assembleProcessors(busMode)) {
    releaseResources();
    return rsReturnCode::layerOverload; }

  resetPlayerSettings();
  setupDspSettingsFor(region, fs, busMode, iv);

  // new:
  if(!modSources.empty())
    prepareToPlay1((Processor**) &modSources[0], (int) modSources.size(), key, vel, fs);
  if(!effectChain.processors.empty())
    prepareToPlay1((Processor**) &effectChain.processors[0], (int) effectChain.processors.size(), 
      key, vel, fs); 
  // This should be cleaned up: Perhaps we should get rid of the two subclasses Effect and 
  // Modulator of Processor and let modSources and effectChain both be vectors of Processor. Then
  // prepareToPlay1 can take a refecrence to a std::vector<Processor*> and inside it, we could use
  // a range-based loop. See also comment below class Modulator.

  // old:
  //effectChain.prepareToPlay(key, vel, fs);
  // Should be replaced by:
  //   prepareToPlay(modSources,  key, vel, fs);
  //   prepareToPlay(effectChain, key, vel, fs);
  // prepareToPlay should take a std::vector<Processor> as input. The rationale for preparing the
  // modulators first is that their initial output may already affect the initial parameters of
  // the effects? ...but does that matter? Aren't the params recomputed in processFrame anyway?
  // ...perhaps it doesn't matter, but the order feels right this way anyway. This prepareToPlay
  // function may be a member of SamplePlayer.


  return rsReturnCode::success;
  // Overload should actually not happen in therory (as by the sfz spec, and unlimited number of 
  // layers is available), but in practice, it may happen in extreme situations like triggering a
  // whole lot of layers at once or in very short succession while already being close to the 
  // limit, such that we don't have enough pre-allocated players and/or dsp objects available and 
  // the required additional allocation of more is not fast enough. In such a case, we don't want 
  // to play the region at all. The caller should clean up and discard the RegionPlayer object.
}

bool RegionPlayer::hasFinished()
{
  //int numFrames = stream->getNumFrames();
  //int tmp = stream->getNumFrames() - 1;
  //if( sampleTime >= stream->getNumFrames() )  // old
  if(sampleTime >= endTime)                   // new
    return true;

  // ToDo:
  // -Make sure that this function is inlined - it's called per sample.
  // -Check also, if the amplitude envelope has reached its end and/or filter and interpolator 
  //  ringout is finished. Maybe we need a function dspChain.getRingoutTime() that we can call
  //  on initialization. Maybe endTime should not be a member after all. not sure yet
  // -Maybe, if we allow the frequency envelope to go negative, we could also move 
  //  backward through the sample, so having reached the end of the stream may not actually be an
  //  appropriate condition. Or maybe, we should allow more general time-warping envelopes. 
  //  We'll see

  return false;
}

bool RegionPlayer::assembleProcessors(bool busMode)
{
  const Region* reg = region;
  const Group*  grp = reg->getGroup();
  const Global* glb = grp->getInstrument();

  // The DSPs for which the region itself defines settings/opcodes are always needed:
  bool ok = SamplePlayer::assembleProcessors(
    reg->getOpcodeTypeChain(), reg->getModulationSettings());
  if(!ok)
    return false;

  // If we are not in busMode, the enclosing group and/or enclosing instrument settings act as
  // fallback values for the region so we may require additional DSPs to apply these opcodes
  // to the region, too:
  if(!busMode) {
    if(!augmentOrCleanProcessors(grp->getOpcodeTypeChain())) {
      return false; }
    if(!augmentOrCleanProcessors(glb->getOpcodeTypeChain())) {
      return false; }}

  return true;
  // OK - everything went well so we report success. If, on the other hand, false is returned, it
  // means we do not have enough processors of the required types available. In this case, the 
  // caller should roll back and discard the whole RegionPlayer object. We either play a region 
  // correctly or not at all. This is an error condition that could conceivably arise in normal 
  // usage (because we did not pre-allocate enough DSPs), so we should be able to handle it 
  // gracefully. It should actually not happen, i.e. we should make sure to always pre-allocate 
  // enough DSPs but it may be impractical to ensure in a 100% bulletproof manner. But let's 
  // try at least to make that an exception that occurs only in extreme scenarios.
}

void RegionPlayer::resetPlayerSettings()
{
  // Initialize all values and DSP objects to default values (maybe factor out):
  sampleTime = 0.0;
  increment  = 1.0;
  loopStart  = 0.0;
  loopEnd    = 0.0;
  loopMode   = LoopMode::no_loop;
  offset     = 0.f;
  //tune       = 0.f;
  //transpose  = 0.f;
  endTime    = (float)stream->getNumFrames();
  // Maybe use -1? That may require updating the unit tests. But maybe it's appropriate to use 
  // numFrames when assuming linear interpolation. I think, for general interpolators, we should 
  // use endTime = numFrames - 1 + kernelSize/2. Test this with very high downshifting factors and
  // maybe with a sample that is just 1 sample long with a value of 1. We should see the 
  // interpolation kernel as output.

  // What about key?
  // key = 0;   // uncommenting breaks unit tests - figure out why and document
}

void RegionPlayer::setupDspSettingsFor(const Region* r, double fs, bool busMode, 
  PlayStatus* iv)
{
  // To set up the settings, we call setupDspSettings 3 times to:
  //   (1) set up the general instrument-wide settings
  //   (2) set up group specific settings (this may override instrument settings)
  //   (3) set up region specific settings (this may override group and/or instrument settings)
  // but only if not in bus-mode in which case the group and instrument settings apply to separate
  // DSP objects on the groups which map to sub-busses and the instrument which maps to the final
  // master bus:
  if(!busMode) {
    setupDspSettings(region->getGroup()->getInstrument()->getSettings(), fs, this, busMode, iv);  
    setupDspSettings(region->getGroup()->getSettings(), fs, this, busMode, iv); }
  setupDspSettings(region->getSettings(), fs, this, busMode, iv);
  // Works for accumulative setting but not for those that should awlays override

  // Compute the final member variables from the intemediates:
  setupFromIntemediates(*iv, fs);
  // maybe call this only when not in busMode because in busMode, it will get called again later 
  // and we don't want to call it twice. It will not do anything wrong, but it does the computation
  // at the first time for nothing....
}

void RegionPlayer::setupFromIntemediates(const PlayStatus& iv, double fs) // fix typo!
{
  // This may get called twice on noteOn in busMode -> verify and try to avoid the first call

  // Compute the per-sample increment:

  // old:
  float pitch_keycenter = region->getSettingValue(Opcode::PitchKeyCenter, -1, false);
  double pitchOffset = 0.01 * iv.pitch_keytrack * (double(key) - pitch_keycenter) 
    + iv.transpose + 0.01 * iv.tune; 

  // new:
  //double pitchOffset = 0.01 * iv.pitch_keytrack * (double(key) - iv.pitch_keycenter) 
  //  + iv.transpose + 0.01 * iv.tune;

  increment = pow(2.0, pitchOffset / 12.0) * stream->getSampleRate() / fs;
  // The formula using rsPitchOffsetToFreqFactor is too imprecise: when we have a pitchOffset of 
  // exactly -12, for example, we want the increment be multiplied by exactly 0.5, but using this
  // function, the factor comes out as 0.49999..., so we use the more precise (and more 
  // expensive) call to pow. It's not per-sample here code anyway, so we may afford that.

  // Sanity-check the loop boundaries:
  loopStart = RAPT::rsMax(loopStart, 0.0);
  loopEnd   = RAPT::rsClip(loopEnd, loopStart, (double)stream->getNumFrames());

  // ToDo:
  // -Take into account pitch_veltrack...maybe we need to take an addition vel parameter or we
  //  add such a fields to PlayStatus...which should be extended to a "MidiStatus" anyway.
  // -Verify the formulas

  // -Can we avoid the inquiry for the rootKey? This may be a bit expensive. Maybe the RegionPlayer
  //  should have a rootKey member? That may be messy to deal with in busMode. Maybe try to treat
  //  the pitch_keycenter opcode somehow like the transpose and pitch_keytrack opcode. What would 
  //  actually be a meaningful "accumulative" behavior of a setting like that? Maybe accumulate 
  //  differences to the default or something? Maybe keep a rootKeyShift member that is by default 
  //  zero, accumulates th differences of the pitch_keycenter opcode with respce to 60 (the 
  //  default)? Try this and write unit tests..
}

void RegionPlayer::setupPlayerSetting(const PlaybackSetting& s, double sampleRate, 
  RegionPlayer* rp, PlayStatus* iv)
{
  RAPT::rsAssert(rp == this);
  // For RegionPlayer objects like this, this function is supposed to be called only for the object
  // itself. The rp pointer is needed here only to conform to the baseclass interface. For 
  // GroupPlayers and InstrumPlayers, the pointer is supposed to hold the RegionPlayer to which 
  // this setting should be applied accumulatively in busMode. In default mode, GroupPlayer and 
  // InstrumPlayer play no role at all.

  float val = s.getValue();
  int   N   = s.getIndex();
  using OC  = Opcode;
  switch(s.getOpcode())
  {
  // Player:
  //case OC::PitchKeyCenter: { iv->pitch_keycenter = val;           } break; // was done by caller before
  case OC::PitchKeyTrack:  { iv->pitch_keytrack  = val;           } break;
  case OC::Transpose:      { iv->transpose       = val;           } break;
  case OC::Tune:           { iv->tune            = val;           } break;
  case OC::Delay:          { sampleTime    = -val * sampleRate;   } break;
  case OC::Offset:         { offset        =  val;                } break;
  case OC::LoopMode:       { loopMode      = (LoopMode)(int) val; } break;
  case OC::LoopStart:      { loopStart     =  val;                } break;
  case OC::LoopEnd:        { loopEnd       =  val;                } break;

  // Tracking:
  case OC::ampN_veltrack: { 
    iv->ampN_veltrack[N] = val;          } break;
    // ToDo: make sure that ampN_veltrack has large enough size!!!
  }

  // ToDo: maybe handle PitchKeyCenter just like transpose and tune...hmmm.not sure...

}
// ...actually, if we can assume that all values start at neutral values, we can always accumulate
// and the distinction between this implementation and the one in SampleBusPlayer becomes identical
// and can therfore be moved into the baseclass..or..well - not quite: SampleBusPlayer applies 
// everything to rp and we to "this" - but "this" equals rp...hmmm..but no - that accumulation 
// business thwarts the override behavior! We receive calls from (maybe) the instrum, then (maybe) 
// the group, then (maybe) the region. Each of these is optional but when a call happens, it should 
// override the current setting. 

// Maybe tune and transpose do not need to be members of RegionPlayer if we give the function 
// another parameter which is a pointer to a struct that holds all the temporary intermediate
// values. This will become more relevant when more and more parameters are controlled by opcodes
// which all influence a final value

//=================================================================================================
// SampleBusPlayer

void SampleBusPlayer::setupPlayerSetting(const PlaybackSetting& s, double fs, 
  RegionPlayer* rp, PlayStatus* iv)
{
  // We are supposedly a higher level player object like GroupPlayer or InstrumentPlayer but the 
  // setting in question reaches through to an underlying (embedded) RegionPlayer object and must 
  // be accumulated into the respective value there. This applies to all kinds of settings that 
  // cannot be meaningfully achieved by post-processing in the effect chain but instead must be 
  // applied directly at the playback source. Some of them (like delay) *could* be achieved 
  // also by post-processing with a suitable effect but it's more efficient to do it at the source.
  // Others like all tuning related stuff indeed need to be done at the source.

  RAPT::rsAssert(rp != nullptr);
  float val = s.getValue();
  int   N   = s.getIndex();
  using OC  = Opcode;
  switch(s.getOpcode())
  {
  // Player:
  /*
  case OC::PitchKeyCenter: 
  { 
    // keycenter needs different handling - we can't just add values to accumulate
    float delta = val - iv->pitch_keycenter;
    iv->pitch_keycenter += delta;
    // todo: verify, if that works as intended - define pitch_keycenter only for the group and see, 
    // if that works as intended - i'm not sure, if that makes sense at all
  } break;
  */

  case OC::PitchKeyTrack: { iv->pitch_keytrack += val;        } break;
  case OC::Transpose:     { iv->transpose      += val;        } break;
  case OC::Tune:          { iv->tune           += val;        } break;
  case OC::Delay:         { rp->sampleTime     += -val * fs;  } break;
  case OC::Offset:        { rp->offset         += val;        } break;

    // Tracking:
  case OC::ampN_veltrack: { 
    iv->ampN_veltrack[N] += val; } break;
    // !!!!!  TODO: make sure that ampN_veltrack has large enough size !!!!!!
    // we may have to resize the array on sfz-load
  }

  // Maybe the default branch should call rp->setupPlayerSetting(s, sampleRate, val). That would 
  // revert the behavior to override mode which may make most sense for most settings - certainly
  // for stuff like loop_mode, loop_start, etc. - there is no meaningful way for this setting to
  // accumulate

  // hmmm...the PitchKeyCenter stuff doesn't work. the probelm is that in busMode, we call
  // RegionPlayer::setupPlayerSetting for the region first and *then* call 
  // SampleBusPlayer::setupPlayerSetting for the Group and Global, so if we would just set the
  // value here, we would override the region's settings - exactly the opposite of what should 
  // happen. Maybe we could fix this always calling in RegionPlayer::setupDspSettingsFor all 3
  // setup functions without using the if(!busMode) conditional. Then, for any setting that should 
  // always override (and never accumulate), we just leave out the respective branch here. ...We
  // really need to call all 3 in any case to make sure that a pitch_keycenter defined in a group 
  // becomes effective as fallback value, even when in busMode...but what if for some opcode the
  // region has no definition? Then the group-value will be applied twice, i.e. it will accumulate
  // to itself...maybe we need a second pass. seems like in busMode, we need the following call 
  // sequence:
  //   (1) set all parameters to their defaults
  //   (2) set instrum settings that should always override
  //   (3) set group settings that should always override
  //   (4) set region settings (all of them)
  //   (5) set group settings that should accumulate
  //   (6) set isntrum settings that should accumulate
  // the order of 5 and 6 can be swapped because accumulation is commutative
}

bool SampleBusPlayer::setGroupOrInstrumToPlay(const SfzInstrument::HierarchyLevel* thingToPlay,
  uchar key, uchar vel, double sampleRate, RegionPlayer* rp, bool busMode, 
  PlayStatus* intermediates)
{
  RAPT::rsAssert(busMode == true);
  // It makes no sense to use a GroupPlayer when not in busMode. Maybe remove the parameter

  //PlayStatus dummy;  
  // I think, using a dummy is wrong an this is what and breaks the test. Yes - this seems to be 
  // the case. I guess, we need to take a PlayStatus parameter that must be owned 
  // somehwere higher up


  if(thingToPlay == grpOrInstr) {
    setupDspSettings(grpOrInstr->getSettings(), sampleRate, rp, busMode, intermediates);
    return true;  }
    // This is not a new group or restart of the whole InstrumPlayer so we may only need to set up
    // those settings that affect the RegionPlayer, i.e. offset, delay, inc, etc. The other 
    // settings are actually already all set up. Maybe split out a setupPlayerSettings such that we
    // can call only that and don't need to loop through all the settings that don't change...but 
    // maybe that's too complicated to do because the settings are not ordered by type.

  // A GroupPlayer needs to play back a new group or the InstrumPlayer was triggered anew. In such 
  // a case, we need to assemble the DSP chain first:
  disassembleProcessors();
  grpOrInstr = thingToPlay;
  if(grpOrInstr != nullptr) {
    if(!assembleProcessors(busMode)) {
      grpOrInstr = nullptr;
      return false;   }
    setupDspSettings(grpOrInstr->getSettings(), sampleRate, rp, busMode, intermediates);
    effectChain.prepareToPlay(key, vel, sampleRate); 
    rp->setupFromIntemediates(*intermediates, sampleRate); // We need to do this again
  }
  return true;
}

bool SampleBusPlayer::assembleProcessors(bool busMode)
{
  RAPT::rsAssert(busMode == true);
  // If we are not in busMode, this function should actually not even get called because only in
  // busMode, the Group- or InstrumentPlayer's own DSP chain is used. We need to take the busMode
  // parameter anyway because this function is an override.

  return SamplePlayer::assembleProcessors(
    grpOrInstr->getOpcodeTypeChain(), grpOrInstr->getModulationSettings());
    // We need only to take into account the group's DSP settings. The instrument's DSP settings
    // can safely be ignored if we are in busMode (which is supposed to be always the case) because 
    // in busMode, the InstrumentPlayer will take care of the instrument's DSP settings
}

//=================================================================================================
// GroupPlayer

void GroupPlayer::processFrame(float* L, float* R)
{
  *L = *R = 0.f;
  for(size_t i = 0; i < regionPlayers.size(); i++) {
    float tmpL, tmpR;
    regionPlayers[i]->processFrame(&tmpL, &tmpR);
    *L += tmpL; 
    *R += tmpR; }
  effectChain.processFrame(L, R);
}

void GroupPlayer::releaseResources()
{
  SampleBusPlayer::releaseResources();
  regionPlayers.clear();
}

void GroupPlayer::addRegionPlayer(RegionPlayer* newPlayer)
{
  RAPT::rsAssert(!RAPT::rsContains(regionPlayers, newPlayer));
  regionPlayers.push_back(newPlayer);
}

void GroupPlayer::removeRegionPlayer(RegionPlayer* player)
{
  RAPT::rsAssert(RAPT::rsContains(regionPlayers, player)); // ToDo: add and use rsContainsOnce
  RAPT::rsRemoveFirstOccurrence(regionPlayers, player);
}

void InstrumPlayer::addRegionPlayer(RegionPlayer* newPlayer)
{
  //RAPT::rsError("not yet implemented");
  // we may nee to accumulate into the region players delay, pitch etc. variables. a partial
  // DSP setup only for those opcodes thataffect the source

}


}}