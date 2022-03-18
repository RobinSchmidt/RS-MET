namespace rosic {
namespace Sampler {

//=================================================================================================
// Move into SamplerEffects.cpp, maybe give declarations in .h:

Processor* getProcessor(std::vector<Processor*>& processors, OpcodeType type, int index)
{
  RAPT::rsAssert(index >= 1 || index == -1);
  index = RAPT::rsMax(index-1, 0);
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int)processors.size(); i++) {
    Processor* dsp = processors[i];
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++;
    }
  }
  return nullptr;
}

/** Returns the total number of processors in the given array. */
size_t getNumProcessorsOfType(const std::vector<Processor*>& processors, OpcodeType type)
{
  size_t count = 0;
  for(size_t i = 0; i < processors.size(); i++) {
    if(processors[i]->getType() == type)
      count++;
  }
  return count;
}

int findProcessorIndex(const std::vector<Processor*>& processors, OpcodeType type, int index)
{
  // ToDo: 
  // -Check, if index < 0 and if so, modify it to processors.size() + abs(index) to use indices 
  //  -1,-2,-3 for the hardwired modulators for amp, cutoff, pitch (maybe -4 for cutoff2)

  int count = 0;
  for(size_t i = 0; i < processors.size(); ++i) {
    if(processors[i]->getType() == type) {
      count++;
      if(count == index)
        return (int)i; }}
  return -1;
}

Processor* findProcessor(const std::vector<Processor*>& processors, OpcodeType type, int index)
{
  int i = findProcessorIndex(processors, type, index);
  if(i == -1)
    return nullptr;
  else
    return processors[i]; 
}

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
// move into rsArrayTools...maybe we need also a version that takes a std::vector


//=================================================================================================
// SamplePlayer

bool SamplePlayer::augmentOrCleanProcessors(const std::vector<OpcodeType>& dspTypeChain)
{
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
      if(getNumProcessorsOfType(effectChain, opType) >= sfzIndex)
        continue;

      // OK - now we actually need to grab another effect of given type from the pool:
      Processor* p = getEffect(opType);  // use a general getProcessor function
      if(p)
      {
        p->setParametersToDefaults(sfzIndex);
        effectChain.push_back(p);
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
      if(getNumProcessorsOfType(modSources, opType) >= sfzIndex)
        continue;
      Processor* p = getModulator(opType);  // use a general getProcessor function
      if(p) {
        p->setParametersToDefaults(sfzIndex);
        modSources.push_back(p);     }
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
  // Maybe remove the calls to setParametersToDefaults - instead we should perhaps reset them to
  // defaults when the object is reposited - but: when the default value depends on the index as
  // in the eqN_freq opcode, we don't really know the index when the object is just sitting in the
  // pool, so maybe it's indeed more appropriate to do it here. Maybe do both? The rationale is to
  // have all objects in the pool in a well defined state. That's not actually important though but
  // it somehow seems cleaner.
}

bool SamplePlayer::assembleModulations(const std::vector<ModulationSetting>& modSettings)
{
  RAPT::rsAssert(dspPool);

  //RAPT::rsAssert(modMatrix.empty(), "Someone has not cleaned up the modMatrix");
  // no - we don't expect it to be empty here anymore because the function may get called 3
  // times during assembly (with the region, group, global modSettings repectively).
  // But maybe the caller shoould have such an assert

  for(size_t i = 0; i < modSettings.size(); i++)
  {
    ModulationSetting ms = modSettings[i];

    // Determine pointer to modulation source Processor:
    int j = findProcessorIndex(modSources, ms.getSourceType(), ms.getSourceIndex());
    if(j == -1)
      continue;
    Processor* srcProc = modSources[j];
    // The j == -1 case may occur when a mod-routing exsists from some source like an LFO to some 
    // parameter but the source itself does not exist because there is no frequency setting for the 
    // LFO, example: lfo3_cutoff2 exists either on global, group or region level but on neither of 
    // these levels does an lfo3_freq setting exist (or any other setting that would us cause to 
    // put lfo3 into the modSources). In such a case, we just skip this loop iteration.

    // Determine pointers to modulation target Processor and Parameter and append them to the 
    // respective arrays, if they are not already present there:
    Processor* tgtProc;
    if(SfzCodeBook::isModSourceSetting(ms.getTargetType())) 
      tgtProc = findProcessor(modSources,  ms.getTargetType(), ms.getTargetIndex()); // Receiver is another modulator
    else
      tgtProc = findProcessor(effectChain, ms.getTargetType(), ms.getTargetIndex()); // Receiver is an effect
    RAPT::rsAssert(tgtProc);
    Parameter* param = tgtProc->getParameter(ms.getTargetOpcode());
    RAPT::rsAssert(param);
    RAPT::rsAppendIfNotAlreadyThere(modTargetProcessors, tgtProc);
    RAPT::rsAppendIfNotAlreadyThere(modTargetParams,     param); 

    // Figure out, if a suitable connection already exists in out modMatrix, if so, update its 
    // depth and mode. Such a thing happens when for example a group and region setting exists for
    // a given connection and the region setting overrides the group setting:
    bool conUpdated = false;  // flag that we set when we updated a connection
    for(size_t k = 0; k < modMatrix.size(); k++) {
      ModulationConnector* mc = modMatrix[k];
      if(mc->getSourceProcessor() == srcProc && mc->getTargetParam() == param) {
        mc->setDepth(ms.getDepth());
        mc->setMode( ms.getMode());
        conUpdated = true;
        break;      }}     // Done! Leave loop early.

    // If we could not update an existing connection, we must actually grab a new one from the pool
    // and set it up from scratch:
    if(!conUpdated) {
      ModulationConnector* mc = dspPool->grabConnector();
      if(mc == nullptr)
        return false;      // Not enough connectors available in pool
      mc->setSource(srcProc);
      mc->setSourceIndex(j);
      mc->setTarget(tgtProc, param);
      mc->setDepth(ms.getDepth());
      mc->setMode( ms.getMode());
      modMatrix.push_back(mc); }

    int dummy = 0;
  }

  return true;


  /*
  // Check, if enough connectors are available:
  if(dspPool->getNumIdleConnectors() < (int)modSettings.size())
    return false;

  // Do the assembly:
  for(size_t i = 0; i < modSettings.size(); i++)
  {
    // todo: 
    // -figure out, if the mod-matrix already contains a connection between the source and target
    //  and if so, juts update its depth and mod parameters instead of adding a new connection
    //  -maybe reorganize the loop such that we call grabConnector only when needed, i.e. when we
    //   can't just update an existing connection

    ModulationSetting    ms = modSettings[i];
    ModulationConnector* mc = dspPool->grabConnector();
    RAPT::rsAssert(mc); // we already verified that enough are available

    // Determine pointer to modulation source and its index in our modSources array and set it up 
    // in the connector:
    int j = findProcessorIndex(modSources, ms.getSourceType(), ms.getSourceIndex());
    RAPT::rsAssert(j >= 0);
    Processor* src = modSources[j];
    RAPT::rsAssert(src);
    mc->setSourceIndex(j);
    mc->setSource(src);

    // Determine pointer to modulation target (Processor and Parameter) and set it up in the
    // connector:
    Processor* prc;
    if(SfzCodeBook::isModSourceSetting(ms.getTargetType())) 
      prc = findProcessor(modSources,  ms.getTargetType(), ms.getTargetIndex()); // Receiver is another modulator
    else
      prc = findProcessor(effectChain, ms.getTargetType(), ms.getTargetIndex()); // Receiver is an effect
    RAPT::rsAssert(prc);
    Parameter* param = prc->getParameter(ms.getTargetOpcode());
    RAPT::rsAssert(param);
    mc->setTarget(prc, param);
    RAPT::rsAppendIfNotAlreadyThere(modTargetProcessors, prc);
    RAPT::rsAppendIfNotAlreadyThere(modTargetParams,     param);  // new
    //RAPT::rsAppend(modTargetParams, param); // old
    
    // Set up modulation depth and mode and add the connection to the modMatrix:
    mc->setDepth(ms.getDepth());
    mc->setMode(ms.getMode());
    modMatrix.push_back(mc);
    // here, we need to figure out, if such a connection already exists and if so, just overwrite 
    // its depth value instead of adding a new connection
  }

  return true;
  */


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
  // i think, we cannot call this here - we must call it after *all* 3 calls of 
  // augmentOrCleanProcessors have been run through - otherwise we ma encounter a situation where
  // we try to connect a parameter with a not yet existing group or global modulator that only 
  // comes into existence after the group and global modulators have been added

  return true;
}

void SamplePlayer::disassembleProcessors()
{
  for(size_t i = 0; i < effectChain.size(); i++)
    dspPool->repositEffect(effectChain[i]);
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
  Processor* dsp = getProcessor(effectChain, s.getTargetOpcodeType(), s.getIndex());
  if(dsp != nullptr)
    dsp->setParameter(s.getOpcode(), s.getValue());
  else
    RAPT::rsError("No processor available for DSP opcode");
    // We could not find a suitable processor in our dspChain to which the given setting could be
    // applied. If this happens, something went wrong (i.e. we have a bug) in assembleDspChain or 
    // dspChain.getProcessor.
}

void SamplePlayer::setupModSourceSetting(const PlaybackSetting& s)
{
  Processor* dsp = getProcessor(modSources, s.getTargetOpcodeType(), s.getIndex());
  if(dsp != nullptr)
    dsp->setParameter(s.getOpcode(), s.getValue());
  else
    RAPT::rsError("No processor available for DSP opcode");
  // This function almost duplicates setupProcessorSetting -> refactor to get rid of the 
  // duplication. Maybe setupProcessorSetting should receive a pointer to effectChain.processors 
  // or to modSources.
}

void SamplePlayer::setupModRoutingSetting(const PlaybackSetting& s)
{
  RAPT::rsError("Not yet implemented");
}

void SamplePlayer::setupDspSettings(const std::vector<PlaybackSetting>& settings,
  RegionPlayer* rp, bool busMode)
{
  SfzCodeBook* cb = SfzCodeBook::getInstance();
  for(size_t i = 0; i < settings.size(); i++)
  {
    PlaybackSetting s = settings[i];
    Opcode op = s.getOpcode();
    if(     cb->isEffectSetting(op))     { setupProcessorSetting( s    ); }
    else if(cb->isModSourceSetting(op))  { setupModSourceSetting( s    ); }
    else if(cb->isPlayerSetting(op))     { setupPlayerSetting(    s, rp); }

    //else if(cb->isModRoutingSetting(op)) { setupModRoutingSetting(s);   }
    // This branch seems to never get hit...maybe check in the caller(s)...maybe *it* needs to call
    // a different function because the modRoutings settings are not held in the same settings 
    // array as all the other settings...yes, i think, that's the case..

  }
  // Try to refactor stuff in a way that lets use get rid of the branching and treat all cases
  // uniformly...I'm not yet sure, if that's possible in any meaningful way, though. Maybe first
  // try to avoid to pass all the additional parameters to setupPlayerSetting or if it really needs
  // them, pass them to setupProcessorSetting, too even it it doesn't make use of them
  // Why has setupPlayerSetting a different signature, i.e. takes the additional rp parameter? Can
  // we get rid of this?

  // I think , we need also a branch setupModRouting
}

//=================================================================================================
// RegionPlayer

rsReturnCode RegionPlayer::setRegionToPlay(const Region* regionToPlay,
  const AudioFileStream<float>* sampleStream, uchar key, uchar vel, bool busMode)
{
  releaseResources(); // actually, it should not hold any at this point - or should it?
  region = regionToPlay;
  if(region == nullptr)
    return rsReturnCode::nothingToDo;
  stream = sampleStream;
  return prepareToPlay(key, vel, busMode);
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

  //double sampleRate = 44100.0; // preliminary - should be inquired from PlayStatus


  int numSources = (int) modSources.size();

  //std::vector<float> modBuffer;
  //modBuffer.resize(2*numSources);
  // old

  // preliminary - should become a member of PlayStatus, the size should be pre-allocated according
  // to the maximum number of modulators that a region/Group/Instrument has

  std::vector<float>& modBuffer = playStatus->modBuffer;
  // new

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

    float m  = modBuffer[2*i];                  // modulator output
    // Preliminary - we currently only use the 1st channel output maybe use something like 
    // rsVector2D<float> for m and do some appropriate casting.

    float c  = con->getContribution(m, u);      // compute modulation contribution
    par->applyModulation(c);                    // accumulate the contribution


    con->getTargetProcessor()->setDirty();
    // ToDo: dirtification should be more sophisticated: Detect whether the modulatedValue is 
    // actually different than st the previous sample. Maybe collect all modulatedValues in a 
    // buffer before modulation and then in the loop over the modTargetProcessors, call 
    // updateCoeffs only if needed...maybe par->applyModulation should return a bool to 
    // indicate, if the value did actually change ...maybe something like:
    //   if(par->applyModulation(c))
    //     con->getTargetProcessor()->setDirty(true);

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
    if(modTargetProcessors[i]->isDirty())
      modTargetProcessors[i]->updateCoeffs(playStatus->sampleRate);
    // check, if there could be duplicates in modTargetProcessors - if so, avoid that in the 
    // assembly of the modulations
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
  processFrame1(effectChain, L, R);


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
  // later. Actually , this should probably be done in SamplerEngine::preAllocateDspMemory
}



rsReturnCode RegionPlayer::prepareToPlay(uchar key, uchar vel, bool busMode)
{
  RAPT::rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  RAPT::rsAssert(stream != nullptr);   // Ditto.

  this->key = key;

  if(!assembleProcessors(busMode)) {
    releaseResources();
    return rsReturnCode::layerOverload; }

  resetPlayerSettings();
  setupDspSettingsFor(region, busMode);

  double fs = playStatus->sampleRate;  // todo: use float
  prepareToPlay1(modSources, key, vel, fs);
  prepareToPlay1(effectChain, key, vel, fs); 
  // The rationale for preparing the modSources first is that their initial output may already 
  // affect the initial parameters of the effects(?) ...but does that matter? Aren't the params 
  // recomputed in processFrame anyway?...perhaps it doesn't matter, but the order feels right 
  // this way anyway. This prepareToPlay function may be a member of SamplePlayer.

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
  if(!areProcessorsEmpty()) {    // Sanity check
    RAPT::rsError("Someone has not cleaned up after finishing playback!");
    disassembleProcessors();  }  // ...so we do it here. But this should be fixed elsewhere!

  const Region* reg = region;
  const Group*  grp = reg->getGroup();
  const Global* ins = grp->getInstrument();

  // old:
  /*
  // The DSPs for which the region itself defines settings/opcodes are always needed:
  bool ok = SamplePlayer::assembleProcessors(
    reg->getOpcodeTypeChain(), reg->getModulationSettings());
  if(!ok)
    return false;

  // If we are not in busMode, the enclosing group and/or enclosing instrument settings act as
  // fallback values for the region so we may require additional DSPs to apply these opcodes
  // to the region, too:
  if(!busMode)  // maybe rename to mixMode
  {
    ok = augmentOrCleanProcessors(grp->getOpcodeTypeChain());    if(!ok) return false;
    ok = augmentOrCleanProcessors(grp->getOpcodeTypeChain());    if(!ok) return false;



    ok = assembleModulations(     grp->getModulationSettings()); if(!ok) return false;
    ok = assembleModulations(     glb->getModulationSettings()); if(!ok) return false;
    // Maybe factor out into a function assembleFallbackProcessors()
  }
  // I think, here we need to do more to set up the modulation connection settings. Not the 
  // modulator parameter though. Things like lfoN_freq do already work but things like lfoN_cutoff
  // don't. The latter one is a modulation routing or conection setting which has a different
  // quality than the modulator's internal parameters....
  // But this is wrong: when doing it in this order, the group and global modulation settings will
  // override the region settings and not the other way around. we need to augment also
  // setupDspSettingsFor
  */

  // new:
  bool ok = true;

  // Assemble processors (modulators and effects):
  ok &= augmentOrCleanProcessors(reg->getOpcodeTypeChain());
  if(!busMode) { // maybe rename to mixMode
    ok &= augmentOrCleanProcessors(grp->getOpcodeTypeChain());
    ok &= augmentOrCleanProcessors(ins->getOpcodeTypeChain());  }
  if(!ok) {
    this->disassembleProcessors();
    return false; }
  // I think, here, the order of the 3 calls doesn't matter

  // Assemble modulation connections:
  if(!busMode) {
    ok &= assembleModulations(ins->getModulationSettings());
    ok &= assembleModulations(grp->getModulationSettings()); }
  ok &= assembleModulations(reg->getModulationSettings());
  if(!ok) {
    this->disassembleProcessors();  // function also disassembles the modMatrix
    return false; }
  // ..but here, the order of the 3 calls does matter...i think
  // ...but maybe we should make it consistent anyway - always call in the order: ins, grp, reg









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

void RegionPlayer::setupDspSettingsFor(const Region* r, bool busMode)
{
  // To set up the settings, we call setupDspSettings 3 times to:
  //   (1) set up the general instrument-wide settings
  //   (2) set up group specific settings (this may override instrument settings)
  //   (3) set up region specific settings (this may override group and/or instrument settings)
  // but only if not in bus-mode in which case the group and instrument settings apply to separate
  // DSP objects on the groups which map to sub-busses and the instrument which maps to the final
  // master bus:
  if(!busMode) {
    setupDspSettings(region->getGroup()->getInstrument()->getSettings(), this, busMode);  
    setupDspSettings(region->getGroup()->getSettings(), this, busMode); }
  setupDspSettings(region->getSettings(), this, busMode);
  // Works for accumulative setting but not for those that should awlays override



  // Compute the final member variables from the intemediates:
  setupFromIntemediates();
  // maybe call this only when not in busMode because in busMode, it will get called again later 
  // and we don't want to call it twice. It will not do anything wrong, but it does the computation
  // at the first time for nothing....
}

void RegionPlayer::setupFromIntemediates() // fix typo!..actually this function needs a better name anyway
{
  // This may get called twice on noteOn in busMode -> verify and try to avoid the first call

  PlayStatus& iv = *playStatus;  // maybe rename
  double fs = iv.sampleRate; 

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

void RegionPlayer::setupPlayerSetting(const PlaybackSetting& s, RegionPlayer* rp)
{
  RAPT::rsAssert(rp == this);
  // For RegionPlayer objects like this, this function is supposed to be called only for the object
  // itself. The rp pointer is needed here only to conform to the baseclass interface. For 
  // GroupPlayers and InstrumPlayers, the pointer is supposed to hold the RegionPlayer to which 
  // this setting should be applied accumulatively in busMode. In default mode, GroupPlayer and 
  // InstrumPlayer play no role at all.

  PlayStatus* iv = playStatus; // rename
  double sampleRate = iv->sampleRate;



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

void SampleBusPlayer::setupPlayerSetting(const PlaybackSetting& s, RegionPlayer* rp)
{
  // We are supposedly a higher level player object like GroupPlayer or InstrumentPlayer but the 
  // setting in question reaches through to an underlying (embedded) RegionPlayer object and must 
  // be accumulated into the respective value there. This applies to all kinds of settings that 
  // cannot be meaningfully achieved by post-processing in the effect chain but instead must be 
  // applied directly at the playback source. Some of them (like delay) *could* be achieved 
  // also by post-processing with a suitable effect but it's more efficient to do it at the source.
  // Others like all tuning related stuff indeed need to be done at the source.

  PlayStatus* iv = playStatus; 
  double fs = playStatus->sampleRate;


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
  uchar key, uchar vel, RegionPlayer* rp, bool busMode)
{
  RAPT::rsAssert(busMode == true);
  // It makes no sense to use a GroupPlayer when not in busMode. Maybe remove the parameter

  //PlayStatus dummy;  
  // I think, using a dummy is wrong an this is what and breaks the test. Yes - this seems to be 
  // the case. I guess, we need to take a PlayStatus parameter that must be owned 
  // somehwere higher up


  if(thingToPlay == grpOrInstr) {
    setupDspSettings(grpOrInstr->getSettings(), rp, busMode);
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
    setupDspSettings(grpOrInstr->getSettings(), rp, busMode);
    prepareToPlay1(effectChain, key, vel, playStatus->sampleRate); 
    rp->setupFromIntemediates(); // We need to do this again
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
  processFrame1(effectChain, L, R);
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

/*



*/