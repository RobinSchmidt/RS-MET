namespace rosic {
namespace Sampler {

//=================================================================================================
// rsSamplerEngine::SignalProcessorChain

void SignalProcessorChain::processFrame(float* L, float* R)
{
  for(size_t i = 0; i < processors.size(); i++)
    processors[i]->processFrame(L, R);
}

void SignalProcessorChain::processBlock(float* L, float* R, int N)
{
  for(int n = 0; n < N; n++)
    processFrame(L, R);
}

size_t SignalProcessorChain::getNumProcessors(DspType type) const
{
  size_t count = 0;
  for(size_t i = 0; i < processors.size(); i++) {
    if(processors[i]->getType() == type)
      count++;
  }
  return count;
}

SignalProcessor* SignalProcessorChain::getProcessor(DspType type, int index)
{
  RAPT::rsAssert(index >= 1 || index == -1);
  index = RAPT::rsMax(index-1, 0);
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int)processors.size(); i++) {
    SignalProcessor* dsp = getProcessor(i);
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++; }}
  return nullptr;
}

/*
void SignalProcessorChain::resetState()
{
for(size_t i = 0; i < processors.size(); i++)
processors[i]->resetState();
}
*/

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

bool SamplePlayer::augmentOrCleanDspChain(const std::vector<DspType>& dspTypeChain)
{
  for(int i = 0; i < (int)dspTypeChain.size(); i++)
  {
    // Figure out the type of the DSP that may need to be added to the chain:
    DspType dspType = dspTypeChain[i];

    // Figure out, if we actually need to add another DSP to the chain. If not, there's nothing
    // more to do in this iteration:
    int sfzIndex = rsCount(&dspTypeChain[0], i, dspType) + 1;
    if(dspChain.getNumProcessors(dspType) >= sfzIndex)
      continue;

    // OK - now we actually need to grab another DSP of given type from the pool:
    SignalProcessor* dsp = getProcessor(dspType);
    if(dsp)
    {
      dsp->resetSettings(sfzIndex);
      dspChain.addProcessor(dsp);
    }
    else {
      disassembleDspChain();
      return false;
      // Not enough DSPs of desired type are available in the pool so we roll back any partially 
      // built chain and report failure. 
    }
  }
  return true;
  // When we arrive here, we have successfully finished the loop which means that we either could
  // add enough DSPs of the desired types to the chain or they were already present before. In 
  // both cases, the dspChain is now in the required state so we can report success.

}

bool SamplePlayer::assembleDspChain(const std::vector<DspType>& dspTypes)
{
  if(!dspChain.isEmpty()) {
    RAPT::rsError("Someone has not cleaned up after finishing playback!");
    disassembleDspChain();  }  // ...so we do it here. But this should be fixed elsewhere!
  if(!augmentOrCleanDspChain(dspTypes)) 
    return false;
  return true;
}

void SamplePlayer::disassembleDspChain()
{
  for(int i = 0; i < dspChain.getNumProcessors(); i++)
    dspPool->processorPool.repositProcessor(dspChain.getProcessor(i));
  dspChain.clear();
}

void SamplePlayer::setupProcessorSetting(const PlaybackSetting& s)
{
  SignalProcessor* dsp = dspChain.getProcessor(s.getTargetDspType(), s.getIndex());
  if(dsp != nullptr)
    dsp->setParameter(s.getType(), s.getValue());
  else
    RAPT::rsError("No processor available for DSP opcode");
    // We could not find a suitable processor in our dspChain to which the given setting could be
    // applied. If this happens, something went wrong (i.e. we have a bug) in assembleDspChain or 
    // dspChain.getProcessor.
}

void SamplePlayer::setupDspSettings(const std::vector<PlaybackSetting>& settings,
  double sampleRate, RegionPlayer* rp, bool busMode, PlayerIntermediates* iv)
{
  SfzCodeBook* codebook = SfzCodeBook::getInstance();
  for(size_t i = 0; i < settings.size(); i++)
  {
    PlaybackSetting s = settings[i];
    Opcode op = s.getType();
    if(     codebook->isDspSetting(op))    { setupProcessorSetting(s);              }
    else if(codebook->isPlayerSetting(op)) { setupPlayerSetting(s, sampleRate, rp, iv); }
  }
}

//=================================================================================================
// RegionPlayer

rsReturnCode RegionPlayer::setRegionToPlay(const Region* regionToPlay,
  const AudioFileStream<float>* sampleStream, uchar key, uchar vel, double fs, bool busMode, 
  PlayerIntermediates* iv)
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

  // Update our sampleTime counter:
  sampleTime += increment;
  if(loopMode == LoopMode::loop_continuous)
  {
    if(sampleTime >= loopEnd)
      sampleTime -= (loopEnd - loopStart);
  }

  // Apply the DSP chain:
  dspChain.processFrame(L, R);


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

  disassembleDspChain();

  // Return the modulators to the pool:
  // ...something to do...
  //modulators.clear();

  // Return the mod-connections to the pool:
  // ...something to do...
  //modMatrix.clear();

  // Invalidate our pointers:
  stream = nullptr;
  region = nullptr;

  //key = 0;  // shouldn't we do this
}

void RegionPlayer::allocateMemory()
{
  modulators.reserve(8);
  modMatrix.reserve(32);
  dspChain.reserve(8);
  // These values are ad-hoc arbitrarily chosen. Maybe give the function some parameters to choose
  // how many of each should be pre-allocated. It should be enough to avoid having to allocate more
  // later
}

rsReturnCode RegionPlayer::prepareToPlay(uchar key, uchar vel, double fs, bool busMode, 
  PlayerIntermediates* iv)
{
  RAPT::rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  RAPT::rsAssert(stream != nullptr);   // Ditto.

  this->key = key;

  if(!assembleDspChain(busMode))
  {
    releaseResources();
    return rsReturnCode::layerOverload;
  }
  if(!setupModulations()) {
    releaseResources();
    return rsReturnCode::layerOverload;
  }
  resetPlayerSettings();

  setupDspSettingsFor(region, fs, busMode, iv);

  // todo: setup modulators and modulation connections

  dspChain.prepareToPlay(key, vel, fs);
  // modulators.prepareToPlay(fs)

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

bool RegionPlayer::assembleDspChain(bool busMode)
{
  // The DSPs for which the region itself defines settings/opcodes are always needed:
  bool ok = SamplePlayer::assembleDspChain(region->getProcessingChain());
  if(!ok)
    return false;

  // If we are not in busMode, the enclosing group and/or enclosing instrument settings act as
  // fallback values for the region so we may require additional DSPs to apply these opcodes
  // to the region, too:
  if(!busMode) {
    if(!augmentOrCleanDspChain(region->getGroup()->getProcessingChain()))
      return false;
    if(!augmentOrCleanDspChain(region->getGroup()->getInstrument()->getProcessingChain()))
      return false; }

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

bool RegionPlayer::setupModulations()
{
  // ToDo:
  // -Build the set of modulators, similar as in buildProcessingChain
  // -Wire up the modulation connections

  // ...something to do...

  return true;
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
  PlayerIntermediates* iv)
{
  // To set up the settings, we call setupDspSettings 3 times to:
  //   (1) set up the general instrument-wide settings
  //   (2) set up group specific settings (this may override instrument settings)
  //   (3) set up region specific settings (this may override group and/or instrument settings)
  // but only if not in bus-mode in which case the group and instrument settings apply to separate
  // DSP objects on the groups which map to sub-busses and the instrument which maps to the final
  // master bus:
  if(!busMode)
    setupDspSettings(region->getGroup()->getInstrument()->getSettings(), fs, this, busMode, iv);
  if(!busMode)
    setupDspSettings(region->getGroup()->getSettings(), fs, this, busMode, iv);
  setupDspSettings(region->getSettings(), fs, this, busMode, iv);

  // Compute the final member variables from the intemediates:
  setupFromIntemediates(*iv, fs);
  // maybe call this only when not in busMode because in busMode, it will get called again later 
  // and we don't want to call it twice. It will not do anything wrong, but it does the computation
  // at the first time for nothing....
}

void RegionPlayer::setupFromIntemediates(const PlayerIntermediates& iv, double fs) // fix typo!
{
  // This may get called twice on noteOn in busMode -> verify and try to avoid the first call

  // Compute the per-sample increment:
  double rootKey = region->getSettingValue(Opcode::PitchKeyCenter, -1, false);
  double pitchOffset = double(key) - rootKey + iv.transpose + 0.01 * iv.tune;
  increment = pow(2.0, pitchOffset / 12.0) * stream->getSampleRate() / fs;
  // The formula using rsPitchOffsetToFreqFactor is too imprecise: when we have a pitchOffset of 
  // exactly -12, for example, we want the increment be multiplied by exactly 0.5, but using this
  // function, the factor comes out as 0.49999..., so we use the more precise (and more 
  // expensive) call to pow. It's not per-sample here code anyway, so we may afford that.

  // Sanity-check the loop boundaries:
  loopStart = RAPT::rsMax(loopStart, 0.0);
  loopEnd   = RAPT::rsClip(loopEnd, loopStart, (double)stream->getNumFrames());

  // ToDo:
  // -Take into account pitch_keytrack, pitch_veltrack...maybe we need to take key,vel params or we
  //  add such fields to PlayerIntermediates...which should be extended to a "MidiStatus" anyway.
  // -Can we avoid the inquiry for the rootKey? This may be a bit expensive. Maybe the RegionPlayer
  //  should have a rootKey member? That may be messy to deal with in busMode. Maybe try to treat
  //  the pitch_keycenter opcode somehow like the transpose opcode. What would actually be a 
  //  meaningful "accumulative" behavior of a setting like that? Maybe accumulate differences to
  //  the default or something? Maybe keep a rootKeyShift member that is by default zero, 
  //  accumulates th differences of the pitch_keycenter opcode with respce to 60 (the default)?
  //  Try this and write unit tests..
}

void RegionPlayer::setupPlayerSetting(const PlaybackSetting& s, double sampleRate, 
  RegionPlayer* rp, PlayerIntermediates* iv)
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
  switch(s.getType())
  {
  // Player:
  //case TP::PitchKeyCenter: { rootKey    = val; } break;  // done by caller
  case OC::Transpose:     { iv->transpose =  val;                } break;
  case OC::Tune:          { iv->tune      =  val;                } break;
  case OC::Delay:         { sampleTime    = -val * sampleRate;   } break;
  case OC::Offset:        { offset        =  val;                } break;
  case OC::LoopMode:      { loopMode      = (LoopMode)(int) val; } break;  
  case OC::LoopStart:     { loopStart     =  val;                } break;
  case OC::LoopEnd:       { loopEnd       =  val;                } break;

  // Tracking:
  case OC::ampN_veltrack: { 
    iv->ampN_veltrack[N] = val;          } break;
    // ToDo: make sure that ampN_veltrack has large enough size!!!
  }
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
  RegionPlayer* rp, PlayerIntermediates* iv)
{
  // We are supposedly a higher level player object like GroupPlayer or InstrumentPlayer but the 
  // setting in question reaches through to an underlying (embedded) RegionPlayer object and must 
  // be accumulated into the respective value there. This applies to all kinds of settings that 
  // cannot be meaningfully achieved by post-processing in the DSP chain but instead must be 
  // applied directly at the playback source. Some of them (like delay) *could* be achieved 
  // also by post-processing with a suitable DSP but it's more efficient to do it at the source.
  // Others like all tuning related stuff indeed needs to be done at the source.

  RAPT::rsAssert(rp != nullptr);
  double val = (double)s.getValue();
  int    N   = s.getIndex();
  using OC   = Opcode;
  switch(s.getType())
  {
  // Player:
  case OC::Transpose: { iv->transpose  += val;        } break;
  case OC::Tune:      { iv->tune       += val;        } break;
  case OC::Delay:     { rp->sampleTime += -val * fs;  } break;
  case OC::Offset:    { rp->offset     += float(val); } break;

    // Tracking:
  case OC::ampN_veltrack: { 
    iv->ampN_veltrack[N] += (float)val; } break;
    // ToDo: make sure that ampN_veltrack has large enough size!!!
  }

  // Maybe the default branch should call rp->setupPlayerSetting(s, sampleRate, val). That would 
  // revert the behavior to override mode which may make most sense for most settings - certainly
  // for stuff like loop_mode, loop_start, etc. - there is no meaningful way for this setting to
  // accumulate
}

bool SampleBusPlayer::setGroupOrInstrumToPlay(const rsSamplerData::OrganizationLevel* thingToPlay,
  uchar key, uchar vel, double sampleRate, RegionPlayer* rp, bool busMode, 
  PlayerIntermediates* intermediates)
{
  RAPT::rsAssert(busMode == true);
  // It makes no sense to use a GroupPlayer when not in busMode. Maybe remove the parameter

  //PlayerIntermediates dummy;  
  // I think, using a dummy is wrong an this is what and breaks the test. Yes - this seems to be 
  // the case. I guess, we need to take a PlayerIntermediates parameter that must be owned 
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
  disassembleDspChain();
  grpOrInstr = thingToPlay;
  if(grpOrInstr != nullptr) {
    if(!assembleDspChain(busMode)) {
      grpOrInstr = nullptr;
      return false;   }
    setupDspSettings(grpOrInstr->getSettings(), sampleRate, rp, busMode, intermediates);
    dspChain.prepareToPlay(key, vel, sampleRate); 
    rp->setupFromIntemediates(*intermediates, sampleRate); // We need to do this again
  }
  return true;
}

bool SampleBusPlayer::assembleDspChain(bool busMode)
{
  RAPT::rsAssert(busMode == true);
  // If we are not in busMode, this function should actually not even get called because only in
  // busMode, the Group- or InstrumentPlayer's own DSP chain is used. We need to take the busMode
  // parameter anyway because this function is an override.

  return SamplePlayer::assembleDspChain(grpOrInstr->getProcessingChain());
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
  dspChain.processFrame(L, R);
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