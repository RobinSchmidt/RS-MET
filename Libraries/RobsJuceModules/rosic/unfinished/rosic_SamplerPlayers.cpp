namespace rosic {
namespace Sampler {

//=================================================================================================
// rsSamplerEngine::SignalProcessorChain

void SignalProcessorChain::processFrame(rsFloat64x2& inOut)
{
  for(size_t i = 0; i < processors.size(); i++)
    processors[i]->processFrame(inOut);
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

SignalProcessor* SignalProcessorChain::getProcessor(
  DspType type, int index)
{
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int)processors.size(); i++) {
    SignalProcessor* dsp = getProcessor(i);
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++;
    }
  }
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

bool SamplePlayer::addDspsIfNeeded(const std::vector<DspType>& dspTypeChain)
{
  for(int i = 0; i < (int)dspTypeChain.size(); i++)
  {
    // Figure out the type of the DSP that may need to be added to the chain:
    DspType dspType = dspTypeChain[i];

    // Figure out the index within the type-chain (i.e. placeholder chain for the actual DSPs) of 
    // the potentially to-be-added DSP because whether or not we will actually need to add 
    // another DSP to the chain will depend on that index and the number of alike DSPs already 
    // present in the chain. Furthermore, the default value may depend on the index as well:
    int index = 1;
    for(int j = i-1; j >= 0; j--) {
      if(dspTypeChain[j] == dspType)
        index++;
    }
// Maybe factor that out into some dspTypeChain.getSfzIndex(arrayIndex) method. In this call,
// the arrayIndex should be the array-index of a DSP of given type within the chain and the 
// return value should be the sfz-index. Example: 
//   typeChain == { Filter, Equalizer, Filter, Filter, Equalizer, Equalizer, Filter }
//   array-index:     0         1        2       3         4          5        6
//   sfz-index:       1         1        2       3         2          3        4
// Then, when we pass 4 (the array index), it should return 2 because array index 4 refers
// to the 2nd equalizer in the desired sfz-effect chain, i.e. the equalizer that is 
// controlled with opcodes eq2_gain, eq2_freq, eq2_bw. This kind of index-mapping is a bit 
// confusing and needs good documentation and unit-tests. But TypeChain is actually not a 
// class but just a std::vector of DspType. Maybe make it a class..hmm...not sure if that's
// worth it. Maybe make a free helper function getSfzDspIndex(const TypeChain& c, int i).
// Could be a static method of SamplerData. Maybe use size_t: init index to 0 and start loop
// index j at i

// Figure out, if we actually need to add another DSP to the chain. If not, there's nothing
// more to do in this iteration:
    if(dspChain.getNumProcessors(dspType) >= index)
      continue;

    // OK - now we actually need to grab another DSP of given type from the pool:
    SignalProcessor* dsp = getProcessor(dspType);
    if(dsp)
    {
      dsp->resetSettings(index);
      dspChain.addProcessor(dsp);
    }
    else {
      return false;
      // Not enough DSPs of desired type are available in the pool so we report failure. In such a
      // case, it is the job of the caller to roll back any partially built chain, if needed.
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
    disassembleDspChain();
  }  // ...so we do it here. But this should be fixed elsewhere!
  if(!addDspsIfNeeded(dspTypes)) {
    disassembleDspChain();    // addDspsIfNeeded may have built a partial chain which we then..
    return false;
  }           // ..need to clean up here
  return true;

  // The fact that addDspsIfNeeded may build a partial chain without cleaning it up itself if
  // it can't build the full chain is a bit messy. But just putting the clean-up responsibility
  // into this function isn't ideal either because then the name doesn't really fit. Maybe rename
  // it addDspsIfNeededCleanUpWhenFail...but that's a monster of a name...we'll see.... maybe
  // augmentOrCleanDspChain
}

void SamplePlayer::disassembleDspChain()
{
  for(int i = 0; i < dspChain.getNumProcessors(); i++)
    dspPool->processorPool.repositProcessor(dspChain.getProcessor(i));
  dspChain.clear();
}

void SamplePlayer::setupProcessorSetting(const PlaybackSetting& s)
{
  using SD = rsSamplerData::PlaybackSetting;
  // ToDo: We need to call the static member function getTargetProcessorType of that class. That
  // function should be moved elsewhere. We want a sort of database to retrieve all sort of info
  // about opcodes, including to what type of processor they apply

  // Internal helper function to retrieve a pointer to the proccessor within our dspChain to which 
  // the setting applies. It may at some point be dragged out of this function if it turns out to 
  // be useful in other places as well:
  auto getProcessorFor = [this](const PlaybackSetting& s)
  {
    DspType dspType = SD::getTargetProcessorType(s.getType());

    int i = RAPT::rsMax(s.getIndex(), 0);
    // This is still wrong. It works only when s.getIndex() returns -1 or 0. -1 is the default 
    // encoding "not applicable".

    // maybe introduce a getMappedIndex() function that does the right thing, rename index to 
    // getSfzIndex to make it clear that this is the index that occurs in the sfz-file. The other
    // one could then be named getChainIndex or something. ok - for the time beign, let's do it 
    // directly here:
    i = s.getIndex();
    if(i >   0) i--;    // map from 1-based counting to 0-based
    if(i == -1) i = 0;  // map code for "no index" to 1st
    SignalProcessor* dsp = dspChain.getProcessor(dspType, i);
    return dsp;
  };
  // Maybe wrap this whole business into a member function of the dspChain. This class knows best
  // how to map sfz-indices to the actual processor within the chain


  SignalProcessor* dsp = getProcessorFor(s);
  if(dsp != nullptr)
    dsp->setParameter(s.getType(), s.getValue());
  else
    RAPT::rsError("No processor available for DSP opcode");
  // We could not find a suitable processor in our dspChain to which the given setting could be
  // applied. If this happens, something went wrong (i.e. we have a bug) in buildDspChain or 
  // getProcessorFor.
}

void SamplePlayer::setupDspSettings(const std::vector<PlaybackSetting>& settings,
  double sampleRate, bool busMode)
{


  RAPT::rsError("not yet implemented");
}



//=================================================================================================
// RegionPlayer

rsReturnCode RegionPlayer::setRegionToPlay(const Region* regionToPlay,
  const AudioFileStream<float>* sampleStream, double fs, bool busMode)
{
  releaseResources(); // actually, it should not hold any at this point - or should it?
  region = regionToPlay;
  if(region == nullptr)
    return rsReturnCode::nothingToDo;
  stream = sampleStream;
  return prepareToPlay(fs, busMode);
}

rsFloat64x2 RegionPlayer::getFrame()
{
  float L, R;                        // left and right output

  // Negatively initialized sampleTime implements delay. If we are still "waiting", we just 
  // increment the time and return 0,0. Actual output will be produced as soon as sampleTime 
  // reaches zero. 
  if(sampleTime < 0.0)
  {
    sampleTime += 1.0;
    if(sampleTime >= 0.0)
    {
      sampleTime += offset;   // old
      //sampleTime += offset - 1; // new

      //stream->getFrameStereo((float)sampleTime, &L, &R);
      //return this->amp * rsFloat64x2(L, R); 
    }
    return rsFloat64x2(0.0, 0.0);
  }
  stream->getFrameStereo((float)sampleTime, &L, &R);  // try to avoid the conversion to float
  // -implement better interpolation methods (sinc, elephant, ...)
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
  rsFloat64x2 out(L, R);
  dspChain.processFrame(out);
  return this->amp * out;


  //return this->amp * rsFloat64x2(L, R); 
}

void RegionPlayer::processBlock(rsFloat64x2* y, int N)
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

rsReturnCode RegionPlayer::prepareToPlay(double fs, bool busMode)
{
  RAPT::rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  RAPT::rsAssert(stream != nullptr);   // Ditto.

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
  setupDspSettingsFor(region, fs, busMode);
  // todo: move fs before the override parameters for consistency

  // todo: setup modulators and modulation connections

  dspChain.prepareToPlay(fs);
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

  // todo:
  // -check also, if the amplitude envelope has reached its end
  // -hmm - maybe, if we allow the frequency envelope to go negative, we could also move 
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
    if(!addDspsIfNeeded(region->getGroup()->getProcessingChain())) {
      disassembleDspChain();
      return false;
    }
    if(!addDspsIfNeeded(region->getGroup()->getInstrument()->getProcessingChain())) {
      disassembleDspChain();
      return false;
    }
  }

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
  amp        = 1.0;
  sampleTime = 0.0;
  increment  = 1.0;
  offset     = 0.f;

  endTime    = (float)stream->getNumFrames();
  // Maybe use -1? That may require updating the unit tests. But maybe it's appropriate to use 
  // numFrames when assuming linear interpolation. I think, for general interpolators, we should 
  // use endTime = numFrames - 1 + kernelSize/2. Test this with very high downshifting factors and
  // maybe with a sample that is just 1 sample long with a value of 1. We should see the 
  // interpolation kernel as output.

  loopStart  = 0.f;
  loopEnd    = 0.f;
  loopMode   = 0;

  //dspChain.resetSettings();
}

void RegionPlayer::setupDspSettingsFor(const Region* r, double fs, bool busMode)
{
  // To set up the settings, we call setupDspSettings 3 times to:
  //   (1) set up the general instrument-wide settings
  //   (2) set up group specific settings (this may override instrument settings)
  //   (3) set up region specific settings (this may override group and/or instrument settings)
  // but only if not in bus-mode in which case the group and instrument settings apply to separate
  // DSP objects on the groups which map to sub-busses and the instrument which maps to the final
  // master bus:
  if(!busMode)
    setupDspSettings(region->getGroup()->getInstrument()->getSettings(), fs, busMode);
  if(!busMode)
    setupDspSettings(region->getGroup()->getSettings(), fs, busMode);
  setupDspSettings(region->getSettings(), fs, busMode);

  // The code above will have set up the increment assuming that the current key matches the 
  // rootKey of the sample, Now, as final step, we also adjust it according to the difference 
  // between the played key and the sample's rootKey. This is not done in the code above because
  // it makes no sense to accumulate this parameter
  increment *= stream->getSampleRate() / fs;
  double rootKey = region->getSettingValue(Opcode::PitchKeyCenter, -1, false);
  double pitchOffset = double(key) - rootKey;
  //increment *= RAPT::rsPitchOffsetToFreqFactor(pitchOffset); 
  increment *= pow(2.0, pitchOffset / 12.0);
  // The formula using rsPitchOffsetToFreqFactor is imprecise: when we have a pitchOffset of 
  // exactly -12, for example, we want the increment be multiplied by exactly 0.5, but using this
  // function, the factor comes out as 0.49999..., so we use the more precise (and more 
  // expensive) call to pow. It's not per-sample here code anyway, so we may afford that.

  // If there is no delay to be considered, we advance the sample time into the sample by the 
  // desired offset. Otherwise, sampleTime will have been initialized to -delaySamples and we leave
  // it at that. In this case, the offset will be handled in getFrame etc.:
  if(sampleTime >= 0.0)
    sampleTime = offset;
}

void RegionPlayer::setupDspSettings(
  const std::vector<PlaybackSetting>& settings, double fs, bool busMode)
{
  //using PS = PlaybackSetting;
  using TP = Opcode;               // rename to OC

  //double  tuneCoarse = 0.0;        // in semitones
  //double  tuneFine   = 0.0;        // in cents

  // Loop through the settings of the region and for each setting that is present, change the 
  // value from its default to the stored value:
  for(size_t i = 0; i < settings.size(); i++)
  {

    PlaybackSetting setting = settings[i];
    TP type = setting.getType();                  // rename to opcode
    //double val = (double)setting.getValue();

    // todo: 
    // -factor out into virtual function setDspSetting(setting) which is a virtual function
    //  defined in SamplePlayer and we override it here to ctach the RegionPlayer-sepcific
    //  opcodes such as delay

    if(setting.isDspSetting())
    {
      setupProcessorSetting(setting);
      continue;
    }
    if(setting.isPlayerSetting())
    {
      setupPlayerSetting(setting, fs);
      continue;
    }
    // optimize: instead of calling isDspSetting in the setting object, use the setting.getType()
    // and a pointer to an SfzCodebook that has been acquired before the loop (saves some checks
    // in the aqcuiry of the singleton pointer)
    

    /*
    // move these into if(setting.isPlayerSetting())
    switch(type)
    {

    // Pitch settings:
    //case TP::PitchKeyCenter: { rootKey    = val; } break;  // done by caller
    case TP::Transpose:      { tuneCoarse = val; } break;
    case TP::Tune:           { tuneFine   = val; } break;

    //
    case TP::Delay:
    {
      if(busMode) sampleTime += -val * fs;
      else        sampleTime  = -val * fs;
    }  break;

    case TP::Offset:
    {
      if(busMode) offset += float(val);
      else        offset  = float(val);
    }  break;


    }
    */
  }

  /*
  double tune   = tuneCoarse + 0.01 * tuneFine;
  double factor = pow(2.0, tune / 12.0);
  if(busMode) this->increment *= factor;
  else        this->increment  = factor;
  */

}

void RegionPlayer::setupProcessorSetting(const PlaybackSetting& s)
{
  SamplePlayer::setupProcessorSetting(s);
  // preliminary: todo: catch those settings that apply only to the sample-source, i.e. do not
  // apply to a DSP in the chain. they need different handling
}

void RegionPlayer::setupPlayerSetting(const PlaybackSetting& s, double sampleRate)
{
  double tuneCoarse = 0.0;        // in semitones
  double tuneFine   = 0.0;        // in cents
  double val = (double)s.getValue();
  using OC   = Opcode;
  switch(s.getType())
  {
  // Pitch settings:
  //case TP::PitchKeyCenter: { rootKey    = val; } break;  // done by caller
  case OC::Transpose: { tuneCoarse = val;               } break;
  case OC::Tune:      { tuneFine   = val;               } break;
  case OC::Delay:     { sampleTime = -val * sampleRate; } break;
  case OC::Offset:    { offset     = float(val);        } break;
  }
  double tune     = tuneCoarse + 0.01 * tuneFine;
  double factor   = pow(2.0, tune / 12.0);
  this->increment = factor;
}


//=================================================================================================
// GroupPlayer

rsFloat64x2 GroupPlayer::getFrame()
{
  rsFloat64x2 out = 0.0;
  for(size_t i = 0; i < regionPlayers.size(); i++)
    out += regionPlayers[i]->getFrame();
  dspChain.processFrame(out);
  return out;
}

void GroupPlayer::releaseResources()
{
  disassembleDspChain();
  regionPlayers.clear();
  //group = nullptr;   // why is this commented? seems to make sense to set it to null here
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

bool GroupPlayer::setGroupToPlay(const rsSamplerData::Group* groupToPlay, bool busMode)
{
  if(groupToPlay == group)
    return true;               // nothing to do
  disassembleDspChain();
  group = groupToPlay;
  if(group != nullptr) {
    if(!assembleDspChain(busMode)) {
      group = nullptr;
      return false;
    }
    setupDspChain();
  }
  return true;
}

bool GroupPlayer::assembleDspChain(bool busMode)
{
  RAPT::rsAssert(busMode == true);
  // If we are not in busMode, this function should actually not even get called because only in
  // busMode, the GroupPlayer's own DSP chain is used. We need to take the busMode parameter 
  // anyway because this function is an override.

  return SamplePlayer::assembleDspChain(group->getProcessingChain());
  // We need only to take into account the group's DSP settings. The instrument's DSP settings
  // can safely be ignored if we are in busMode (which is supposed to be always the case) because 
  // in busMode, the InstrumentPlayer will take care of the instrument's DSP settings
}

void GroupPlayer::setupDspChain()
{
  RAPT::rsError("Not yet implemented");
}

void GroupPlayer::setupPlayerSetting(const PlaybackSetting& s, double sampleRate)
{
  double val = (double)s.getValue();
  using OC   = Opcode;

  //RegionPlayer* rp = getSourcePlayer(); ...the RegionPlayer object to which this setting
  // applies...how can we obtain a pointer to it? Maybe it needs to be passed in a parameter?
  // ...yeah - that's probably what needs to be done and we can do it in the baseclass

  switch(s.getType())
  {
    // Pitch settings:
  case OC::Transpose:
  {
    // accumulate pitchOffsetToFreqFactor(val) into the increment of rp
  } break;

  case OC::Tune:
  { 
    // accumulate pitchOffsetToFreqFactor(0.01 * val) into the increment of rp
  } break;

  case OC::Delay:     
  { 
    //sampleTime = -val * sampleRate; 
    // accumulate -val * sampleRate into sampleTime of rp

  } break;
  case OC::Offset:    
  { 
    //offset = float(val);
    // accumulate val into offset of rp
  } 
  break;
  }


  //RAPT::rsError("Not yet implemented");
}

//=================================================================================================
// InstrumPlayer

bool InstrumPlayer::assembleDspChain(bool busMode)
{
  RAPT::rsAssert(busMode == true); // see comment in GroupPlayer::assembleDspChain
  return SamplePlayer::assembleDspChain(instrum->getProcessingChain());
}

void InstrumPlayer::setupPlayerSetting(const PlaybackSetting& s, double sampleRate)
{
  RAPT::rsError("Not yet implemented");
}














}}