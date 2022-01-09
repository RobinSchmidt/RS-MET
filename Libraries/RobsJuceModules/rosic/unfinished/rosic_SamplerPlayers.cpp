namespace rosic { namespace Sampler {

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
      count++; }
  return count;
}

SignalProcessor* SignalProcessorChain::getProcessor(
  DspType type, int index)
{
  int count = 0;  // counts, how many DSPs of given type we have iterated over - why not size_t?
  for(int i = 0; i < (int) processors.size(); i++) {
    SignalProcessor* dsp = getProcessor(i);
    if(dsp->getType() == type) {
      if(count == index)
        return dsp;
      else
        count++;   }}
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
// RegionPlayer

rsReturnCode RegionPlayer::setRegionToPlay(const Region* regionToPlay, 
  double fs, bool groupSettingsOverride, bool regionSettingsOverride)
{  
  releaseDspObjects();
  region = regionToPlay;
  if(region == nullptr)
    return rsReturnCode::nothingToDo;
  stream = getSampleStreamFor(region);
  return prepareToPlay(fs, groupSettingsOverride, regionSettingsOverride);
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

void RegionPlayer::releaseDspObjects()
{
  RAPT::rsAssert(dspPool, "This pointer should be assigned soon after creation");

  // Return the DSP objects to the pool:
  for(int i = 0; i < dspChain.getNumProcessors(); i++)
    dspPool->processorPool.repositProcessor(dspChain.getProcessor(i));
  dspChain.clear();

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

void rsSamplerEngine::RegionPlayer::allocateMemory()
{
  modulators.reserve(8);
  modMatrix.reserve(32);
  dspChain.reserve(8);
  // These values are ad-hoc arbitrarily chosen. Maybe give the function some parameters to choose
  // how many of each should be pre-allocated. It should be enough to avoid having to allocate more
  // later
}

rsReturnCode RegionPlayer::prepareToPlay(
  double fs, bool groupsOverride, bool regionsOverride)
{
  RAPT::rsAssert(isPlayable(region));  // This should not happen. Something is wrong.
  RAPT::rsAssert(stream != nullptr);   // Ditto.

  // The DSPs required by the group opcodes need to be in the RegionPlayer, if group-opcodes are
  // merely fallback values for absent region values which is indicated by the regionsOverride 
  // flag. So, of that flag is true, it means the group opcodes specify parameters for DSPs on the
  // RegionPlayer which means we must add the group DSPs here:
  bool withGroupDsps = regionsOverride;

  // Likewise, the DSPs required by the global instrument opcodes also need to be integrated into 
  // the RegionPlayer, iff these settings are merely fallback values (as opposed to values for
  // independent DSPs that are applied to the whole instrument):
  bool withInstrumDsps = groupsOverride;
  // At the moment, either both of these flags are true or both are false and i'm not yet sure if 
  // it will ever make sense to have these adjustable seperately. The behavior when one is false 
  // and the other true might be complicated to implement and/or understand for the user. And 
  // within the rsSamplerEngine, they are always both true. It's only in subclass rsSamplerEngine2
  // where we may switch ot the other kind of behavior (accumulation instead of fallback)


  if(!buildProcessingChain(withGroupDsps, withInstrumDsps)) 
  {
    releaseDspObjects();
    return rsReturnCode::layerOverload; 
  }
  if(!setupModulations()) {
    releaseDspObjects();
    return rsReturnCode::layerOverload; }
  resetPlayerSettings(); 
  setupDspSettingsFor(region, fs, groupsOverride, regionsOverride);
  // todo: move fs before the override parameters for consistency

  // todo: setup modulators and modulation connections

  dspChain.prepareToPlay(fs);
  // modulators.prepareToPlay(fs)

  // ToDo:
  // -resetDspState should reset only the state of the sample-player and be renamed accordingly
  //  because for the DSP objects that call comes to early, i.e. before setup() for the core
  //  objects is called. This may lead to wrong reset behavior if this behavior depends on the
  //  settings (which is the case for the Filter).
  // -The objects in the dspChain should reset themselves in prepareToPlay

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
  if( sampleTime >= endTime )                   // new
    return true;

  // todo:
  // -check also, if the amplitude envelope has reached its end
  // -hmm - maybe, if we allow the frequency envelope to go negative, we could also move 
  //  backward through the sample, so having reached the end of the stream may not actually be an
  //  appropriate condition. Or maybe, we should allow more general time-warping envelopes. 
  //  We'll see

  return false;
}

SignalProcessor* RegionPlayer::getProcessor(DspType type)
{
  RAPT::rsAssert(dspPool, "This pointer should be assigned soon after creation");
  return dspPool->processorPool.grabProcessor(type);
}

bool rsSamplerEngine::RegionPlayer::buildProcessingChain(bool withGroupDsps, bool withInstrumDsps)
{
  RAPT::rsAssert(dspChain.isEmpty(), "Someone has not cleaned up after finishing playback!");
  dspChain.clear(); // ...so we do it here. But this should be fixed elsewhere!

  using TypeChain = std::vector<DspType>; // maybe rename to DspTypeArray...or get rid

  // Adds the DSPs of the given types to the chain of actual DSP objects, if needed. Adding a 
  // particular DSP is needed, if no suitable such DSP is already there where "suitable" means:
  // "with right type and index":
  auto addDspsIfNeeded = [this](const TypeChain& dspTypeChain)
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
          index++; }
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
        dspChain.clear();
        return false; 
        // Not enough DSPs of desired type are available in the pool. We roll back the partial 
        // chain that we may have built already and report failure.
      }
    }
    return true;
    // When we arrive here, we have successfully finished the loop which means that we either could
    // add enough DSPs of the desired types to the chain or they were already present before. In 
    // both cases, the dspChain is now in the required state so we can report success.
  };

  // Add all the required DSPs (This is ugly! Can we do this in a more elegant way?):
  bool ok;
  ok = addDspsIfNeeded(region->getProcessingChain()); 
  if(!ok) return false;

  if(withGroupDsps)
  {
    ok = addDspsIfNeeded(region->getGroup()->getProcessingChain());
    if(!ok) return false;
  }
  if(withInstrumDsps)
  {
    ok = addDspsIfNeeded(region->getGroup()->getInstrument()->getProcessingChain());
    if(!ok) return false;
  }

  return true;
  // If false is returned, it means we do not have enough processors of the required types 
  // available. In this case, the caller should roll back and discard the whole RegionPlayer 
  // object. We either play a region correctly or not at all. This is an error condition that could
  // conceivably arise in normal usage (because we did not pre-allocate enough DSPs), so we should 
  // be able to handle it gracefully. It should actually not happen, i.e. we should make sure to 
  // always pre-allocate enough DSPs - but that may be impractical to ensure in a 100% airtight 
  // manner. But let's try at least to make that an exception that occurs only in extreme 
  // scenarios.
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

void RegionPlayer::setupDspSettingsFor(
  const Region* r, double fs, bool groupsOverride, bool regionsOverride)
{
  // To set up the settings, we call setupDspSettings 3 times to:
  //   (1) set up the general instrument-wide settings
  //   (2) set up group specific settings (this may override instrument settings)
  //   (3) set up region specific settings (this may override group and/or instrument settings)
  // but only if the respective flags are true. The flag groupsOverride means that the instrument
  // settings act as fallback settings for group-settings, so the instrument settings must be 
  // applied before applying group settings. Likewise for the regionsOverride flag: it means the 
  // group settings act as fallback values for region settings and must be set up before those:
  if(groupsOverride)
    setupDspSettings(region->getGroup()->getInstrument()->getSettings(), fs, true);
  if(regionsOverride)
    setupDspSettings(region->getGroup()->getSettings(), fs, groupsOverride);
  setupDspSettings(region->getSettings(), fs, regionsOverride);

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
  const std::vector<PlaybackSetting>& settings, double fs, bool overrideOldSetting)
{
  using PS = PlaybackSetting;
  using TP = Opcode;               // rename to OC

  double  amp        = 1.0;  // raw factor, computed "volume" opcode which is given in dB
  double  pan        = 0.0;  // -100...+100
  double  tuneCoarse = 0.0;  // in semitones
  double  tuneFine   = 0.0;  // in cents
  PanRule panRule    = PanRule::linear;
  //int    offset     = 0;
  bool   onTop       = !overrideOldSetting; // maybe get rid 

  // Loop through the settings of the region and for each setting that is present, change the 
  // value from its default to the stored value:
  for(size_t i = 0; i < settings.size(); i++)
  {

    PlaybackSetting setting = settings[i];
    TP type = setting.getType();                  // rename to opcode
    double val = (double) setting.getValue();
    switch(type)
    {
    // Amp settings:
    case TP::Volume:  { amp      = RAPT::rsDbToAmp(val); } break;
    case TP::Pan:     { pan      = val;                  } break;
    case TP::PanLaw:  { panRule  = (PanRule)(int)val;    } break;

    // Pitch settings:
    //case TP::PitchKeyCenter: { rootKey    = val; } break;  // done by caller
    case TP::Transpose:      { tuneCoarse = val; } break;
    case TP::Tune:           { tuneFine   = val; } break;

    //
    case TP::Delay:
    { 
      if(onTop) sampleTime += -val * fs; 
      else      sampleTime  = -val * fs; 
    }  break;

    case TP::Offset:
    { 
      if(onTop) offset += float(val); 
      else      offset  = float(val); 
    }  break;



    // Filter settings:
    case TP::filN_type:    { setupProcessorSetting(setting); } break;
    case TP::cutoffN:      { setupProcessorSetting(setting); } break;
    case TP::resonanceN:   { setupProcessorSetting(setting); } break;

    // Waveshaper settings:
    case TP::DistShape:  { setupProcessorSetting(setting); } break;
    case TP::DistDrive:  { setupProcessorSetting(setting); } break;
    case TP::DistOffset: { setupProcessorSetting(setting); } break;


    // ToDo: order the opcodes in the enum according to their type, such that we can here write
    // something like: 
    //   if( type > Opcode::dspStart && type < Opcode:dspEnd )
    //     setupProcessorSetting(setting);
    // instead of listing them all one-by-one, we should have tags: playerStart/End, dspStart/End
    // modulatorsStart/End, modConnectionsStart/End

    // Equalizer settings:
    case TP::eqN_gain: { setupProcessorSetting(setting); } break;
    case TP::eqN_freq: { setupProcessorSetting(setting); } break;
    case TP::eqN_bw:   { setupProcessorSetting(setting); } break;

      // ToDo: Get rid of this repetition! Maybe these can all be caught in a "default" branch?
      // but then what about the modulation settings? i think, we need a 3-way branch based on
      // whether a setting applies to the sample-player, dsp-chain or a modulator and then branch
      // further, if needed

    }
  }
  double tune   = tuneCoarse + 0.01 * tuneFine;
  double factor = pow(2.0, tune / 12.0);
  if(onTop) this->increment *= factor;
  else      this->increment  = factor;


  // From the computed local amp/pan/panRule variables, compute the amp member (which is
  // a rsFloat64x2)
  double t1, t2;  // temporaries
  switch(panRule)
  {
  case PanRule::linear:
  {
    t1 = (pan/200.0) + 0.5; // -100..+100 -> 0..1
    t2 = 1.0 - t1;
    if(onTop) this->amp *= 2.0 * amp * rsFloat64x2(t2, t1);
    else      this->amp  = 2.0 * amp * rsFloat64x2(t2, t1);
    // i'm not sure about the factor 2 -> check against sfz+ ..such a factor might be undesirable
    // in "onTop" mode: when 3 panners pan hard-left or right, we'll actually get a boost of 8
  } break;
  case PanRule::sinCos:
  {
    RAPT::rsError("not yet implemented");
  } break;
  }

  // ToDo: factor out the repetitive if(onTop)...into a local helper function 
  // set(double& setting, double value, bool onTop, bool multiplicative = false)
  // ...or maybe the multiplicative bool should be an int with values: 0: additive, 
  // 1: multiplicative, 2: dunno yet, ...

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
  // -Implement position and width opcodes. Maybe we should maintain a 2x2 gain matrix as member
  //  and all the different amp, pan, width, pos, etc. settings accumulate into this matrix. But 
  //  always compare results to sfz+ which serves as reference engine
}

void RegionPlayer::setupProcessorSetting(const PlaybackSetting& s)
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















}}