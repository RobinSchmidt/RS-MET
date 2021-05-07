//=================================================================================================
// Function definitions for the helper classes:

bool AudioFileStreamPreloaded::setData(
  float** newData, int numFrames, int numDataChannels, float sampleRate, int numStreamChannels, 
  const std::string& uniqueName)
{
  // Deallocate old and allocate new memory:
  clear();
  int numChannelsMin = rsMin(numDataChannels, numStreamChannels);
  flatData = new float[numChannelsMin*numFrames];
  channelPointers = new float*[numStreamChannels];
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
  return true; // success
}

/*
void AudioFileStreamPreloaded::setNumOutputChannels(int newNumChannels)
{
  int dummy = 0;
}
*/

void AudioFileStreamPreloaded::clear()
{
  numChannels = 0;
  numFrames   = 0; 

  delete[] channelPointers;
  channelPointers = nullptr; 

  delete[] flatData;
  flatData = nullptr;
}

//-------------------------------------------------------------------------------------------------

void SamplePool::clear()
{
  for(size_t i = 0; i < samples.size(); i++)
    delete samples[i];
  samples.clear();
}

//-------------------------------------------------------------------------------------------------
// rsDataSFZ:

int rsDataSFZ::Group::addRegion()
{
  rsDataSFZ::Region* r = new rsDataSFZ::Region;
  //r->group = this;
  r->parent = this;
  regions.push_back(r);
  return ((int) regions.size()) - 1;
}

int rsDataSFZ::Group::getRegionIndex(const rsDataSFZ::Region* region) const
{
  for(size_t i = 0; i < regions.size(); i++)
    if(regions[i] == region)
      return (int) i;
  return -1;
}

rsDataSFZ::Region* rsDataSFZ::Group::getRegion(int i) const
{
  if(i < 0 || i >= (int)regions.size()) {
    rsError("Invalid region index");
    return nullptr; 
  }
  return regions[i];
}

void rsDataSFZ::Group::clearRegions()
{
  for(size_t i = 0; i < regions.size(); i++)
    delete regions[i];
  regions.clear();
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
  float** data, int numFrames, int numChannels, float sampleRate, const std::string& uniqueName)
{
  // todo: 
  // -check, if a sample with the same uniqueName already exists - if so, we have nothing to 
  //  do and may return early with an appropriate code
  // if(isSampleInPool(..)) return ReturnCode::nothingToDo;

  // -maybe the Streamer object should always have two channel pointers, but in case of 
  //  mono-samples, both just point to the same buffer. that may make it easier to handle things
  //  uniformly

  AudioFileStreamPreloaded* stream = new AudioFileStreamPreloaded;
  bool allocOK = stream->setData(data, numFrames, numChannels, sampleRate, 2, uniqueName);
  if(allocOK == false)
    return ReturnCode::memAllocFail;
  return samplePool.addSample(stream);
}

int rsSamplerEngine::addGroup()
{
  Group g;
  groups.push_back(g);
  return ((int) groups.size()) - 1;
}

int rsSamplerEngine::addRegion(int gi, uchar loKey, uchar hiKey)
{
  if(gi < 0 || gi >= (int)groups.size()) {
    rsError("Invalid group index");
    return ReturnCode::invalidIndex; 
  }
  int ri = groups[gi].addRegion();      // region index within its group

  // Add the region to the regionsForKey in between loKey and hiKey
  const Region* region = getRegion(gi, ri);
  for(uchar k = loKey; k <= hiKey; k++)
    addRegionForKey(k, region);

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
  Region* r = getRegion(gi, ri);
  //r->setSampleStream(samplePool.getSampleStream(si));
  r->setCustomPointer(samplePool.getSampleStream(si));
  return ReturnCode::success;
}

int rsSamplerEngine::setRegionSetting(Region* region, PlaybackSetting::Type type, float value)
{
  // todo: 
  // -check, if the passed region is actually a member of any of the groups - maybe
  //  rsAssert(isRegionValid(region)) or something
  //  -maybe have a function findRegion(const Region* region, int* groupIndex, int* regionIndex)
  // -before pushing the setting, we should check, whether the setting is already present and if
  //  so, just change the value and not push a 2nd, conflicting value

  //region.settings.push_back(PlaybackSetting(type, value));

  int gi, ri;
  findRegion(region, &gi, &ri);
  if(gi == -1 || ri == -1)
    return ReturnCode::notFound;
  // todo: simplify to hasRegion(region)

  Region* r = getRegion(gi, ri);
  // That's a bit dirty - it actually has the same effect as casting away the const from the passed
  // region. The caller may assume, that the passed region remins constant when in fact, it 
  // doesn't. Maybe make the region parameter non-const and prevent the client code from modifying
  // regions directly by other means like making all setters private such that only the 
  // rsSamplerEngine friend class can set them

  rsAssert(r == region);  // todo: operate directly on the passed region

  r->settings.push_back(PlaybackSetting(type, value));
  // Preliminary. We need to figure out, if that setting already exists and if so, just change its
  // value instead of pushing another value for the same parameter

  return ReturnCode::success;
  // Maybe we should distinguish between settingAdded and settingModified in the return code. But
  // actually, from the caller's perspective, that information should be irrelevant anyway, so 
  // maybe not.
}

rsSamplerEngine::Region* rsSamplerEngine::getRegion(int groupIndex, int regionIndex)
{
  int gi = groupIndex, ri = regionIndex;
  if(gi < 0 || gi >= (int)groups.size()) {
    rsError("Invalid group index");
    return nullptr; 
  }
  return groups[gi].getRegion(ri);
}

//-------------------------------------------------------------------------------------------------
// Processing:

void rsSamplerEngine::processFrame(float* left, float* right)
{
  rsFloat64x2 out = 0.0;
  for(size_t i = 0; i < activePlayers.size(); i++)
    out += activePlayers[i]->getFrame();
  *left  = (float) out[0];
  *right = (float) out[1];
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

void rsSamplerEngine::addRegionForKey(uchar k, const Region* region)
{
  regionsForKey[k].addRegion(region);
  // What, if the region is already there? We should check that before. It's probably not supposed
  // to happen, but anyway. Well...maybe it is, when the user tweaks loKey/hiKey settings on a GUI.
  // On the other hand, it may actually be useful to be able to duplicate regions, especially on
  // a GUI: duplicate-and-edit could be a common workflow
}

void rsSamplerEngine::findRegion(const rsSamplerEngine::Region* r, int* gi, int* ri)
{
  *gi = -1;
  *ri = -1;
  for(size_t i = 0; i < groups.size(); i++) {
    int j = groups[i].getRegionIndex(r);
    if(j != -1) {
      *gi = (int) i;
      *ri = j;
      return; }}
  rsError("Region not found");
  // A region should always be found in one (and only one group). If we don't find it, it means
  // the caller has passed a pointer to a region object that is not part of this instrument. If 
  // this happens, it indicates a bug at the call site.
}

rsSamplerEngine::RegionPlayer* rsSamplerEngine::getRegionPlayerFor(const Region* r)
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
  rp->setRegionToPlay(r);
  activePlayers.push_back(rp);
  return rp;
}

const AudioFileStream* rsSamplerEngine::getSampleStreamFor(const Region* r)
{
  //return r->getSampleStream();
  return (const AudioFileStream*) r->getCustomPointer();
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
    RegionPlayer* rp = getRegionPlayerFor(r);
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
      numRegions++;
  }
  return ReturnCode::success;
  // Another possibility for the return value would have been to return the number of voices that
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

  // ToDo:
  // -loop through all activePlayers to find those who are playing the given key
  // -mark them for going into release state
  // -later: trigger all regions that should play a release sample for the given key/vel combo


  return ReturnCode::success;
}



//-------------------------------------------------------------------------------------------------
// rsSamplerEngine::RegionPlayer

void rsSamplerEngine::RegionPlayer::setRegionToPlay(const rsSamplerEngine::Region* regionToPlay)
{
  region = regionToPlay;
  stream = getSampleStreamFor(region);
  prepareToPlay();
}

rsFloat64x2 rsSamplerEngine::RegionPlayer::getFrame()
{
  if(sampleTime < 0) {               // Negatively initialized sampleTime implements delay.
    sampleTime++;                    // We just increment the time and return 0,0. Actual output
    return rsFloat64x2(0.0, 0.0); }  // will be produced as soon as sampleTime reaches zero.  

  //float tmp[2];
  //stream->getFrame(sampleTime, tmp);

  float L, R;
  stream->getFrameStereo(sampleTime, &L, &R);



  // more stuff to do:
  // -apply pitch envelope and lfo
  // -implement interpolation
  // -apply the DSP processes


  sampleTime++; 
  // sampleTime should probably be of type double and we should have use an increment that depends 
  // on the current pitch, i.e. determined by pitchKeyCenter, noteFreq, sampleRates of the file and 
  // output.


  return this->amp * rsFloat64x2(L, R); 


  //return rsFloat64x2(tmp[0], tmp[1]);  // preliminary

  //return rsFloat64x2(0.0, 0.0);  // preliminary
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
  //ok &= region->getGroup()->getInstrument() != nullptr;  // uncomment
  return ok;
}

void rsSamplerEngine::RegionPlayer::prepareToPlay()
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
  //setupDspSettings(region->getGroup()->getInstrument()->getSettings()); // uncomment!
  setupDspSettings(region->getGroup()->getSettings());
  setupDspSettings(region->getSettings());

  // return rsSamplerEngine::ReturnCode::success;
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
  sampleTime = 0;
  // ampEnv.setToDefaults(), etc.
  // ...more to do... 
}

void rsSamplerEngine::RegionPlayer::setupDspSettings(const std::vector<PlaybackSetting>& settings)
{
  // Loop through the settings of the region and for each setting that is present, change the 
  // value from its default to the stored value:
  
  using PS = PlaybackSetting;
  using TP = PS::Type;

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
    case TP::Volume:  { amp     = rsDbToAmp(val); } break;
    case TP::Pan:     { pan     = val;            } break;
    case TP::PanRule: { panRule = (int)val;       } break;

      //case TP::FilterCutoff: { flt.setCutoff(val);  } break;

      // ...more to do...

    }
  }


  // ToDo: from the computed local amp/pan/panRule variables, compute the amp member (which is
  // a rsFloat64x2)
  double t1, t2, t3;  // temporaries
  switch(panRule)
  {
  case PS::PanRule::linear:
  {
    t1 = (pan/200.0) + 0.5; // -100..+100 -> 0..1
    t2 = 1.0 - t1;
    this->amp = 2.0 * amp * rsFloat64x2(t1, t2);
  } break;
  case PS::PanRule::sinCos:
  {
    rsError("not yet implemented");
  } break;
  }

 




  // ToDo:
  // -Maybe within the switch statement set up some flags that indicate, if a particular setting is
  //  used. If the flag is false, we may skip the associated DSP process in getFrame/processBlock. 
  //  We may need inquiry functions such as hasFilter, hasAmpEnv, hasPitchEnv, hasFilterEnv, 
  //  hasPitchLFO. But this makes things more complicated, so maybe it's not really a good idea.
}

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