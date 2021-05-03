//-------------------------------------------------------------------------------------------------
// Setup:

// Shortcuts to reduce repetitive verbosity:
#define RS_TMPDEC template<class TSig, class TPar, class TSmp>  // template declarations
#define RS_SMPENG rsSamplerEngine<TSig, TPar, TSmp>             // sampler engine type

RS_TMPDEC int RS_SMPENG::addSampleToPool(
  TSmp** data, int numFrames, int numChannels, TPar sampleRate, const std::string& uniqueName)
{
  // todo: 
  // -check, if a sample with the same uniqueName already exists - if so, we have nothing to 
  //  do and may return early with an appropriate code
  // if(isSampleInPool(..)) return ReturnCode::nothingToDo;

  // -maybe the Streamer object should always have two channel pointers, but in case of 
  //  mono-samples, both just point to the same buffer. that may make it easier to handle things
  //  uniformly


  AudioFileStreamPreloaded* stream = new AudioFileStreamPreloaded;
  int result = stream->setData(data, numFrames, numChannels, sampleRate, uniqueName);
  samplePool.addSample(stream);
  return result;
}

RS_TMPDEC int RS_SMPENG::addGroup()
{
  Group g;
  groups.push_back(g);
  return ((int) groups.size()) - 1;
}

RS_TMPDEC int RS_SMPENG::addRegion(int gi, uchar loKey, uchar hiKey)
{
  if(gi < 0 || gi >= (int)groups.size()) {
    rsError("Invalid group index");
    return ReturnCode::invalidIndex; 
  }
  int ri = groups[gi].addRegion();      // regiin index within its group

  // todo: add the region to the regionsForKey in between loKey and hiKey

  //for(int i = loKey; i <= hiKey; i++)
  //  addRegionForKey(gi, ri);

  return ri;
}

//-------------------------------------------------------------------------------------------------
// Inquiry:

RS_TMPDEC bool RS_SMPENG::shouldRegionPlay(
  const Region* r, const char key, const char vel)
{
  return false; // preliminary
}

//-------------------------------------------------------------------------------------------------
// Processing:

RS_TMPDEC void RS_SMPENG::processFrame(TSig* frame)
{

}

RS_TMPDEC void RS_SMPENG::processBlock(TSig** block, int numFrames)
{

}

RS_TMPDEC void RS_SMPENG::handleMusicalEvent(const rsMusicalEvent<TPar>& ev)
{

}

//=================================================================================================
// Function definitions for the helper classes:

RS_TMPDEC int RS_SMPENG::AudioFileStreamPreloaded::setData(
  TSmp** newData, int numFrames, int numChannels, TPar sampleRate,
  const std::string& uniqueName)
{
  // Deallocate old and allocate new memory:
  clear();
  flatData = new TSmp[numChannels*numFrames];
  channelPointers = new TSmp*[numChannels];
  if(flatData == nullptr || channelPointers == nullptr) {
    clear(); return ReturnCode::memAllocFail; }

  // Copy the new data into the freshly allocated memory:
  for(int c = 0; c < numChannels; c++) {
    channelPointers[c] = &flatData[c*numFrames];
    for(int n = 0; n < numFrames; n++)
      channelPointers[c][n] = newData[c][n]; }
  // Maybe we should have a version of this function which does not need to copy data but instead
  // just takes over ownership of the passed array. But this would need a parameter for the flat 
  // data array, too. We'll see, how this meshes with the wavefile loading functions...

  // Update metadata members and report success:
  this->numChannels = numChannels;
  this->numFrames   = numFrames;
  this->sampleRate  = sampleRate;
  return ReturnCode::success;
}

RS_TMPDEC void RS_SMPENG::AudioFileStreamPreloaded::clear()
{
  numChannels = 0;
  numFrames   = 0; 

  delete[] channelPointers;
  channelPointers = nullptr; 

  delete[] flatData;
  flatData = nullptr;
}

//-------------------------------------------------------------------------------------------------

RS_TMPDEC int RS_SMPENG::Group::addRegion()
{
  RS_SMPENG::Region r;
  r.group = this;
  regions.push_back(r);
  return ((int) regions.size()) - 1;
}

//-------------------------------------------------------------------------------------------------

RS_TMPDEC void RS_SMPENG::SamplePool::clear()
{
  for(size_t i = 0; i < samples.size(); i++)
    delete samples[i];
  samples.clear();
}

#undef RS_TMPDEC
#undef RS_SMPENG





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
-We need a sort of sample-pool containing all our sample-buffers along with their playback 
 parameters.
-In th sfz spec, the performance parameters defined for the whole instrument and for groups work as
 fallback values that can be overriden by respective values on a lower level of the hierarchy. The
 drum-sampler needs them to work in an accumulative fashion. For example a cutoff defined for a 
 region should be used as is and the cutoff for the group should be for a 2nd filter through which 
 the signal of the whole group is passed. That's a different semantic. In the original sfz-spec, 
 the group cutoff would just have been overriden by the region setting (i think -> verify). Maybe 
 we should have switches like: goup/instrumentSettingsAccumulate
-Loop-points and start/end points shall be floating point numbers. When parsing an sfz file, we 
 just split the number string into pre-dot and post-dot part, parse them separately and store the 
 post-dot part in a double (or TPar) and the pre-dot part in an integer. That way, we won't suffer 
 precision loss for numbers bigger pre-dot part. This should be downward compatible with sfz spec 
 (which has only integer loop points (i think -> verify))


Notes:

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

may use SampleBuffer, SamplePlaybackParameters

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

*/