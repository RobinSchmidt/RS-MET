
//-------------------------------------------------------------------------------------------------
// construction/destruction:

ImmediatePlaybackAudioSource::ImmediatePlaybackAudioSource() : audioSampleBuffer(2,1024)
{
  audioSampleBuffer.clear();
  playPosition = 0;
  isPaused     = false;
}

ImmediatePlaybackAudioSource::~ImmediatePlaybackAudioSource()
{
  releaseResources();
}

//-------------------------------------------------------------------------------------------------
// overrides for subclasses AudioSource and PositionableAudioSource:

void ImmediatePlaybackAudioSource::prepareToPlay(int samplesPerBlockExpected, double sampleRate)
{
  int dummy = 0;
}

void ImmediatePlaybackAudioSource::releaseResources()
{

  int dummy = 0;
}

void ImmediatePlaybackAudioSource::getNextAudioBlock(const AudioSourceChannelInfo &bufferToFill)
{
  // make sure that the content of the preListenBuffer (in particular, the allocated memory) does
  // not change during reading it out:
  const ScopedLock scopedLock(lock); 

  // clear the region to be filled in the passed buffer:
  bufferToFill.clearActiveBufferRegion();

  if( isPaused )
    return;

  // determine the number of samples and channels to fill:
  int numChannels = jmin(bufferToFill.buffer->getNumChannels(),audioSampleBuffer.getNumChannels());
  int numSamples  = jmin(bufferToFill.numSamples, audioSampleBuffer.getNumSamples()-playPosition);

  // copy the a chunk of data from the preListenBuffer into an appropriate chunk of the output 
  // buffer (must loop over the channels):
  int c;
  for(c=0; c<numChannels; c++)
  {
    bufferToFill.buffer->copyFrom(c, bufferToFill.startSample, audioSampleBuffer, c, 
      playPosition, numSamples);
  }

  // if the output buffer has more channels than the buffer to be pre-listened, copy the content of 
  // the last channel in the preListenBuffer into the extra channels of the output buffer:
  for(c=numChannels; c<bufferToFill.buffer->getNumChannels(); c++)
  {
    bufferToFill.buffer->copyFrom(c, bufferToFill.startSample, audioSampleBuffer, numChannels-1,
      playPosition, numSamples);
  }

  // increment the read-position for the playback buffer:
  playPosition += numSamples;
}

//-------------------------------------------------------------------------------------------------
// others:

void ImmediatePlaybackAudioSource::startPlayback(AudioSampleBuffer *newBufferToPlay)
{
  // make sure that the content of the preListenBuffer (in particular, the allocated memory) is not
  // read out during we change it:
  const ScopedLock scopedLock(lock); 

  // re-allocate memory and reset the read-position:
  audioSampleBuffer.setSize(newBufferToPlay->getNumChannels(), newBufferToPlay->getNumSamples());
  for(int c=0; c<newBufferToPlay->getNumChannels(); c++)
    audioSampleBuffer.copyFrom(c, 0, *newBufferToPlay, c, 0, newBufferToPlay->getNumSamples());
  playPosition = 0;
}

void ImmediatePlaybackAudioSource::pausePlayback(bool shouldBePaused)
{
  isPaused = shouldBePaused;
}

void ImmediatePlaybackAudioSource::resumePlayback()
{
  isPaused = false;
}



