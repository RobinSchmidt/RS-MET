//#include "romos_VoiceAllocator.h"
//using namespace romos;

VoiceAllocator voiceAllocator;  // definition of the global object

//-------------------------------------------------------------------------------------------------
// construction/destruction:

VoiceAllocator::VoiceAllocator()
{
  numVoices     = 16;
  stealingMode  = STEAL_OLDEST_VOICE;
  retriggerMode = true;
  reset();
}

VoiceAllocator::~VoiceAllocator()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void VoiceAllocator::setNumVoices(int newNumVoices)
{
  if( numVoices >= 1 && numVoices <= maxNumVoices )
  {
    reset();
    numVoices = newNumVoices;
  }
  else
    DEBUG_BREAK;
}

void VoiceAllocator::setVoiceStealingMode(int newVoiceStealingMode)
{
  if( newVoiceStealingMode >= 0 && newVoiceStealingMode < NUM_VOICE_STEALING_MODES )
    stealingMode = newVoiceStealingMode;
  else
    DEBUG_BREAK;
}

void VoiceAllocator::setRetriggerMode(bool shouldRetrigger)
{
  reset();
  retriggerMode = shouldRetrigger;
}

void VoiceAllocator::setNoteOnTriggerFlag(int voiceIndex)
{
  noteOnTriggerFlags |= ( maskForFlags >> voiceIndex);
}

void VoiceAllocator::setNoteOffTriggerFlag(int voiceIndex)
{
  noteOffTriggerFlags |= ( maskForFlags >> voiceIndex);
}

//-------------------------------------------------------------------------------------------------
// event handling:

int VoiceAllocator::noteOn(int key, int velocity)
{
  if( velocity == 0 )
    return noteOff(key);

  int playingIndex = -1;  // index in the array of playing voices
  int voiceToUse   = -1;

  if( retriggerMode == true )
  {
    playingIndex = findAmongPlayingVoices(key);
    if( playingIndex != -1 )
      voiceToUse = playingVoiceIndices[playingIndex];
  }

  if( voiceToUse == -1 )
  {
    if( numPlayingVoices < numVoices )
    {
      voiceToUse   = getFirstFreeVoice();
      playingIndex = numPlayingVoices;
      numPlayingVoices++;
      playingVoiceIndices[playingIndex] = voiceToUse;
    }
    else
    {
      playingIndex = getPlayingIndexToSteal(key);
      voiceToUse   = playingVoiceIndices[playingIndex];
      movePlayingVoiceToBottom(playingIndex);
    }
  }
  else
    movePlayingVoiceToBottom(playingIndex);

  voiceStates[voiceToUse].key                = key;
  voiceStates[voiceToUse].frequency          = RAPT::rsPitchToFreq(key);  // preliminary - use tuning tables later
  voiceStates[voiceToUse].normalizedVelocity = velocity / 127.0;
  voiceStates[voiceToUse].isPlaying          = true;

  setNoteOnTriggerFlag(voiceToUse);

  return voiceToUse;  // maybe make this void
}

int VoiceAllocator::noteOff(int key)
{
  int voiceIndex           = -1;
  int indexInPlayingVoices = findAmongPlayingVoices(key, 0, true);

  if( indexInPlayingVoices != -1 )
  {
    voiceIndex = playingVoiceIndices[indexInPlayingVoices];
    voiceStates[voiceIndex].normalizedVelocity = 0.0;
  }

  setNoteOffTriggerFlag(voiceIndex);
  return voiceIndex;  // make this void
}

void VoiceAllocator::killVoice(int voiceIndex)
{
  removeFromPlayingVoices(voiceIndex);
  resetVoice(voiceIndex);
}

void VoiceAllocator::reset()
{
  noteOnTriggerFlags  = 0;
  noteOffTriggerFlags = 0;
  numPlayingVoices    = 0;
  for(int i = 0; i < maxNumVoices; i++)
  {
    resetVoice(i);
    playingVoiceIndices[i] = -1;
  }
}

void VoiceAllocator::resetVoice(int voiceIndex)
{
  voiceStates[voiceIndex].key                = -1;
  voiceStates[voiceIndex].frequency          = 0.0;
  voiceStates[voiceIndex].normalizedVelocity = 0.0;
  voiceStates[voiceIndex].isPlaying          = false;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

int VoiceAllocator::getFirstFreeVoice()
{
  for(int i = 0; i < numVoices; i++)
  {
    if( voiceStates[i].isPlaying == false )
      return i;
  }
  return -1;
}

int VoiceAllocator::getPlayingIndexToSteal(int /*key*/)
{
  switch( stealingMode )
  {
  case STEAL_OLDEST_VOICE:  return getOldestPlayingVoiceIndex();
  case STEAL_NEWEST_VOICE:  return getNewestPlayingVoiceIndex();

  }
  return -1;
}

int VoiceAllocator::getOldestPlayingVoiceIndex()
{
  return 0;
}

int VoiceAllocator::getNewestPlayingVoiceIndex()
{
  return numPlayingVoices-1;
}

int VoiceAllocator::findAmongPlayingVoices(int key, int searchStartIndex,
  bool ignoreVoicesWithZeroVelocity)
{
  for(int i = searchStartIndex; i < numPlayingVoices; i++)
  {
    if( voiceStates[playingVoiceIndices[i]].key == key )
    {
      if( ignoreVoicesWithZeroVelocity == true
        && voiceStates[playingVoiceIndices[i]].normalizedVelocity == 0.0 )
      {
        // skip
      }
      else
        return i;
    }
  }
  return -1;
}

int VoiceAllocator::getPlayingIndexOfVoice(int voiceIndex)
{
  for(int i = 0; i < numPlayingVoices; i++)
  {
    if( playingVoiceIndices[i] == voiceIndex )
      return i;
  }
  return -1;
}

void VoiceAllocator::removeFromPlayingVoices(int voiceIndex)
{
  rassert(numPlayingVoices > 0); // for debugg

  // \todo: robustify this - when there are more than one VoiceKiller modules, the same voice might
  // be killed multiple times but maybe it's already robust due to
  // "if( indexInPlayingVoices >= 0 )"? ...check this

  int indexInPlayingVoices = getPlayingIndexOfVoice(voiceIndex);
  if( indexInPlayingVoices >= 0 )
  {
    moveUpPlayingVoices(indexInPlayingVoices+1);
    playingVoiceIndices[numPlayingVoices-1] = -1; // heap corruption in the ModularSynth when a note was triggered
    numPlayingVoices--;
  }
}

void VoiceAllocator::movePlayingVoiceToBottom(int playingIndex)
{
  int voiceIndex = playingVoiceIndices[playingIndex];
  moveUpPlayingVoices(playingIndex+1);
  playingVoiceIndices[numPlayingVoices-1] = voiceIndex;
}

void VoiceAllocator::moveUpPlayingVoices(int fromWhichIndex)
{
  for(int i = fromWhichIndex; i < numPlayingVoices; i++)
    playingVoiceIndices[i-1] = playingVoiceIndices[i];
}
