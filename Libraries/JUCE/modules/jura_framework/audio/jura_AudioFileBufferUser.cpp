#include "rojue_AudioFileBufferUser.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AudioFileBufferUser::AudioFileBufferUser(AudioFileBuffer *newBuffer) 
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  bufferToUse         = NULL;
  bufferBeforeLocking = NULL;
  assignAudioFileBuffer(newBuffer);
}

AudioFileBufferUser::AudioFileBufferUser(const AudioFileBufferUser &otherUser)
{
  ScopedLock pointerLockForThis(audioFileBufferPointerLock);
  ScopedLock pointerLockForOther(otherUser.audioFileBufferPointerLock);
  bufferBeforeLocking = NULL;
  assignAudioFileBuffer(otherUser.bufferToUse);
}

AudioFileBufferUser::~AudioFileBufferUser()
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
  {
    ScopedWriteLock dataLock(bufferToUse->audioDataReadWriteLock);
    bufferToUse->deRegisterUser(this);
    bufferToUse = NULL;
  }
}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioFileBufferUser::assignAudioFileBuffer(AudioFileBuffer *newBuffer)
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( newBuffer != bufferToUse )
  {
    // de-regitser from the old assigned buffer:
    if( bufferToUse != NULL )
      bufferToUse->deRegisterUser(this);

    bufferToUse = newBuffer;

    // regitser to the new assigned buffer:
    if( bufferToUse != NULL )
      bufferToUse->registerUser(this);
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int AudioFileBufferUser::getNumSamples() const
{
  int result = 0;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getNumSamples();
  return result;
}

int AudioFileBufferUser::getNumChannels() const
{
  int result = 0;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getNumChannels();
  return result;
}

double AudioFileBufferUser::getFileSampleRate() const
{
  double result = 0.0;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getFileSampleRate();
  return result;
}

double AudioFileBufferUser::getBufferLengthInSeconds() const
{
  double result = 0.0;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getLengthInSeconds();
  return result;
}

const File AudioFileBufferUser::getUnderlyingFile() const
{
  File result = File::nonexistent;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getUnderlyingFile();
  return result;
}

const String AudioFileBufferUser::getFileName() const
{
  String result = String::empty;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getFileName();
  return result;
}

const String AudioFileBufferUser::getFileNameWithoutExtension() const
{
  String result = String::empty;
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToUse != NULL )
    result = bufferToUse->getFileNameWithoutExtension();
  return result;
}

const bool AudioFileBufferUser::isUsingBuffer(AudioFileBuffer *bufferToCheck)
{
  ScopedLock pointerLock(audioFileBufferPointerLock);
  if( bufferToCheck == bufferToUse )
    return true;
  else 
    return false;
}

AudioFileBuffer* AudioFileBufferUser::getUsedBuffer() const
{
  return bufferToUse;
}

//-------------------------------------------------------------------------------------------------
// thread syncronization:

void AudioFileBufferUser::lockUsedBufferPointer()
{
  audioFileBufferPointerLock.enter();
}

void AudioFileBufferUser::unlockUsedBufferPointer()
{
  audioFileBufferPointerLock.exit();
}

void AudioFileBufferUser::acquireAudioFileBufferReadLock()
{
  audioFileBufferPointerLock.enter();
  bufferBeforeLocking = bufferToUse;
  if( bufferToUse != NULL )
    bufferToUse->acquireReadLock();
}

void AudioFileBufferUser::releaseAudioFileBufferReadLock()
{
  jassert( bufferToUse == bufferBeforeLocking ); 
    // don't assign a new buffer in between acquireAudioFileBufferReadLock() and 
    // releaseAudioFileBufferReadLock() - you will then release a lock fo a buffer which is 
    // different from the one, you acquired it for

  if( bufferToUse != NULL )
    bufferToUse->releaseReadLock();
  audioFileBufferPointerLock.exit();
}

void AudioFileBufferUser::acquireAudioFileBufferWriteLock()
{
  audioFileBufferPointerLock.enter();
  bufferBeforeLocking = bufferToUse;
  if( bufferToUse != NULL )
    bufferToUse->acquireWriteLock();
}

void AudioFileBufferUser::releaseAudioFileBufferWriteLock()
{
  jassert( bufferToUse == bufferBeforeLocking ); 
    // don't assign a new buffer in between acquireAudioFileBufferWriteLock() and 
    // releaseAudioFileBufferWriteLock() - you will then release a lock fo a buffer which is 
    // different from the one, you acquired it for

  if( bufferToUse != NULL )
    bufferToUse->releaseWriteLock();
  audioFileBufferPointerLock.exit();
}
