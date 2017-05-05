#ifndef rojue_AudioFileBufferUser_h
#define rojue_AudioFileBufferUser_h

#include "rojue_AudioFileBuffer.h"

namespace rojue
{

  /**

  This class is a baseclass for all classes that use some AudioFileBuffer object. The idea is 
  that the actual AudioFileBuffer is stored in some central repository and users can attach to
  the buffer and use it in some way (for example playing it back or displaying it on the screen). 
  The AudioFileBufferUser baseclass will take care of registering with the AudioFileBuffer that is 
  being used in order to allow the buffer to keep track of its referencers (for example to check in 
  its destructor that it is not referenced anymore and later maybe also for garbage collection and 
  deletion-notification purposes).

  */

  class AudioFileBufferUser
  {

    friend class AudioClipComponent;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AudioFileBufferUser(AudioFileBuffer *newBufferToUse);   

    /** Copy constructor. */
    AudioFileBufferUser(const AudioFileBufferUser &otherUser);

    /** Destructor. */
    virtual ~AudioFileBufferUser(); 

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the AudioFileBuffer object to be used. The baseclass implementation takes care of 
    updating the reference counter in the AudioFileBuffer object. Your subclasses will likely want 
    to override this to do their additional operations, which is well and good - but they should 
    always call the baseclass'es implementation also in their own implementations. */
    virtual void assignAudioFileBuffer(AudioFileBuffer *newBuffer);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of samples in the audioFileBufferToUse. */
    virtual int getNumSamples() const;

    /** Returns the number of channels in the audioFileBufferToUse. */
    virtual int getNumChannels() const;

    /** Retruns the sample-rate of the underlying file. */
    virtual double getFileSampleRate() const;

    /** Returns the length of the audiofile in seconds (assuming playback at the natural samplerate 
    of the file). */
    virtual double getBufferLengthInSeconds() const;

    /** Returns the file on which the audioFileBufferToUse is based on. */
    virtual const File getUnderlyingFile() const;

    /** Returns the name of the file on which the audioFileBufferToUse is based on. */
    virtual const juce::String getFileName() const;

    /** Returns the name of the file on which the audioFileBufferToUse is based on without 
    extension. */
    virtual const juce::String getFileNameWithoutExtension() const;

    /** Returns true if the bufferToCheck is the same as the buffer that is used by this object, 
    false otherwise. */
    virtual const bool isUsingBuffer(AudioFileBuffer *bufferToCheck);

    /** Returns a pointer to the AudioFileBuffer that is used by this object. IMPORTANT: wrap the 
    retrieval and everything you subsequently do with this pointer into locks, i.e. calls to 
    lockUsedBufferPointer()/unlockUsedBufferPointer(). If you are also planning to access the 
    actual audio data, use acquire/realeaseAudioFileBufferReadLock(), if you are planning to change 
    the size of the buffer (which potentially involves memory-reallcoation), use 
    acquire/realeaseAudioFileBufferWriteLock(). */
    virtual AudioFileBuffer* getUsedBuffer() const;

    /** Acquires the lock for the pointer to the underlying AudioFileBuffer and then delegates to
    AudioFileBuffer::getMinMaxSamples(), or assigns both values (minValue, maxValue) to zero in 
    case that this pointer is NULL. */
    virtual void getMinMaxSamples(int channel, const int startSample, const int numSamples,
      double &minValue, double &maxValue) const
    {
      ScopedLock pointerLock(audioFileBufferPointerLock);
      if( bufferToUse == NULL )
      {
        minValue = 0.0;
        maxValue = 0.0;
      }
      else
        bufferToUse->getMinMaxSamples(channel, startSample, numSamples, minValue, maxValue);
    }

    //---------------------------------------------------------------------------------------------
    // thread syncronization:

    /** Enters the lock for the member audioFileBufferPointerLock. */
    virtual void lockUsedBufferPointer();

    /** Exits the lock for the member audioFileBufferPointerLock. */
    virtual void unlockUsedBufferPointer();

    /** Enters the lock for the member audioFileBufferPointerLock (which locks accesses to the 
    pointer *bufferToUse) and acquires a read-lock for the audio-data that is stored in the buffer 
    to which this pointer points (in that order). When holding these locks, you can safely do some 
    work with the buffer but do not assign a new object to the pointer and don't change the 
    buffer's size. */
    virtual void acquireAudioFileBufferReadLock();

    /** Exits the read-lock for the audio-data that is stored in the buffer to which our member 
    *bufferToUse points and then exits the audioFileBufferPointerLock (which locks accesses to the 
    pointer *bufferToUse). The pointer-member *bufferToUse should point to the same object as it 
    did when you called acquireAudioFileBufferReadLock() for this object. */
    virtual void releaseAudioFileBufferReadLock();

    /** Enters the lock for the member audioFileBufferPointerLock (which locks accesses to the 
    pointer *bufferToUse) and acquires a write-lock for the audio-data that is stored in the buffer 
    to which this pointer points (in that order). When holding these locks, you can safely do some 
    work with the buffer (including changing its size) but do not assign a new object to the 
    pointer. */
    virtual void acquireAudioFileBufferWriteLock();

    /** Exits the read-lock for the audio-data that is stored in the buffer to which our member 
    *bufferToUse points and then exits the audioFileBufferPointerLock (which locks accesses to the 
    pointer *bufferToUse). The pointer-member *bufferToUse should point to the same object as it 
    did when you called acquireAudioFileBufferWriteLock() for this object. */
    virtual void releaseAudioFileBufferWriteLock();

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** This lock should used whenever the pointer-member AudioFileBuffer* bufferToUse is 
    assigned or derefeneced. It does not lock the contents of the buffer. */
    CriticalSection audioFileBufferPointerLock;

    /** The AudioFileBuffer which is to be used. */
    AudioFileBuffer *bufferToUse;

    /** A pointer to a buffer before some client code acquired the lock via the respective 
    functions (acquireAudioFileBufferReadLock(), etc.) to facilitate a check that the pointer 
    didn't change when the client code releases the lock again 
    (releaseAudioFileBufferReadLock(), etc.) - this check is for debug purposes only. */
    AudioFileBuffer *bufferBeforeLocking;

  };

}

#endif 
