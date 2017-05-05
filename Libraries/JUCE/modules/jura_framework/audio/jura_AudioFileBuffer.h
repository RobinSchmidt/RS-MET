#ifndef jura_AudioFileBuffer_h
#define jura_AudioFileBuffer_h

class AudioFileBufferUser; // forward declaration

/**

This class represents an audio-file, the data of which is buffered within an object of this
class. This is done via inheriting from juce::AudioSampleBuffer, but in addition to the
functionality of the base-class, this class here stores also information about the underlying
file such as the relative path (with respect to some sample-content directory). It also provides
some infastructure to facilitate the use of this buffer for clips that share audiodata.

@see AudioFileBufferUser, AudioClip, WaveformDisplay

*/

class AudioFileBuffer : protected AudioSampleBuffer
{

  friend class AudioFileBufferUser;
  friend class AudioClip;
  friend class MixsonicContentComponent;
  // \todo: throw this away and aquire the mutex locks in the constructor of AudioClipComponent 
  // instead of MixsonicContentComponent


public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AudioFileBuffer(const File &fileToLoadFrom = File::nonexistent);

  /** Destructor. */
  virtual ~AudioFileBuffer() throw();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Loads the actual audio data from a file and returns whether this was successful. */
  virtual bool loadAudioDataFromFile(const File &fileToLoadFrom, bool showAlertBoxOnFail = false);

  /** When for some reason the path of the underlying was changed (i.e. the underlying file was
  moved), but the content of the file stays the same, you should call this function to update
  the stored path information for the underlying file. This will not trigger a reloading or
  anything of the file from the new location - it just stores the new path. */
  virtual void filePathChanged(const File& newFilePath);

  /** Registers a new user for this buffer (does nothing when the user is already registered). */
  virtual void registerUser(AudioFileBufferUser* userToRegister);

  /** De-registers a user from this buffer (does nothing when the user not registered). */
  virtual void deRegisterUser(AudioFileBufferUser* userToDeRegister);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** returns true when the buffer contains data from a valid audio file. */
  virtual bool isAudioFileValid() const;

  /** Returns the number of samples in the buffer. */
  virtual int getNumSamples() const;

  /** Returns the number of channels in the buffer. */
  virtual int getNumChannels() const;

  /** Retruns the sample-rate of the underlying file. */
  virtual double getFileSampleRate() const;

  /** Returns the length of the audiofile in seconds (assuming playback at the natural samplerate
  of the file). */
  virtual double getLengthInSeconds() const;

  /** Returns the file on which thies buffer is based on. */
  virtual const File getUnderlyingFile() const;

  /** Returns the name of the file on which thies buffer is based on. */
  virtual const juce::String getFileName() const;

  /** Returns the name of the file on which thies buffer is based on without the file
  extension. */
  virtual const juce::String getFileNameWithoutExtension() const;

  /** Copies data from this buffer into another AudioSampleBuffer, making sure, that this
  (source) buffer is accessed in a thread safe manner (but possibly not the target buffer). */
  virtual void copyTo(const int destChannel, const int destStartSample,
    AudioSampleBuffer &destBuffer, const int sourceChannel, const int sourceStartSample,
    int numSamples);

  /** Returns a pointer to a sample in one of the buffers channels - IMPORTANT: if you retrieve
  the pointer via this function, it might be invalidated by other threads during you use it - to
  avoid this, you must call aquireLockForBufferRead() before retrieving the pointer and when
  you're done, call releaseLockForBufferRead(). */
  virtual float* getSampleData(const int channelNumber, const int sampleOffset = 0);

  /** Returns a sample in a given channel at a given sample-number in a thread-safe manner (i.e.
  no other thread will change the size of this buffer during this call, give that all functions
  that attempt to change the size acquire the write-lock. */
  virtual float getSample(int channel, int sampleNumber) const
  {
    ScopedReadLock scopedLock(audioDataReadWriteLock);

    channel = jmin(channel, AudioSampleBuffer::getNumChannels()-1);
    if(sampleNumber >= AudioSampleBuffer::getNumSamples())
      return 0.f;

    //float* pointer = AudioSampleBuffer::getSampleData(channel, sampleNumber);
    const float* pointer = AudioSampleBuffer::getReadPointer(channel, sampleNumber);
    float  result  = *pointer;
    return result;
  }

  /** Returns the minimum and maximum samples in a given channel in a given interval in a
  thread-safe manner (i.e. no other thread will change the size of this buffer during this call,
  given that all functions that attempt to change the size acquire the write-lock.
  @see AudioSampleBuffer::findMinMax() */
  virtual void getMinMaxSamples(int channel, const int startSample, const int numSamples,
    double &minValue, double &maxValue) const
  {
    ScopedReadLock scopedLock(audioDataReadWriteLock);

    //float minFloat, maxFloat;
    //AudioSampleBuffer::findMinMax(channel, startSample, numSamples, minFloat, maxFloat);
    //minValue = (double)minFloat;
    //maxValue = (double)maxFloat;

    Range<float> r = AudioSampleBuffer::findMinMax(channel, startSample, numSamples);
    minValue = (double) r.getStart();
    maxValue = (double) r.getEnd();

    // maybe use getMinMaxSamplesWithoutLock here - avoid duplication
  }

  /** Same as getMinMaxSamples() but without acquiring the read-lock before - meant to be called
  inside a loop when you have acquired the lock already outside the loop (this is IMPORTANT -
  don't call this function when you have not acquired the lock!). */
  virtual void getMinMaxSamplesWithoutLock(int channel, const int startSample, const int numSamples,
    double &minValue, double &maxValue) const
  {
    //float minFloat, maxFloat;
    //AudioSampleBuffer::findMinMax(channel, startSample, numSamples, minFloat, maxFloat);
    //minValue = (double)minFloat;
    //maxValue = (double)maxFloat;

    Range<float> r = AudioSampleBuffer::findMinMax(channel, startSample, numSamples);
    minValue = (double) r.getStart();
    maxValue = (double) r.getEnd();
  }

  //---------------------------------------------------------------------------------------------
  // thread synchronization:

  /** Acquires the lock for accessing the audio-data in this buffer. The 'read'-term is not to be
  taken literally, you may also write, but don't change the size of the buffer (don't re-allocate
  memory). */
  virtual void acquireReadLock();

  /** Releases the lock for accessing the audio-data in this buffer. */
  virtual void releaseReadLock();

  /** Aquires the lock for accessing the audio-data in this buffer and also possibly changing the
  size of the buffer (re-allocating memory). */
  virtual void acquireWriteLock();

  /** Releases the lock for accessing and/or changing the size of this buffer. */
  virtual void releaseWriteLock();

  //---------------------------------------------------------------------------------------------
  // operators:

  /** Compares two AudioFileBuffers of equality - they regarded as equal when the filenames and
  times of last modification match. */
  bool operator==(const AudioFileBuffer& otherBuffer) const
  {
    ScopedReadLock scopeLock(audioDataReadWriteLock);

    if(theAudioFile.getFileName() != otherBuffer.theAudioFile.getFileName())
      return false;
    if(theAudioFile.getLastModificationTime()!=otherBuffer.theAudioFile.getLastModificationTime())
      return false;
    return true;

    /*
    // an alternative implementation that compares the actual contents of the files sample by
    // sample (this has never been tested yet)
    ScopedLock lock1(audioDataLock);
    ScopedLock lock2(otherBuffer.audioDataLock);
    if( info.numChannels != otherBuffer.info.numChannels )
    return false;
    if( info.numSamples != otherBuffer.info.numSamples )
    return false;
    for(int c=0; c<info.numChannels; c++)
    {
    float* ownSampleData   =             getSampleData(c, 0);
    float* otherSampleData = otherBuffer.getSampleData(c, 0);
    for(int n=0; n<info.numSamples; n++)
    {
    if( ownSampleData[n] != otherSampleData[n] )
    return false;
    }
    }
    // filenames match and file contents are sample-by-sample the same - we regard the two
    // buffers as equal:
    return true;
    */

    /*
    // this old implementation that compares the underlying file objects itself is not suitable
    // for Mixsonic because it will give false negatives for samples that are dragged into the
    // arrangement from the sample content directory and from the project directory:
    if( theAudioFile == otherBuffer.theAudioFile  )
    return true;
    else
    return false;
    */
  }

  /** Compares two AudioFileBuffers of inequality - they regarded as inequal when they point to
  different files. */
  bool operator!=(const AudioFileBuffer& otherBuffer) const
  {
    return !(*this == otherBuffer);
  }

protected:

  /** The root directory in which the audio files are generally stored. */
  juce::String audioFileRootDirectory;

  /** The audiofile to which this buffer refers. \todo: store this in AudioFileInfo */
  File theAudioFile;

  /** Some information about the audiofile such as sample-rate etc. */
  AudioFileInfo info;

  /** An array of objects which use this buffer - this is used for tracing the references to this
  buffer. It's purpose is currently only to provide a means for garbage collection (so a simple
  reference counter could have done the job also) - but knowing our referencors explicitly allows
  for extending it later to support callbacks to the users (for example when the buffer is edited
  destructively and users of this buffer must be notified about that). */
  juce::Array<AudioFileBufferUser*, CriticalSection> users;

  /** We use two separate mutex locks for accessing the audio data from the audio-thread and from
  the graphic thread in order to not not have them block each other because there is no good
  reason why not two threads should be able to read the same data. Only when the data is changed
  (i.e. a new audiofile is being loaded) we will aquire both locks. */
  //CriticalSection lockForAudioThread, lockForDrawingThread;

  /** Always aquire at least a read-lock before you access the data in the inherited
  AudioSampleBuffer via pointers as this whole inherited object might be modified by other
  threads. */
  ReadWriteLock audioDataReadWriteLock;


  juce_UseDebuggingNewOperator;
};

#endif  