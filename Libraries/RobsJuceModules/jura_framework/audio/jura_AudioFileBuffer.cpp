
//-------------------------------------------------------------------------------------------------
// construction/destruction:

AudioFileBuffer::AudioFileBuffer(const File &fileToLoadFrom) : AudioSampleBuffer(1, 0)
{
  ScopedWriteLock scopedLock(audioDataReadWriteLock);
  audioFileRootDirectory   = getApplicationDirectory();
  if( fileToLoadFrom.existsAsFile() )  // was formerly: if( fileToLoadFrom != File() )
    loadAudioDataFromFile(fileToLoadFrom);
}

AudioFileBuffer::~AudioFileBuffer() throw()
{
  ScopedWriteLock scopedLock(audioDataReadWriteLock);
  users.getLock().enter();

  // make sure, that no user still refers to this buffer after it's deleted - actually, ideally
  // this should already be the case when you delete an AudioFileBuffer object:
  jassert( users.size() == 0 );
  while( users.size() > 0 )
  {
    //AudioFileBufferUser* user = users[0];
    users[0]->assignAudioFileBuffer(nullptr);
    // ... or maybe we should call some sort of deletion notification callback?
  }

  users.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// inquiry:

bool AudioFileBuffer::isAudioFileValid() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return info.isValidAudioFile;
}

int AudioFileBuffer::getNumSamples() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return info.numSamples;
}

int AudioFileBuffer::getNumChannels() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return info.numChannels;
}

double AudioFileBuffer::getFileSampleRate() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return info.sampleRate;
}

double AudioFileBuffer::getLengthInSeconds() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return info.numSamples / info.sampleRate;
}

const File AudioFileBuffer::getUnderlyingFile() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return theAudioFile;
}

const String AudioFileBuffer::getFileName() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return getUnderlyingFile().getFileName();
}

const String AudioFileBuffer::getFileNameWithoutExtension() const
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  return getUnderlyingFile().getFileNameWithoutExtension();
}

void AudioFileBuffer::copyTo(const int destChannel, const int destStartSample,
                             AudioSampleBuffer &destBuffer, const int sourceChannel,
                             const int sourceStartSample, int numSamples)
{
  ScopedReadLock scopedLock(audioDataReadWriteLock);
  destBuffer.copyFrom(destChannel, destStartSample, *this, sourceChannel, sourceStartSample,
    numSamples);
}

//float* AudioFileBuffer::getSampleData(const int channelNumber, const int sampleOffset)
//{
//  return AudioSampleBuffer::getWritePointer(channelNumber, sampleOffset);
//}

//-------------------------------------------------------------------------------------------------
// thread synchronization:

void AudioFileBuffer::acquireReadLock()
{
  audioDataReadWriteLock.enterRead();
}

void AudioFileBuffer::releaseReadLock()
{
  audioDataReadWriteLock.exitRead();
}

void AudioFileBuffer::acquireWriteLock()
{
  audioDataReadWriteLock.enterWrite();
}

void AudioFileBuffer::releaseWriteLock()
{
  audioDataReadWriteLock.exitWrite();
}

//-------------------------------------------------------------------------------------------------
// others:

bool AudioFileBuffer::loadAudioDataFromFile(const File &fileToLoadFrom, bool showAlertBoxOnFail)
{
  ScopedWriteLock scopedLock(audioDataReadWriteLock);

  // retrieve meta-information about the file:
  info = AudioFileInfo(fileToLoadFrom);
  if( !info.isValidAudioFile )
  {
    if( showAlertBoxOnFail )
      showAudioFileInvalidErrorBox(fileToLoadFrom.getFileName());
    theAudioFile = File();
    return false;
  }

  // file is valid, read it in:
  AudioFormatManager formatManager;
  formatManager.registerBasicFormats();
  AudioFormatReader* reader = formatManager.createReaderFor(fileToLoadFrom);
  if( reader != NULL )
  {
    // we need some temporary int arrays because JUCE's AudioFormatReader reads the data only
    // into such arrays:
    int*  intDataBuffer = new int[(int) (reader->numChannels * reader->lengthInSamples)];
    if( intDataBuffer == NULL )
    {
      if( showAlertBoxOnFail )
        showMemoryAllocationErrorBox("AudioFileBuffer::loadAudioDataFromFile");
      delete reader;
      return false;
    }

    int** intChannelPointers = new int*[reader->numChannels+1];
    if( intChannelPointers == NULL )
    {
      if( showAlertBoxOnFail )
        showMemoryAllocationErrorBox("AudioFileBuffer::loadAudioDataFromFile");
      delete[] intDataBuffer;
      delete   reader;
      return false;
    }

    // assign the pointer-pointers to their targets, and zero-terminate the array of channel
    // pointers:
    unsigned int c;
    for(c=0; c<reader->numChannels; c++)
      intChannelPointers[c] = &(intDataBuffer[0 + c*reader->lengthInSamples]);
    intChannelPointers[reader->numChannels] = NULL;

    // read the integer data into the temporary arrays:
    //reader->read(intChannelPointers, 0, (int)reader->lengthInSamples); // old - for juce 1.46
    reader->read(intChannelPointers, reader->numChannels, 0, (int)reader->lengthInSamples, true);

    // now we allocate the memory for our actual member data, and convert the integer sample data
    // to float and store it in the inherited AudioSampleBuffer object
    AudioSampleBuffer::setSize(reader->numChannels, (int) reader->lengthInSamples);
    double normalizer  = 1.0/2147483648.0;  // = 1/(2^31) = 1/fullScale @ 32 bit int
    int n;
    float *writePointer;
    for(c=0; c<reader->numChannels; c++)
    {
      writePointer = AudioSampleBuffer::getWritePointer(c, 0);
      for(n=0; n<reader->lengthInSamples; n++)
        writePointer[n] = (float) (normalizer * (double) (intChannelPointers[c][n]));
    }

    delete[] intDataBuffer;
    delete[] intChannelPointers;
    delete   reader;
  }

  // update the aoociated file:
  theAudioFile = fileToLoadFrom;

  return true;
}

void AudioFileBuffer::filePathChanged(const File &newFilePath)
{
  ScopedWriteLock scopedLock(audioDataReadWriteLock);
  theAudioFile = newFilePath;
}

void AudioFileBuffer::registerUser(AudioFileBufferUser *userToRegister)
{
  //ScopedWriteLock scopedLock(audioDataReadWriteLock);
  users.getLock().enter();
  users.addIfNotAlreadyThere(userToRegister);
  users.getLock().exit();
}

void AudioFileBuffer::deRegisterUser(AudioFileBufferUser *userToDeRegister)
{
  //ScopedWriteLock scopedLock(audioDataReadWriteLock); // unnecesarry and causes dealock on dropping a clip onto the arrangement
  users.getLock().enter();
  users.removeFirstMatchingValue(userToDeRegister);
  users.getLock().exit();
}
