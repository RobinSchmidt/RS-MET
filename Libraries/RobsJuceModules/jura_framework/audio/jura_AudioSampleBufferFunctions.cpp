
void copyAudioSampleBufferContents(AudioSampleBuffer source, AudioSampleBuffer destination)
{
  destination.setSize(source.getNumChannels(), source.getNumSamples(), false, true, false);
  for(int c=0; c<source.getNumChannels(); c++)
    destination.copyFrom(c, 0, source, c, 0, source.getNumSamples());
}

AudioSampleBuffer* resampleBuffer(AudioSampleBuffer* inBuffer, double inSampleRate, 
  double outSampleRate)
{
  // calculate increment by which we step through the inBuffer:
  double increment = inSampleRate / outSampleRate;

  // obtain the lengths of the in- and out-buffers:
  int inLength    = inBuffer->getNumSamples();
  int numChannels = inBuffer->getNumChannels();
  int outLength   = (int) ceil( (double) inLength / increment );

  // create a new buffer with the desired target-length:
  AudioSampleBuffer* outBuffer = new AudioSampleBuffer(numChannels, outLength);
 
  // outer loop over the channels:
  for(int c=0; c<numChannels; c++)
  {
    // obtain pointers to the in- and out-buffers for the current channel
    const float* readBuf  = inBuffer->getReadPointer(  c, 0);
    float*       writeBuf = outBuffer->getWritePointer(c, 0);

    // initializations:
    double nRead     = 0.0;
    double nReadFrac = 0.0;
    int    nReadInt  = 0;
    double value     = 0.f;

    // inner loop over the samples:
    for(int nWrite=0; nWrite<outLength; nWrite++)
    {
      nReadInt  = (int) floor(nRead);
      nReadFrac = nRead - (double) nReadInt;

      // this crappy linear interpolation is only preliminary....implement windowed-sinc later
      jassert( nReadInt+1 < inLength );
      if( nReadInt+1 < inLength ) // should always be, actually
        value = (1.0-nReadFrac)*readBuf[nReadInt] + nReadFrac*readBuf[nReadInt+1];
      else
        value = 0.0;

      writeBuf[nWrite] = (float) value;

      nRead += increment;
    }

  }

  return outBuffer;
}

void saveAudioSampleBufferToFile(AudioSampleBuffer* bufferToSave, File fileToSaveTo,                                 
                                 double sampleRateToUse, int bitsPerSample, 
                                 bool warnIfFileExists)
{
  if( warnIfFileExists )
  {
    if( fileToSaveTo.existsAsFile() )
    {
      if( showOverwriteAudioFileWarningBox(fileToSaveTo.getFileName()) == true )
      {
        // user has confirmed to overwrite - we clear the existing file:
        fileToSaveTo.deleteFile();
      }
      else
        return;
    }
  }
  else
  {
    if( fileToSaveTo.existsAsFile() )
      fileToSaveTo.deleteFile();
  }

  AudioFormat *audioFormat = NULL;
  if( fileToSaveTo.hasFileExtension(String("wav")) )
    audioFormat = new WavAudioFormat();
  else if( fileToSaveTo.hasFileExtension(String("flac")) )
    audioFormat = new FlacAudioFormat();
  else if( fileToSaveTo.hasFileExtension(String("aiff")) )
    audioFormat = new AiffAudioFormat();
  else
    return;

  // create a FileOutputStream to write into:
  //FileOutputStream* outputStream = fileToSaveTo.createOutputStream();  // old
  FileOutputStream* outputStream = fileToSaveTo.createOutputStream().get();  // new

  // create a writer for the stream:
  AudioFormatWriter *writer = audioFormat->createWriterFor(outputStream, sampleRateToUse, 
    bufferToSave->getNumChannels(), bitsPerSample, StringPairArray(), 0);

  // if we have a writer, use it, otherwise clean up and shown an error message box:
  if( writer != NULL )
  {
    // create a temporary instance of ImmediatePlaybackAudioSource for the writer to use:
    ImmediatePlaybackAudioSource source;
    source.startPlayback(bufferToSave);

    // write:
    writer->writeFromAudioSource(source, bufferToSave->getNumSamples());

    // delete the writer:
    delete writer; // the destructor will also delete the outputStream
  }
  else
  {
    // as there was no writer to take care of the deletion of the the outputStream, we have to do 
    // it ourselves:
    delete outputStream;

    // show an error reporting message box:
    showAudioFileWriteErrorBox(fileToSaveTo.getFullPathName());
  }

  // clear the created AudioFormat object:
  delete audioFormat;
}