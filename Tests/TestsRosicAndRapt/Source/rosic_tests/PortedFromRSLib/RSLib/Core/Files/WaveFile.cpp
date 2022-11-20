// move them into the class:
const static char expectedChunkID[]          = "RIFF";
const static char expectedFormatDescriptor[] = "WAVE";
const static char expectedSubChunk1ID[]      = "fmt ";
const static char expectedSubChunk2ID[]      = "data";

//-------------------------------------------------------------------------------------------------
// class rsInputWaveFile:

bool rsInputWaveFile::readHeaderAndResetPosition()
{
  fseek(filePointer, 0, SEEK_SET);
  memset(&header, 0, sizeof(header));
  positionInData = 0;
  if( fread(&(header), sizeof(WaveFileHeader), 1, filePointer) != 1 )
    return false;
  else
    return true;
}

bool rsInputWaveFile::isFileFormatSupported()
{
  bool result = true;

  // check RIFF section of header:
  result &= (memcmp(expectedChunkID, header.riffChunkDescriptor.chunkID, 4) == 0);
  result &= (memcmp(expectedFormatDescriptor, header.riffChunkDescriptor.formatDescriptor, 4) 
             == 0);

  // check 'fmt'-subchunk section of header:
  result &= (memcmp(expectedSubChunk1ID, header.formatSubChunk.subChunk1ID, 4) == 0);
  result &= header.formatSubChunk.sampleFormat  == PCM_LINEAR;
  result &= header.formatSubChunk.bitsPerSample == 16;

  // check 'data'-subchunk section of header:
  result &= (memcmp(expectedSubChunk2ID, header.dataSubChunk.subChunk2ID, 4) == 0);

  return result;
}

bool rsInputWaveFile::openForRead()
{
  bool success = rsFileStream::openForRead();
  success &= readHeaderAndResetPosition();
  success &= isFileFormatSupported();
  if( !success )
    rsError("File format not supported");
    // should we close the file in this case?
  return success;
}

int rsInputWaveFile::read16BitInt(rsInt16 *buffer, int maxElems)  
{
  // \todo rename maxElemes into something more descriptive

  rsUint32 afterDataRead;  // maybe rename to newPosition
  int numBytes;
  int numElems;

  numBytes = maxElems * 2;
  afterDataRead = positionInData + numBytes;
  if(afterDataRead > header.dataSubChunk.totalNumBytesOfData)
  {
    // Don't read more samples than are marked available in header
    numBytes = header.dataSubChunk.totalNumBytesOfData - positionInData;
    //assert(numBytes >= 0);
  }

  numBytes        = (int)fread(buffer, 1, numBytes, filePointer);
  positionInData += numBytes;
  numElems        = numBytes / 2;

  // \todo: if the end of the wavefile is reached, fill the rest of the buffer with zeros

  return numElems;
}

int rsInputWaveFile::readAndConvertToFloat(float *buffer, int maxElems)
{
  if( header.formatSubChunk.bitsPerSample == 16 )
  {
    rsInt16 *shortBuffer = new rsInt16[maxElems];
    int   numRead      = read16BitInt(shortBuffer, maxElems);
    convert16BitToFloat(shortBuffer, buffer, numRead);
    delete[] shortBuffer;
    return numRead;
  }
  else
    return 0;
}

void rsInputWaveFile::convert16BitToFloat(rsInt16 *inBuffer, float *outBuffer, int length)
{
  const double scaler = 1.0 / 32768.0;  // should perhaps use 32767
  for(int i = 0; i < length; i++)
    outBuffer[i] = (float) (scaler * (double) inBuffer[i]);
}

int rsInputWaveFile::isEndOfFileReached() const
{
  return( positionInData == header.dataSubChunk.totalNumBytesOfData || feof(filePointer) );
}

//-------------------------------------------------------------------------------------------------
// class OutputWaveFile:

rsOutputWaveFile::rsOutputWaveFile(const rsString& absolutePath, int sampleRate, 
  int bitsPerSample, int numChannels, int sampleFormat) :
rsWaveFile(absolutePath)
{
  openForWrite();
  positionInData = 0;
  createPreliminaryHeader(sampleRate, bitsPerSample, numChannels, sampleFormat);
  writeHeader();
}

void rsOutputWaveFile::createPreliminaryHeader(const rsUint32 sampleRate, const rsUint32 bits, 
  const rsUint32 channels, const rsUint32 sampleFormat)
{
  // fill in the 'riff' part:
  memcpy(&(header.riffChunkDescriptor.chunkID), expectedChunkID, 4);   // chunkID = 'RIFF'
  header.riffChunkDescriptor.chunkSize = 0;                            // chunkSize unknown so far
  memcpy(&(header.riffChunkDescriptor.formatDescriptor), expectedFormatDescriptor, 4);  // 'WAVE'

  // fill in the 'format' part:
  memcpy(&(header.formatSubChunk.subChunk1ID), expectedSubChunk1ID, 4); // subChunk1ID = 'fmt '
  header.formatSubChunk.subChunk1Size       = 0x10;
  //header.formatSubChunk.sampleFormat        = PCM_LINEAR;
  header.formatSubChunk.sampleFormat        = sampleFormat;
  header.formatSubChunk.numChannels         = (short) channels;
  header.formatSubChunk.sampleRate          = sampleRate;
  header.formatSubChunk.bitsPerSample       = (short) bits;
  header.formatSubChunk.bytesPerSampleFrame = (short) (bits * channels / 8);
  header.formatSubChunk.bytesPerSecond      = header.formatSubChunk.bytesPerSampleFrame*sampleRate;
  header.formatSubChunk.sampleRate          = sampleRate;

  // fill in the 'data' part:
  memcpy(&(header.dataSubChunk.subChunk2ID), expectedSubChunk2ID, 4);  // copy 'data' to data_field
  header.dataSubChunk.totalNumBytesOfData = 0;                         // data_len unknown so far
}

void rsOutputWaveFile::finalizeHeaderAndCloseFile()
{
  header.riffChunkDescriptor.chunkSize    = positionInData + 36;
  header.dataSubChunk.totalNumBytesOfData = positionInData;
  writeHeader();
  close();
}

void rsOutputWaveFile::writeHeader()
{
  fseek(filePointer, 0, SEEK_SET);
  int success = (int)fwrite(&header, sizeof(header), 1, filePointer);
  if(success != 1)
    rsError("Error while writing to a wav file.");
  fseek(filePointer, 0, SEEK_END);
}

void rsOutputWaveFile::write16Bit(const short *buffer, int numElems)
{
  // should we assert that numChannles == 1?
  if(numElems < 1) return; // nothing to do
  rsAssert(header.formatSubChunk.bitsPerSample == 16);
  int res = (int)fwrite(buffer, 2, numElems, filePointer);
  if(res != numElems)
    rsError("Error while writing to a wav file.");
  positionInData += 2 * numElems;
}

void rsOutputWaveFile::write24Bit(const rsInt32* buffer, int numElems)
{
  // should we assert that numChannles == 1?
  if(numElems < 1) return; // nothing to do
  rsAssert(header.formatSubChunk.bitsPerSample == 24);
  rsInt8* buf8 = new rsInt8[3*numElems];
  copyBytes4to3((const rsInt8*) buffer, 4*numElems, buf8);
  int res = (int)fwrite(buf8, 3, numElems, filePointer);
  if(res != numElems)
    rsError("Error while writing to a wav file.");
  positionInData += 3 * numElems;
  delete[] buf8;
}

void rsOutputWaveFile::write(const float *buffer, int numElems)
{
  // should we assert that numChannles == 1?
  if(getBitsPerSample() == 16)
  {
    short* shortBuffer = new short[numElems];
    convertFloatTo16BitInt(buffer, shortBuffer, numElems);
    write16Bit(shortBuffer, numElems);
    delete[] shortBuffer;
  }
  else if(getBitsPerSample() == 24)
  {
    rsInt32* buf32 = new rsInt32[numElems];
    convertFloatTo24BitInt(buffer, buf32, numElems);
    write24Bit(buf32, numElems);
    delete[] buf32;
  }
  else if(getBitsPerSample() == 32)
  {
    if(getSampleFormat() == SampleFormat::PCM_LINEAR)
    {
      // todo:
      //rsInt32* buf32 = new rsInt32[numElems];
      //convertFloatTo32BitInt(buffer, buf32, numElems);
      //write32Bit(buf32, numElems);
      //delete[] buf32;
    }
    if(getSampleFormat() == SampleFormat::IEEE_FLOAT)
    {
      int res = (int)fwrite(buffer, 4, numElems, filePointer);
      if(res != numElems)
        rsError("Error while writing to a wav file.");
      positionInData += 4 * numElems;
    }
  }
  else if(getBitsPerSample() == 64 && getSampleFormat() == SampleFormat::IEEE_FLOAT)
  {

  }
  else
    rsError("Required bit format not supported.");

  // ToDo: 
  // -implement 8,32 bit integer and 64 bit floating point formats
  // -implement handling of multichannel data...maybe this function here can be used for that, but
  //  we need to interleave the data - maybe the caller should do that and here, we may only 
  //  provide a convenience function for that ...maybe it should not allocate memeory but 
  //  interleave -> write -> de-interleave everything in place, maybe the function should have a
  //  boolean parameter isInterleaved
  // -maybe implement A-law and mu-law
  //  https://web.archive.org/web/20110719132013/http://hazelware.luggle.com/tutorials/mulawcompression.html
  //  ...but no - that's a lossy compression that only applies to 16 bit data, so we don't really
  //  want to encode files like this - but it may be the case that we occasionally need to decode
  //  them, so maybe it should be integrated for completeness
}

void rsOutputWaveFile::write(const double* buffer, int numElems)
{
  if(getBitsPerSample() == 64 && getSampleFormat() == SampleFormat::IEEE_FLOAT)
  {
    int res = (int)fwrite(buffer, 8, numElems, filePointer);
    if(res != numElems)
      rsError("Error while writing to a wav file.");
    positionInData += 8 * numElems;
  }
  else
  {
    float* floatBuffer = new float[numElems];
    for(int i = 0; i < numElems; i++)
      floatBuffer[i] = (float)buffer[i];
    write(floatBuffer, numElems);
    delete[] floatBuffer;
    // It's silly to first convert double to float and then float to short or whatever - convert 
    // directly to target format...but that may introduce code duplication...maybe templatize the
    // relevant method, so it can be called for double and float
  }
}

void rsOutputWaveFile::convertFloatTo16BitInt(
  const float *inBuffer, rsInt16 *outBuffer, int length)
{
  for(int i = 0; i < length; i++)
  {


    //int tmp = (int) (32768.f * inBuffer[i]); // or should the factor be 32767.f?
    //int tmp = (((int) ((inBuffer[i] * 32767) + 32768.5f)) - 32768);
    //outBuffer[i] = (short) rsLimitToRange(tmp, -32768, 32767);

    outBuffer[i] = float32ToInt16(inBuffer[i]);
    // so, to get the float back, we would need to do:
    // floatVal = (((float)(intVal+32768)) - 32768.5f) / 32767.f;
    // -> test, if this works for an int -> float -> int roundtrip simulating a load/save roundtrip
    // we can actually do an exhaustive test for all possible 2^16 values here, for 24 bit, an 
    // exhaustive test may be practical to do rarely
  }
  // ToDo: use round, see
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=552377&start=45
  // or use:
  // (((int) ((x * 32767) + 32768.5f)) - 32768)
  // as mystran suggests, the formula is also here:
  // https://www.cs.cmu.edu/~rbd/papers/cmj-float-to-int.html

  // Using the factor 32768 is a power of 2 so it has the advantage of not introducing any 
  // additional rounding error - but if the maximum sample in a normalized buffer happens to be
  // negative, it will be clipped. Maybe the normalization should take into account the question,
  // whether the max-sample is negative and use different factors in the different cases. Look at 
  // how other wav-writing libraries handle this.
  //
  // See:
  // https://github.com/libsndfile/libsndfile
  // https://github.com/libsndfile/libsndfile/blob/master/src/wav.c
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=520701&p=7323214
  //
}

void rsOutputWaveFile::convertFloatTo24BitInt(
  const float* inBuffer, rsInt32* outBuffer, int length)
{
  for(int i = 0; i < length; i++)
  {
    outBuffer[i] = float32ToInt24(inBuffer[i]);

    ////int tmp = (int) (8388608.f * inBuffer[i]);    // 8388608 = 2^23
    //int tmp = (((int) ((inBuffer[i] * 8388607) + 8388608.5f)) - 8388608);
    //outBuffer[i] = rsLimitToRange(tmp, -8388608, 8388607);
  }
}

int rsOutputWaveFile::copyBytes4to3(const rsInt8* x, int N, rsInt8* y)
{
  rsAssert(N % 4 == 0, "N must be divisible by 4");
  N /= 4;          // N is now the number of samples, not the number of bytes anymore
  int ix = 0;      // index into x array
  int iy = 0;      // index into y array
  for(int n = 0; n < N; n++) {
    y[iy+0] = x[ix+0];
    y[iy+1] = x[ix+1];
    y[iy+2] = x[ix+2];
    ix += 4;
    iy += 3; }
  return 3*N;      // 3*N is the number of bytes written into y

  // Notes:
  // -I think, initializing ix = 0 works only on little endian machines. On big endian machines, we
  //  may need ix = 1 or ix = 3 -> try that. Maybe we should have a conditional compilation like 
  //  #ifdef RS_BIG_ENDIAN ix = 1 or something
  // -Algorithm can be generalized to copyElementsMtoN that takes as parameters input and output 
  //  spacings M,N and initial offsets for ix, iy (deault 0) -> maybe implement that in a general 
  //  way in rsArrayTools. I think, it can be used in place, if N < M (verify!). Or maybe it should
  //  be named copyWithGaps, copyGapped(T* x, int sx, T* y, int sy, int N, int ix = 0, iy = 0)
  //  sx,,sy are the spacings, ix, iy the initial offsets/indices (0 <= ix < sx, ...)
}

/*

ToDo:
-don't use char, short, long - instead use rsInt8, rsInt16, rsInt32 consistently
-use members of type std::vector for the buffers, resize as necessary -> avoid unecessary 
 allocations
-write a unit test 
 -test int -> float -> int roundtrips between float32 and 16 and 24 bit integer data
-avoid memory allocations in each call to write* -> have std::vector members for various buffers 
 and resize them in write* as necessary
-when saving the same file with 16 and 24 bit then converting 24 to 16 and subtracting it from
 the 16bit file, the result is not zero - there's a noise at -90dB (at least when doing it in 
 Audition) ...is that expected? -> do experiments and unit tests (maybe just for the conversions, 
 without necessarily writing to files), also implement dithering ..but not in this class
 https://en.wikipedia.org/wiki/Dither
 https://en.wikipedia.org/wiki/Noise_shaping
 https://dsp.stackexchange.com/questions/15436/audio-signal-dither-and-noise-shaping
 https://ieeexplore.ieee.org/document/969564

 LEAST SQUARES THEORY AND DESIGN OF OPTIMAL NOISE SHAPING FILTERS:
 https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.5.6708&rep=rep1&type=pdf


References:
https://en.wikipedia.org/wiki/WAV
https://www.aelius.com/njh/wavemetatools/doc/riffmci.pdf
http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
https://sites.google.com/site/musicgapi/technical-documents/wav-file-format



See also:
Write a WAV file from scratch - C++ Audio Programming https://www.youtube.com/watch?v=qqjvB_VxMRM&list=PLbqhA-NKGP6Cwrh-FUMFpc_pOQbardIRh


*/