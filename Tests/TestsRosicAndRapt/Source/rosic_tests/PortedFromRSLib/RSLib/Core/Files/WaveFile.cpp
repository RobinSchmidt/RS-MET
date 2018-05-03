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
  const double scaler = 1.0 / 32768.0;
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
                                   int bitsPerSample, int numChannels) :
rsWaveFile(absolutePath)
{
  openForWrite();
  positionInData = 0;
  createPreliminaryHeader(sampleRate, bitsPerSample, numChannels);
  writeHeader();
}

void rsOutputWaveFile::createPreliminaryHeader(rsUint32 sampleRate, rsUint32 bits, 
                                               rsUint32 channels)
{
  // fill in the 'riff' part:
  memcpy(&(header.riffChunkDescriptor.chunkID), expectedChunkID, 4);   // chunkID = 'RIFF'
  header.riffChunkDescriptor.chunkSize = 0;                            // chunkSize unknown so far
  memcpy(&(header.riffChunkDescriptor.formatDescriptor), expectedFormatDescriptor, 4);  // 'WAVE'

  // fill in the 'format' part:
  memcpy(&(header.formatSubChunk.subChunk1ID), expectedSubChunk1ID, 4); // subChunk1ID = 'fmt '
  header.formatSubChunk.subChunk1Size       = 0x10;
  header.formatSubChunk.sampleFormat        = PCM_LINEAR;
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
  if(numElems < 1) return; // nothing to do
  //assert(header.formatSubChunk.bitsPerSample == 16);
  int res = (int)fwrite(buffer, 2, numElems, filePointer);
  if(res != numElems)
    rsError("Error while writing to a wav file.");
  positionInData += 2 * numElems;
}

void rsOutputWaveFile::write(const float *buffer, int numElems)
{
  short *shortBuffer = new short[numElems];
  convertFloatTo16BitInt(buffer, shortBuffer, numElems);
  write16Bit(shortBuffer, numElems);
  delete[] shortBuffer;
}

void rsOutputWaveFile::convertFloatTo16BitInt(const float *inBuffer, rsInt16 *outBuffer, 
                                              int length)
{
  int tmp;
  for(int i=0; i<length; i++)
  {
    tmp          = (int) (32768.0f * inBuffer[i]);
    outBuffer[i] = (short) rsLimitToRange(tmp, -32768, 32767);
  }
}
