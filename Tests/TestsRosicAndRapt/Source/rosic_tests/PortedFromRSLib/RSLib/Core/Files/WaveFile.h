#ifndef RS_WAVEFILE_H
#define RS_WAVEFILE_H

namespace RSLib
{

  /** Baseclass for the WaveInputFile and WaveOutputFile classes. This is just a dumb data class 
  that holds a wavefile header and a pointer to the actual data and some type definitions.

  \todo maybe absorb the reader and writer class into this single class and provide functions 
  like openForWrite, etc.

   */

  class RSLib_API rsWaveFile : public rsFileStream
  {

  public:

    rsWaveFile(const rsString& absolutePath = rsString::empty) : rsFileStream(absolutePath) {}

     
    /** \name Inquiry */

    int getBitsPerSample() const { return header.formatSubChunk.bitsPerSample; }

    enum sampleFormats
    {
    
      PCM_LINEAR = 1      // ...other formats are currently not supported
      // todo:
      //IEEE_FLOAT = 3,
      //ALAW = 6,
      //MULAW = 7
    };
    // see http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html

    int getSampleFormat() const { return header.formatSubChunk.sampleFormat; }
    // one of the values from enum sampleFormats


  protected:



    /** The "RIFF" chunk descriptor part of the header. */
    typedef struct
    {
      char chunkID[4];              // the ASCII letters "RIFF"
      int  chunkSize;               // size of the rest of the entire file
      char formatDescriptor[4];     // the ASCII letters "WAVE"
    } RiffChunkDescriptor;

    /** The "fmt" sub-chunk part of the header. */
    typedef struct
    {
      char  subChunk1ID[4];         // the ASCII letters "fmt_"
      int   subChunk1Size;          // rest of the size of the subchunk, 16 for PCM
      short sampleFormat;           // 1: linear PCM, 3: float, other values: nonlinear quantization
      short numChannels;            // mono: 1, stereo: 2, etc.
      int   sampleRate;             // 8000, 44100, 96000, etc.
      int   bytesPerSecond;         // == sampleRate * numChannels * bitsPerSample/8
      short bytesPerSampleFrame;    // == numChannels * bitsPerSample/8
      short bitsPerSample;          // 8, 16, etc.
    } FormatSubChunk;
    // todo: support also 24, 32 bit int and 32, 64 bit float

    /** The data sub-chunk minus the actual data. */
    typedef struct
    {
      char   subChunk2ID[4];         // the ASCII letters "data"
      rsUint32 totalNumBytesOfData;    // == numSamples * numChannels * bitsPerSample/8
    } DataSubChunk;

    /** The complete header consisiting of the 3 sub-headers defined above. */
    typedef struct
    {
      RiffChunkDescriptor riffChunkDescriptor;
      FormatSubChunk      formatSubChunk;
      DataSubChunk        dataSubChunk;
    } WaveFileHeader;

    // data members:
    WaveFileHeader header;     // holds the header as struct
    rsUint32 positionInData;   // current position inside the actual raw data (in bytes)

  };

  //===============================================================================================

  /** Class for reading wavefiles. Currently, the use of this class is limited to 16-bit files on 
  little-endian machines. */

  class RSLib_API rsInputWaveFile : public rsWaveFile
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Opens the file for reading. */
    rsInputWaveFile(const rsString& absolutePath = rsString::empty) 
      : rsWaveFile(absolutePath) 
    { 
      openForRead(); 
    }


    /** \name Inquiry */

    virtual rsUint32 getSampleRate()       const { return header.formatSubChunk.sampleRate; }
    virtual rsUint32 getNumBitsPerSample() const { return header.formatSubChunk.bitsPerSample; }
    virtual rsUint32 getDataSizeInBytes()  const { return header.dataSubChunk.totalNumBytesOfData;}
    virtual rsUint32 getNumChannels()      const { return header.formatSubChunk.numChannels; }
    virtual rsUint32 getNumSampleFrames()  const 
    { 
      return getDataSizeInBytes() / header.formatSubChunk.bytesPerSampleFrame; 
    }

    /** Checks whether the end of the file is reached. */
    virtual int isEndOfFileReached() const;


    /** \name Reading */

    /** Initializes the file-stream for reading and reports if this was successful. */
    virtual bool openForRead();

    /** Starting from the current position, this function reads the given number of sample-frames 
    from the file and stores them in the passed array. If the file does not contain that many 
    samples anymore, it will read to the end of file. The return value informs how many 
    sample-frames were read. */
    virtual int read16BitInt(rsInt16 *buffer, int maxElems);

    /** Similar to read16Bit, but reads into a float buffer therby doing the conversion. */
    virtual int readAndConvertToFloat(float *buffer, int maxElems);

  protected:

    /** \name Internal Functions */

    /** Converts a buffer of 16 bit integer samples in the range [-32768, 32767] to floating point 
    numbers in the range [-1..+1[. */
    virtual void convert16BitToFloat(rsInt16 *inBuffer, float *outBuffer, int length);

    /** Returns true when we can read the given file format, false otherwise (assumes that we 
    already have read the header from the file). */
    virtual bool isFileFormatSupported();

    /** Reads the header data-structure from the file and reports if this was successful. */
    virtual bool readHeaderAndResetPosition();

  };

  //===============================================================================================

  /** Class for writing wavefiles. Currently, the use of this class is limited to 16-bit files on 
  little-endian machines. */

  class RSLib_API rsOutputWaveFile : public rsWaveFile
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Opens the file for writing. */
    rsOutputWaveFile(const rsString& absolutePath, int sampleRate, int bitsPerSample, 
      int numChannels);

    /** Destructor. Finalizes the header and closes the file. */
    virtual ~rsOutputWaveFile() { finalizeHeaderAndCloseFile(); }





    /** \name Writing */

    /** Writes 16 bit integer data into the file. */
    virtual void write16Bit(const short *buffer, int numElems);
    // use rsInt16

    /** Writes 24 bit integer data into the file. Assumes that only the lowest 3 bytes of the 4 
    byte integers are used. Note that on little endian systems, "lowest" means "first", but the 
    caller actually does not need to know or care about this. */
    virtual void write24Bit(const rsInt32 *buffer, int numElems);

    /** Writes float data into the file, thereby converting the format if necessary. */
    virtual void write(const float *buffer, int numElems);

    virtual void write(const double *buffer, int numElems);

    /** Supplements the missing information about the filesize to the preliminary header and 
    closes the file. */
    virtual void finalizeHeaderAndCloseFile();



  protected:

    /** \name Internal Functions */


    /** Creates the preliminray header which does not yet contain valid data for the file 
    length. */
    virtual void createPreliminaryHeader(const rsUint32 sampleRate, const rsUint32 bitsPerSample, 
      const rsUint32 numChannels);

    /** Writes the header into the file. */
    virtual void writeHeader();

    /** Converts a buffer of floating point numbers in the range [-1..+1] to 16 bit integer 
    samples in the range [-32768, 32767]. */
    virtual void convertFloatTo16BitInt(const float *inBuffer, rsInt16 *outBuffer, int length);

    virtual void convertFloatTo24BitInt(const float *inBuffer, rsInt32 *outBuffer, int length);

    /** Given a buffer x of N bytes in groups of 4, this copies the content into a buffer of bytes 
    in groups of 3, leaving out the most significant byte in the source buffer. This is used for
    converting 32 integers into 24 bit integers represented by 3 bytes each. The return value is 
    the number of bytes written, which is 3*N/4. N must be divisble by 4. */
    virtual int copyBytes4to3(const rsInt8* x, int N, rsInt8* y);
 

  };

}

#endif
