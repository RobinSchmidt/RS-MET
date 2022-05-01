#ifndef rosic_AudioStream_h
#define rosic_AudioStream_h

namespace rosic
{

template<class T>   // maybe it should go into rapt
class AudioStream 
{

public:

  virtual ~AudioStream() {}

  /** For random access. Writes the sample frame with given index into the given destination. */
  virtual void getFrame(int sampleIndex, T* destination) const = 0;

  /** Function for speficically handling stereo signals to allow handling that important, common
  special case more efficiently than with the more general implementation. */
  virtual void getFrameStereo(int sampleIndex, T* left, T* right) const = 0;


  /** Computes a stereo output frame at possibly non-integer sample positions. The baseclass 
  implementation does this by simple linear interpolation. Subclasses may override this for using
  better quality interpolation and/or optimization purposes. If the samplePosition is outside the 
  valid range 0 <= samplePosition <= N-1, the signal there is taken to be zero. At positions in
  between -1..0 and N-1..N, linear interpolation will take place between the first (or last) actual
  sample and zero. */
  virtual void getFrameStereo(double samplePosition, T* left, T* right) const
  {
    int i = (int) samplePosition;         // integer part
    T   f = T(samplePosition - T(i));     // fractional part

    T xL0(0), xR0(0), xL1(0), xR1(0);     
    // Is this good as default? Maybe we should use the very first and very last sample here. This
    // may work better for looping.

    if(isValidFrameIndex(i  ))  getFrameStereo(i,   &xL0, &xR0);
    if(isValidFrameIndex(i+1))  getFrameStereo(i+1, &xL1, &xR1);

    *left  = (T(1)-f) * xL0 + f * xL1;
    *right = (T(1)-f) * xR0 + f * xR1;

    // ToDo:
    // -Make a function that takes sampleIndex and frac as separate arguments
    // -Make an unsafe version that avoids the isValidFrameIndex checks
    // -Write a function that fills a whole buffer instead of producing one sample at a time. That
    //  function should figure out safe bounds once and make use of the unsafe version per sample.
    // -Try to optimize using modf: https://en.cppreference.com/w/cpp/numeric/math/modf
  }



  /** Sets the number of output channels for this object. By default, this number will be equal to
  the number of channels in the data, as set by the setData call. However, the number of desired
  output channels to be produced may actually be different from that. For example, we may want to
  produce stereo output from mono data by just copying the single output into both channels. Let
  M be the number of channels in the data and N be the number of output channels of the stream,
  then we will produce: out[n] = data[m % N] where n is the output channel index and m % N the data
  channel index. */
  //virtual void setNumOutputChannels(int newNumChannels) = 0;

  /** Returns the number of sample frames in this stream. */
  int getNumFrames() const { return numFrames; }


  double getSampleRate() const { return sampleRate; }

  bool isValidFrameIndex(int i) const { return i >= 0 && i < numFrames; }

  //virtual bool hasFinished(int sampleTime)

  /** Subclasses may want to override this for optimizing the blockwise access. */
  /*
  virtual void getBlock(int startIndex, int length, TSmp** destination)
  {
  for(int i = 0; i < length; i++)
  getFrame(startIndex + i, destination[i]);
  // assumes interleaved storage in memory
  }
  */

  void clear() { sampleRate  = 1; numChannels = 0; numFrames   = 0; }

protected:

  double sampleRate  = T(44100);  // or maybe init to 0? or 1? or 666?
  int numChannels = 0;
  int numFrames   = 0;         // maybe use -1 to encode "unknown"? would that be useful?

};
// maybe move out of this class - this may be useful in other contexts, too - maybe templatize


/*
class rsAudioFileInfo  // maybe rename to rsFileInfo for use in more general context
{

public:

  bool operator==(const rsAudioFileInfo& rhs) const
  {
    bool same = true;
    same &= fileName     == rhs.fileName;
    same &= extension    == rhs.extension;
    same &= path         == rhs.path;
    same &= rootDirIndex == rhs.rootDirIndex;
    return same;
  }

protected:

  std::string fileName;  // without filename extension (e.g. Piano_A4)
  std::string extension; // filename extension (e.g. wav, flac)
  std::string path;      // relative path from a predefined root directory
  int rootDirIndex = 0;  // index of root directory (among a couple of predefined choices)
  // rootDirIndex stores the root-directory to which the path is relative, but just as an integer
  // index that selects between various pre-defined root-directories that should exist at the
  // rsSamplerEngine level. For example: 0: instrument directory (containing e.g. Piano.sfz),
  // 1: factory sample directory, 2: user sample directory, etc. (but not hardcoded - meanings of
  // the indices should be flexible). This makes it more reasonably possible to uniquely identify
  // samples. It's totally possible to have samples in an instrument with same relative paths and
  // filenames but with respect to different root directories. Yes - that would be weird, but the
  // engine should neverless be able to handle such situations.
};
*/

//=================================================================================================

template<class T>
class AudioFileStream : public AudioStream<T>
{

public:

  // ToDo: add == operator based on fileName, extension, path (and maybe (meta)data such as 
  // sampleRate, numChannels, numFrames, ...)

  // inquiry: existsOnDisk (idea: we should allow to work with samples that are not stored on 
  // disk but rather pre-rendered programmatically into RAM at runtime)

  //virtual void getFrame(int sampleIndex, TSmp* destination) {}

  // todo: getBlock

  /** Returns true, iff the given other stream object uses the same audio file as this. */
  bool usesSameFile(const AudioFileStream<T>* otherStream) const
  {
    //return fileInfo == otherStream->fileInfo;

    bool same = true;
    //same &= fileName     == otherStream->fileName;
    //same &= extension    == otherStream->extension;
    same &= path         == otherStream->path;
    same &= rootDirIndex == otherStream->rootDirIndex;
    return same;
  }

  const std::string& getPath() const { return path; }
  //std::string getPath() const { return path + fileName + extension; }
  // todo: 
  // -don't split the path into 3 fields -> have only a path field which is the full path
  // -return a const reference to that here

  void clear() { AudioStream<T>::clear(); path.clear(); rootDirIndex = 0; }

protected:

  //rsAudioFileInfo fileInfo;


   // factor out into class rsFileInfo
  //std::string fileName;  // without filename extension (e.g. Piano_A4)
  //std::string extension; // filename extension (e.g. wav, flac)
  std::string path;      // relative path from a predefined root directory

  int rootDirIndex = 0;  // index of root directory (among a couple of predefined choices)
  // rootDirIndex stores the root-directory to which the path is relative, but just as an integer
  // index that selects between various pre-defined root-directories that should exist at the 
  // rsSamplerEngine level. For example: 0: instrument directory (containing e.g. Piano.sfz), 
  // 1: factory sample directory, 2: user sample directory, etc. (but not hardcoded - meanings of
  // the indices should be flexible). This makes it more reasonably possible to uniquely identify 
  // samples. It's totally possible to have samples in an instrument with same relative paths and 
  // filenames but with respect to different root directories. Yes - that would be weird, but the 
  // engine should neverless be able to handle such situations.
  // Actually, i don't like the idea that this class should be coupled to the sampler engine. It
  // should be usable in a much more general context. Maybe get rid of the rootDirIndex and replace
  // it by an actual string representing the root dir

  // todo: maybe keep just the path which should represents the full relative path. It doesn't seem
  // to be a good idea to split it into 3 parts
};

//=================================================================================================

// maybe rename to AudioFileStreamRAM, another subclass can be named AudioFileStreamDFD,
// or ...StreamFromMemory/FromDisk

template<class T>
class AudioFileStreamPreloaded : public AudioFileStream<T>
{

public:


  virtual ~AudioFileStreamPreloaded() { clear(); }


  /** Sets the new sample data. The function will *not* make the object take ownership over the 
  passed data nor will it refer to it later. Instead, it will be copy everything into an internal
  buffer. The numFrames parameter gives the number of sample frames, numDataChannels gives the 
  number of channels in the passed input data and numStreamChannels gives the number of desired
  channels in the stream. If these two numbers are different. It beahves as follows: If the input 
  has more channels than the stream, the strem will just use the first few of the input channels. 
  If the stream has more channels than the input, then the input channels will be mapped to the 
  outputs by using wrapping around the channelNumber as needed (i.e. using modulo). The important
  special cas here is when the input data is mono - the same data will just be broadcasted to all
  output channels. Returns true, iff everything went alright and false if it failed to allocate the
  required memory. */
  bool setData(T** newData, int numFrames, int numDataChannels, double sampleRate,
    int numStreamChannels, const std::string& path);
  // todo: 
  // -we need to distiguish between the number of channels in the data and the desired number of 
  //  output channels (done?)
  // -maybe use double for the sampleRate
  // -include fileName etc. later, too
  // -Maybe have variations of this function taking ownership or just refering to data owned 
  //  elsewhere. Give the different versions appropriate names: setDataByCopying, 
  //  setDataByReferring, setDataByTakingOwnership. We'll need a boolean flag to store in which 
  //  mode the data was passed so we know what the appropriate thing to do is when we get new data
  //  or the object gets destroyed.
  // -Maybe the numChannels parameter should come before numFrames because it's the first index 
  //  into the 2D array? -> check conventions used elsewhere. It would be consistent with rsMatrix
  //  when the rows are interpreted as channels.

  //void setNumOutputChannels(int newNumChannels) override;


  void clear();


  void getFrame(int sampleIndex, T* destination) const override
  {
    int n = sampleIndex;
    //RAPT::rsAssert(n >= 0 && n < numFrames, "sampleIndex out of range");
    for(int c = 0; c < this->numChannels; c++)
      destination[c] = channelPointers[c][n];
    // What, if the number of output channels shall be different than the number of of channels
    // in the data? Maybe it's best (most efficient) to ensure that this cannot happen on a 
    // higher level. We can just let our channelPointers all point to the same actual data, for
    // example. When the number of output channels of the sampler engine is changed (which 
    // should happen rarely, if ever), we need to update all channelPointer arrays in all 
    // AudioFileStreamPreloaded objects.
  }
  // this api sucks! pass a float** destinations - but for the sfz player, this is too general
  // we really need to support only mono and stereo. Maybe make a function getFrameStereo


  void getFrameStereo(int sampleIndex, T* left, T* right) const override
  {
    //RAPT::rsAssert(numChannels == 2); // Can be used only for stereo signals
    int n = sampleIndex;
    *left  = channelPointers[0][n];
    *right = channelPointers[1][n];
  }

  //void getFrameStereo(int sampleIndex, float* destination) const

protected:

  T* flatData = nullptr;         // pointer to the sample data
  T** channelPointers = nullptr;  // pointers to the channels
  // If we store the data in interleaved format, the channelPointers will be not needed and 
  // getFrame must be implemented differently. Maybe that's better (more efficient)
};

//=================================================================================================

/** A class to represent a pool of audio samples...tbc...

ToDo: currently, this is used in rsSamplerEngine - in a more general context, we would probably 
need a more complex referencing system for the AudioStream objects - regions would need to be 
something like AudioStreamClients, that register/deregister themselves, etc. Maybe AudioStream 
should be a subclass of some DataStream baseclass. We'll see... */

template<class T>
class SamplePool
{

public:


  ~SamplePool() { clear();  }


  int addSample(const AudioFileStream<T>* newSample)
  {
    // rsAssert(!contains(newSample))
    samples.push_back(newSample);
    return ((int) samples.size()) - 1;
  }
  // should the pool take ownership? ...i think so

  void removeSample(int i)
  {
    RAPT::rsAssert(i >= 0 && i < (int)samples.size());
    delete samples[i];
    RAPT::rsRemove(samples, (size_t)i);
  }


  /** Returns the number of samples in this pool. */
  int getNumSamples() const { return (int)samples.size(); }

  /** Returns true, if the given index i refers to a valid sample index. */
  bool isSampleIndexValid(int i) const { return i >= 0 && i < (int)samples.size(); }

  int findSample(const std::string& path) const;


  const AudioFileStream<T>* getSampleStream(int i) const
  {
    if(!isSampleIndexValid(i)) {
      RAPT::rsError("Invalid sample index");
      return nullptr; }
    return samples[i];
  }


  const std::string& getSamplePath(int i) const { return getSampleStream(i)->getPath();  }

  void clear();

  // todo:
  // setup: removeSample...but if regions refer to it, we need to update them, too by 
  // invalidating their pointers. We either need an observer mechanism (complex and general) or 
  // we allow removal of samples only via a member function of the outlying rsSamplerEngine 
  // class, which also takes care of resetting the sample-streams in all regions that use it
  // (simpler but less general). Maybe to make it safer, we could also introduce a reference
  // counter and check, if it is zero, before a stream objects gets removed
  // inquiry: hasSample

protected:

  std::vector<const AudioFileStream<T>*> samples;

};
// maybe rename to AudioFileStreamPool



}

#endif
