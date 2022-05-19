template<class T>
bool AudioFileStreamPreloaded<T>::setData(
  T** newData, int numFrames, int numDataChannels, double sampleRate, int numStreamChannels, 
  const std::string& path)
{
  RAPT::rsAssert(numFrames > 0 && numDataChannels > 0, "Data must be non-empty");


  // Deallocate old and allocate new memory:
  clear();
  if(numFrames <= 0 || numDataChannels <= 0 || numStreamChannels <= 0) {
    RAPT::rsError("Data should be non-empty"); 
    return false; }
    // I'm not sure what would be the most meaningful way to handle a situation where the caller
    // passes some invalid data but it seems sensible to just reset the object into a cleared 
    // state. If this happens, it probably means we have a bug at a higher level. That's why the
    // assertion is there.

  int numChannelsMin = RAPT::rsMin(numDataChannels, numStreamChannels);
  flatData = new T[numChannelsMin*numFrames];
  channelPointers = new T*[numStreamChannels];
  if(flatData == nullptr || channelPointers == nullptr) {
    clear(); return false; }  // memory allocation failed

  // Copy the new data into the freshly allocated memory:
  for(int c = 0; c < numChannelsMin; c++) 
    for(int n = 0; n < numFrames; n++)
      flatData[c*numFrames + n] = newData[c][n];
  for(int c = 0; c < numStreamChannels; c++) {
    int channelStart = (c % numChannelsMin) * numFrames;
    channelPointers[c] = &flatData[channelStart]; }

  // Update metadata members and report success:
  this->numChannels = numStreamChannels;
  this->numFrames   = numFrames;
  this->sampleRate  = sampleRate;
  this->path        = path;
  return true; // success

  // Maybe we should have a version of this function which does not need to copy data but instead
  // just takes over ownership of the passed array. But this would need a parameter for the flat 
  // data array, too. But maybe the caller has no flat array. We'll see, how this meshes with the
  // wavefile loading functions. For the time being, we copy...
}

/*
void AudioFileStreamPreloaded::setNumOutputChannels(int newNumChannels)
{
  int dummy = 0;
}
*/

template<class T>
void AudioFileStreamPreloaded<T>::clear()
{
  AudioFileStream<T>::clear();
  delete[] channelPointers;
  channelPointers = nullptr; 
  delete[] flatData;
  flatData = nullptr;
}

//-------------------------------------------------------------------------------------------------

template<class T>
int SamplePool<T>::findSample(const std::string& path) const
{
  for(size_t i = 0; i < samples.size(); i++) {
    if(samples[i]->getPath() == path)
      return (int) i; }
  return -1;
  // ToDo: 
  // -Maybe keep the samplePool sorted, so we can use binary search here. That requires to
  //  reorder it when new samples are added, but addition of new samples is costly anyway due to 
  //  disk access, so that probably doesn't really matter.
}

template<class T>
void SamplePool<T>::clear()
{
  for(size_t i = 0; i < samples.size(); i++)
    delete samples[i];
  samples.clear();
}

template<class T>
void SamplePool<T>::copyContent(const SamplePool<T>& other)
{
  clear();
  for(int i = 0; i < other.getNumSamples(); i++)
  {

    int dummy = 0;
  }



  RAPT::rsError("Not yet implemented");
}
