template<class T>
bool AudioFileStreamPreloaded<T>::setData(
  T** newData, int numFrames, int numDataChannels, T sampleRate, int numStreamChannels, 
  const std::string& path)
{
  // Deallocate old and allocate new memory:
  clear();
  int numChannelsMin = rsMin(numDataChannels, numStreamChannels);
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
  // data array, too. We'll see, how this meshes with the wavefile loading functions...
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
  numChannels = 0;
  numFrames   = 0; 

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
}

template<class T>
void SamplePool<T>::clear()
{
  for(size_t i = 0; i < samples.size(); i++)
    delete samples[i];
  samples.clear();
}