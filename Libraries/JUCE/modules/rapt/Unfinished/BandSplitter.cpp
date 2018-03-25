template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::setOmega(TPar newOmega)
{
  w = newOmega;
  TPar c = tan(TPar(0.5) * w);
  c = (c-1) / (c+1);
  b0 = b1 = TPar(0.5) * (1+c); // use just one b variable b01 or something
  a1 = c;  // use a1 directly in computation instead of temporary c
  // Formulas from DAFX, page 40.
}

template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::copySettingsFrom(const rsTwoBandSplitter& s)
{
  w  = s.w; 
  b0 = s.b0; 
  b1 = s.b1;
  a1 = s.a1;
}

template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::copyStateFrom(const rsTwoBandSplitter& s)
{
  x1 = s.x1;
  y1 = s.y1; 
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsMultiBandSplitter<TSig, TPar>::~rsMultiBandSplitter()
{
  clearArrays();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateSplitters();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequency(int bandIndex, TPar newFreq)
{
  splitFreqs[bandIndex] = newFreq;
  splitters[bandIndex]->setOmega(TPar(2*PI)*newFreq/sampleRate);
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequencies(const std::vector<TPar>& newFrequencies)
{
  splitFreqs = newFrequencies;
  updateSplitters();

  //// it's not always necessarry to clear and create everything anew - optimize this later
  //clearArrays();
  //for(int i = 0; i < newFrequencies.size(); i++)
  //  addBand(newFrequencies[i]);
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setNumberOfActiveBands(int newNumber)
{
  rsAssert(newNumber <= getNumBands());
  numActiveBands = newNumber;
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::addBand(TPar splitFrequency)
{
  rsTwoBandSplitter<TSig, TPar>* splitter = new rsTwoBandSplitter<TSig, TPar>;
  splitter->setOmega(TPar(2*PI)*splitFrequency/sampleRate);
  splitters.push_back(splitter);        // later: insert sorted
  splitFreqs.push_back(splitFrequency); // here too

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::insertBand(int index, TPar splitFrequency)
{
  rsTwoBandSplitter<TSig, TPar>* splitter = new rsTwoBandSplitter<TSig, TPar>;
  splitter->setOmega(TPar(2*PI)*splitFrequency/sampleRate);
  rsInsert(splitters,  splitter,       index);
  rsInsert(splitFreqs, splitFrequency, index);
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::removeBand(int index, bool mergeWithRightNeighbour)
{
  rsRemove(splitters,  index);
  rsRemove(splitFreqs, index);
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::copyBandSettings(int src, int dst, bool copyState)
{
  splitters[dst]->copySettingsFrom(*splitters[src]);
  if(copyState)
    splitters[dst]->copyStateFrom(*splitters[src]);
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::reset()
{
  for(size_t i = 0; i < splitters.size(); i++)
    splitters[i]->reset();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setNumberOfBands(int newNumBands)
{
  int oldNumBands = getNumBands();
  if(newNumBands == oldNumBands)
    return; // nothing to do
  if(newNumBands > oldNumBands)
  {
    TPar loFreq = splitFreqs[oldNumBands-2];
    //TPar hiFreq = TPar(0.5) * sampleRate;
    splitFreqs.resize(newNumBands-1);
    for(int k = 0; k < (newNumBands-oldNumBands); k++)
      splitFreqs[oldNumBands-2+k] = loFreq; // preliminary ...
    // ..new bands should have frequencies spread between currently highest and nyquist freq
  }
  else
    splitFreqs.resize(newNumBands-1);
  updateSplitters();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::clearArrays()
{
  for(size_t i = 0; i < splitters.size(); i++)
    delete splitters[i];
  splitters.clear();
  splitFreqs.clear();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::updateSplitters()
{
  // delete superfluous splitters:
  for(size_t i = splitFreqs.size(); i < splitters.size(); i++)
    delete splitters[i];
  splitters.resize(splitFreqs.size());

  // create new splitters as necessary and set them up:
  for(size_t k = 0; k < splitters.size(); k++)
  {
    if(splitters[k] == nullptr)
      splitters[k] = new rsTwoBandSplitter<TSig, TPar>;
    splitters[k]->setOmega(TPar(2*PI)*splitFreqs[k]/sampleRate);
  }
}