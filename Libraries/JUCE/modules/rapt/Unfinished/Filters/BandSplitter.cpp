template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::setOmega(CRPar newOmega)
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

template<class TSig, class TPar>
std::complex<TPar> rsTwoBandSplitter<TSig, TPar>::getLowpassTransferFunctionAt(
  const std::complex<TPar>& z) const
{
  std::complex<TPar> zi = TPar(1)/z;
  return (b0 + b1*zi) / (TPar(1) + a1*zi); 
}

//-------------------------------------------------------------------------------------------------

template<class TSig, class TPar>
rsMultiBandSplitter<TSig, TPar>::~rsMultiBandSplitter()
{
  clearArrays();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSampleRate(CRPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateSplitters();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequency(int bandIndex, CRPar newFreq)
{
  //if(bandIndex <= getNumActiveBands()-1)
  //  return;  // a bit kludgy to avoid crash on MultiComp startup
  rsAssert(bandIndex >= 0 && bandIndex < splitFreqs.size());

  splitFreqs[bandIndex] = newFreq;
  splitters[bandIndex]->setOmega(TPar(2*PI)*newFreq/sampleRate);
}

// given n, this function returns a number r which is power-of-2-minus-1 such that r >= n
inline size_t newPowerOfTwoMinusOne(size_t n)
{
  size_t r = 1;
  while(r-1 < n)
    r *= 2;
  return r-1;
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequencies(const std::vector<TPar>& newFrequencies)
{
  //// old:
  //splitFreqs = newFrequencies;
  //updateSplitters();

  // new (needs testing):
  if(newFrequencies.size() > splitFreqs.size()) {
    size_t newSize = newPowerOfTwoMinusOne(newFrequencies.size());
    splitFreqs.resize(newSize);
  }
  for(size_t i = 0; i < newFrequencies.size(); i++)
    splitFreqs[i] = newFrequencies[i];
  numActiveBands = (int)newFrequencies.size()+1;
  updateSplitters();
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setNumberOfActiveBands(int newNumber)
{
  rsAssert(newNumber <= getNumAvailableBands());
  if(newNumber <= getNumAvailableBands())
    numActiveBands = newNumber;
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::addBand(CRPar splitFrequency)
{
  rsTwoBandSplitter<TSig, TPar>* splitter = new rsTwoBandSplitter<TSig, TPar>;
  splitter->setOmega(TPar(2*PI)*splitFrequency/sampleRate);
  splitters.push_back(splitter);        // later: insert sorted
  splitFreqs.push_back(splitFrequency); // here too

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::insertBand(int index, CRPar splitFrequency)
{
  rsTwoBandSplitter<TSig, TPar>* splitter = new rsTwoBandSplitter<TSig, TPar>;
  splitter->setOmega(TPar(2*PI)*splitFrequency/sampleRate);
  rsInsert(splitters,  splitter,       index);
  rsInsert(splitFreqs, splitFrequency, index);
  numActiveBands++;
}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::removeBand(int index, bool mergeWithRightNeighbour)
{
  numActiveBands--;
  rsRemove(splitters,  index);
  rsRemove(splitFreqs, index);


  // maybe we need to update a splitFreq?
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
  int oldNumBands = getNumAvailableBands();
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
std::complex<TPar> rsMultiBandSplitter<TSig, TPar>::getBandFrequencyResponseAt(
  int bandIndex, CRPar frequency) const
{
  TPar w = TPar(2*PI)*frequency/sampleRate;
  std::complex<TPar> j  = std::complex<TPar>(0, 1); // imaginary unit
  std::complex<TPar> z  = exp(j*w);                 // z-value where to evaluate H(z)
  std::complex<TPar> H = 1;                         // frequency response H(z=e^(j*w))
  int k;
  int numSplitters = numActiveBands - 1; // assume numActiveBands >= 1
  switch(mode)
  {
  case ACCUMULATE_INTO_HIGHPASS: 
  {
    // verify this (maybe it has to be the other way around)
    //for(k = 0; k < bandIndex; k++)
    //  H *= splitters[k]->getLowpassTransferFunctionAt(z);
    //for(k = bandIndex; k < numSplitters; k++)
    //  H *= splitters[k]->getHighpassTransferFunctionAt(z);

    // this certainly looks more plausible:
    for(k = 0; k < bandIndex; k++)
      H *= splitters[k]->getHighpassTransferFunctionAt(z);
    for(k = bandIndex; k < numSplitters; k++)
      H *= splitters[k]->getLowpassTransferFunctionAt(z);

  } break;

  // handle other cases...

  default: H = std::complex<TPar>(0, 0);
  }


  return H;
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
  int k;
  if(splitters.size() < splitFreqs.size())
    for(k = 0; k < splitters.size(); k++)
      delete splitters[k];

  splitters.resize(splitFreqs.size());

  for(k = 0; k < splitters.size(); k++) {
    splitters[k] = new rsTwoBandSplitter<TSig, TPar>;
    splitters[k]->setOmega(TPar(2*PI)*splitFreqs[k]/sampleRate);
  }


  /*
  // old:
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

  numActiveBands = getNumAvailableBands();
  */
}