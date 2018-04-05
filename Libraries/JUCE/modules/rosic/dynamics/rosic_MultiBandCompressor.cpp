rsMultiBandEffect::rsMultiBandEffect()
{
  initBands();

  /*
  tmpL.resize(maxNumBands);
  tmpR.resize(maxNumBands);
  //indices.resize(maxNumBands);
  //initIndices();
  */


  /*
  // init splitter frequencies for 16 bands (15 frequencies):
  vector<double> splitFreqs = { 60, 90, 135, 200, 300, 450, 675, 1000, 1500, 2200, 3000, 4500, 
    6000, 9000, 13000 }; 
  splitterL.setSplitFrequencies(splitFreqs);
  splitterR.setSplitFrequencies(splitFreqs);
  setNumberOfBands(16);
  */
}

// setup:

void rsMultiBandEffect::setSampleRate(double newSampleRate)
{
  splitter.setSampleRate(newSampleRate);
  splitterL.setSampleRate(newSampleRate);
  splitterR.setSampleRate(newSampleRate);
}

/*
void rsMultiBandEffect::setNumberOfBands(int newNumber)
{
  numBands = newNumber;
  splitterL.setNumberOfActiveBands(newNumber);
  splitterR.setNumberOfActiveBands(newNumber);
}
*/

void rsMultiBandEffect::insertBand(int i, double splitFreq)
{
  numBands++;
  splitter.insertBand(i, splitFreq);
  splitterL.insertBand(i, splitFreq);
  splitterR.insertBand(i, splitFreq);
  tmp.resize(numBands);
  tmpL.resize(numBands);
  tmpR.resize(numBands);
}

void rsMultiBandEffect::removeBand(int i, bool mergeRight)
{
  numBands--;
  splitter.removeBand(i, mergeRight);
  splitterL.removeBand(i, mergeRight);
  splitterR.removeBand(i, mergeRight);
  tmp.resize(numBands);
  tmpL.resize(numBands);
  tmpR.resize(numBands);
}

void rsMultiBandEffect::setSplitMode(int newMode)
{
  splitter.setSplitMode(newMode);
  splitterL.setSplitMode(newMode);
  splitterR.setSplitMode(newMode);
}

void rsMultiBandEffect::setSplitFrequency(int bandIndex, double newFrequency)
{
  splitter.setSplitFrequency(bandIndex, newFrequency);
  splitterL.setSplitFrequency(bandIndex, newFrequency);
  splitterR.setSplitFrequency(bandIndex, newFrequency);
}

void rsMultiBandEffect::initBands()
{
  numBands = 1;

  splitter.setNumberOfActiveBands(numBands);
  splitterL.setNumberOfActiveBands(numBands);
  splitterR.setNumberOfActiveBands(numBands);

  tmp.resize(numBands);
  tmpL.resize(numBands);
  tmpR.resize(numBands);
}

double rsMultiBandEffect::getDecibelsAt(int index, double frequency)
{
  return 0; // not yet implemented
}

// processing:

void rsMultiBandEffect::reset()
{
  splitter.reset();
  splitterL.reset();
  splitterR.reset();
}

/*
void rsMultiBandEffect::initIndices()
{
  fillWithIndex(&indices[0], (int)indices.size());
}
*/
//=================================================================================================

rsMultiBandCompressor::rsMultiBandCompressor() 
{
  initBands();
}

rsMultiBandCompressor::~rsMultiBandCompressor()
{
  clearCompressors();
}

void rsMultiBandCompressor::setSampleRate(double newSampleRate)
{
  rsMultiBandEffect::setSampleRate(newSampleRate);
  for(size_t k = 0; k < compressors.size(); k++)
    compressors[k]->setSampleRate(newSampleRate);
}

void rsMultiBandCompressor::setThreshold(int bandIndex, double newThreshold)
{
  compressors[bandIndex]->setThreshold(newThreshold);
}

void rsMultiBandCompressor::setRatio(int bandIndex, double newRatio)
{
  compressors[bandIndex]->setRatio(newRatio);
}

void rsMultiBandCompressor::setAttackTime(int bandIndex, double newAttackTime)
{
  compressors[bandIndex]->setAttackTime(newAttackTime);
}

void rsMultiBandCompressor::setReleaseTime(int bandIndex, double newReleaseTime)
{
  compressors[bandIndex]->setReleaseTime(newReleaseTime);
}

void rsMultiBandCompressor::insertBand(int i, double splitFreq)
{
  rsMultiBandEffect::insertBand(i, splitFreq);
  RAPT::rsInsertValue(compressors, new rosic::Compressor, i);
  // todo: set up compressor setting to the same values as compressor at i (or i-1)
}

void rsMultiBandCompressor::removeBand(int i, bool mergeRight)
{
  rsMultiBandEffect::removeBand(i, mergeRight);
  delete compressors[i];
  RAPT::rsRemove(compressors, i);
  // maybe adjust splitFreq of splitter i or i-1?
}

void rsMultiBandCompressor::initBands()
{
  rsMultiBandEffect::initBands();
  clearCompressors();
  compressors.resize(getNumberOfBands());      // should be initially 1
  for(size_t k = 0; k < compressors.size(); k++)
    compressors[k] = new rosic::Compressor;
}

void rsMultiBandCompressor::reset()
{
  rsMultiBandEffect::reset();
  for(size_t k = 0; k < compressors.size(); k++)
    compressors[k]->reset();
}

void rsMultiBandCompressor::clearCompressors()
{
  for(size_t k = 0; k < compressors.size(); k++)
    delete compressors[k];
  compressors.resize(0);
}