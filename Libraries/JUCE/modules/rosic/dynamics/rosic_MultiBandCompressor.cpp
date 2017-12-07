rsMultiBandCompressor::rsMultiBandCompressor() 
{
  tmpL.resize(maxNumBands);
  tmpR.resize(maxNumBands);
  compressors.resize(maxNumBands);
  for(int k = 0; k < maxNumBands; k++)
    compressors[k] = new rosic::Compressor;

  // init splitter frequencies for 16 bands (15 frequencies):
  vector<double> splitFreqs = { 60, 90, 135, 200, 300, 450, 675, 1000, 1500, 2200, 3000, 4500, 
    6000, 9000, 13000 }; 
  splitterL.setSplitFrequencies(splitFreqs);
  splitterR.setSplitFrequencies(splitFreqs);
  setNumberOfBands(16);
}

rsMultiBandCompressor::~rsMultiBandCompressor()
{
  for(int k = 0; k < maxNumBands; k++)
    delete compressors[k];
}

void rsMultiBandCompressor::setSampleRate(double newSampleRate)
{
  splitterL.setSampleRate(newSampleRate);
  splitterR.setSampleRate(newSampleRate);
  for(int k = 0; k < maxNumBands; k++)
    compressors[k]->setSampleRate(newSampleRate);
}

void rsMultiBandCompressor::setNumberOfBands(int newNumber)
{
  numBands = newNumber;
  splitterL.setNumberOfActiveBands(newNumber);
  splitterR.setNumberOfActiveBands(newNumber);
}

void rsMultiBandCompressor::setSplitFrequency(int bandIndex, double newFrequency)
{
  splitterL.setSplitFrequency(bandIndex, newFrequency);
  splitterR.setSplitFrequency(bandIndex, newFrequency);
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

void rsMultiBandCompressor::reset()
{
  splitterL.reset();
  splitterR.reset();
  for(int k = 0; k < maxNumBands; k++)
    compressors[k]->reset();
}
