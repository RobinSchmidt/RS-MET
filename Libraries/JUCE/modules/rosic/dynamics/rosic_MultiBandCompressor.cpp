rsMultiBandCompressor::rsMultiBandCompressor() 
{
  tmpL.resize(maxNumBands);
  tmpR.resize(maxNumBands);
  compressors.resize(maxNumBands);
  for(int k = 0; k < maxNumBands; k++)
    compressors[k] = new rosic::Compressor;
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

}

void rsMultiBandCompressor::reset()
{
  splitterL.reset();
  splitterR.reset();
  for(int k = 0; k < maxNumBands; k++)
    compressors[k]->reset();
}

