//-------------------------------------------------------------------------------------------------
// Lifetime:

SpectralProcessor::SpectralProcessor(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: OverlapAddProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  transformer.setBlockSize(paddingFactor*blockSize);
  transformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);
  maxSpectrumSize = maxBlockSize * maxPaddingFactor / 2;
  spectrum = new Complex[maxSpectrumSize];
}

SpectralProcessor::~SpectralProcessor()
{
  delete[] spectrum;
}

//-------------------------------------------------------------------------------------------------
// Processing:

void SpectralProcessor::processBlock(double *block, int blockSize)
{
  int spectrumSize = blockSize/2;                    // /2 because we only get positive frequencies
  transformer.setBlockSize(blockSize);
  //Complex *spectrum = new Complex[spectrumSize];
  transformer.transformRealSignal(block, spectrum);
  processSpectrum(spectrum, spectrumSize);
  transformer.transformSymmetricSpectrum(spectrum, block);
  //delete[] spectrum;

  // ToDo: 
  // -Get rid of the allocation for spectrum by keeping the buffer as member variable which should
  //  be allocated in the constructor.
}



