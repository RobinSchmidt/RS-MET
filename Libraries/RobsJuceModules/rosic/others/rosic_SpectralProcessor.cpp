//-------------------------------------------------------------------------------------------------
// Lifetime:

SpectralProcessor::SpectralProcessor(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: OverlapAddProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  transformer.setBlockSize(maxPaddingFactor*maxBlockSize); // Allocate enough memory for worst case
  transformer.setBlockSize(paddingFactor*blockSize);       // Set up actual blocks size
  transformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);
  spectrum = new Complex[getMaxSpectrumSize()];
}

SpectralProcessor::~SpectralProcessor()
{
  delete[] spectrum;
}

//-------------------------------------------------------------------------------------------------
// Processing:

void SpectralProcessor::processBlock(double *block, int blockSize)
{
  int spectrumSize = blockSize/2;
  transformer.setBlockSize(blockSize);
  //Complex *spectrum = new Complex[spectrumSize];
  transformer.transformRealSignal(block, spectrum);
  processSpectrum(spectrum, spectrumSize);
  transformer.transformSymmetricSpectrum(spectrum, block);
  //delete[] spectrum;

  // ToDo: 
  // -Get rid of the allocation for spectrum by keeping the buffer as member variable which should
  //  be allocated in the constructor. ...done - but we need to do this also in 
  //  SpectralEnvelopeProcessor and maybe in FormantShifter
}



