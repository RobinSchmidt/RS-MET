//#include "rosic_SpectralProcessor.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SpectralProcessor::SpectralProcessor(int maxBlockSize, int maxOverlapFactor, int maxPaddingFactor)                                                   
: OverlapAddProcessor(maxBlockSize, maxOverlapFactor, maxPaddingFactor)
{
  transformer.setBlockSize(paddingFactor*blockSize);
  transformer.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO);
}

SpectralProcessor::~SpectralProcessor()
{

}

//-------------------------------------------------------------------------------------------------
// internal fucntions:

void SpectralProcessor::processBlock(double *block, int blockSize)
{
  int spectrumSize = blockSize/2;  // because we only get positive frequencies
  transformer.setBlockSize(blockSize);
  Complex *spectrum = new Complex[spectrumSize];
  transformer.transformRealSignal(block, spectrum);
  processSpectrum(spectrum, spectrumSize);
  transformer.transformSymmetricSpectrum(spectrum, block);
  delete[] spectrum;
}



