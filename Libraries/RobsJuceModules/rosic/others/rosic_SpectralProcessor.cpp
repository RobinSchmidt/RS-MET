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
  transformer.transformRealSignal(block, spectrum);

  //rsPlotComplexArrays(spectrumSize/2, (double*) spectrum);  // for debug
  //rsPlotComplexArrays(spectrumSize/2, (double*) spectrum, "Input spectrum");  // for debug

  processSpectrum(spectrum, spectrumSize);
  transformer.transformSymmetricSpectrum(spectrum, block);

  //rsPlotComplexArrays(spectrumSize/2, (double*) spectrum);  // for debug
  //rsPlotComplexArrays(spectrumSize/2, (double*) spectrum, "Output spectrum");  // for debug
  // ToDo: add titles - needs some API adaption for the plot function. We may need two separate 
  // functions for pltting one and two complex arrays


  // ToDo: 
  // -Get rid of the allocation for spectrum by keeping the buffer as member variable which should
  //  be allocated in the constructor. ...done - but we need to do this also in 
  //  SpectralEnvelopeProcessor and maybe in FormantShifter
}



