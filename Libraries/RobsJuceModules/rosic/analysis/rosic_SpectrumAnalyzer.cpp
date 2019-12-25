//#include "rosic_SpectrumAnalyzer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SpectrumAnalyzer::SpectrumAnalyzer()
{
  sampleRate    = 44100.0;
  numChannels   = 2;
  blockSize     = 8192;
  fftSize       = 1*blockSize;
  windowType    = RAPT::rsWindowFunction::WindowType::hanningZZ; // maybe use a flat-top window
  sampleCounter = 0;
  midSideMode   = false;

  fourierTransformer.setBlockSize(blockSize);

  // init the pointer to pointer for the magnitudes:
  for(int c=0; c<maxNumChannels; c++)
  {
    magnitudePointer[c] = &(magnitudeBuffer[c][0]);
  }

  setSampleRate(sampleRate);
  clearBuffers();
  makeWindow();
}

SpectrumAnalyzer::~SpectrumAnalyzer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void SpectrumAnalyzer::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  int k;
  for(k=0; k<blockSize; k++)
    frequencyBuffer[k] = k*sampleRate / (double) blockSize;
  for(k=blockSize; k<maxBlockSize; k++)
    frequencyBuffer[k] = 0.0;

  // a hack to avoid log-of zero in logarithmic plots:
  frequencyBuffer[0] = 1.0;
}

void SpectrumAnalyzer::setBlockSize(int newBlockSize)
{
  // check new blocksize for validity:
  if( newBlockSize >= 2 && RAPT::rsIsPowerOfTwo(newBlockSize) )
  {
    blockSize = newBlockSize;
    fftSize   = blockSize;

    clearBuffers();
    makeWindow();

    int k;
    for(k=0; k<blockSize; k++)
      frequencyBuffer[k] = k*sampleRate / (double) blockSize;
    for(k=blockSize; k<maxBlockSize; k++)
      frequencyBuffer[k] = 0.0;

    // a hack to avoid log-of zero in logarithmic plots:
    frequencyBuffer[0] = 1.0;
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SpectrumAnalyzer::clearBuffers()
{
  for(int c=0; c<numChannels; c++)
  {
    for(int n=0; n<maxBlockSize; n++)
    {
      inBuffer[c][n]        = 0.0;
      tmpBuffer[c][n]       = 0.0;
      magnitudeBuffer[c][n] = 0.0;
    }
  }
}

void SpectrumAnalyzer::updateDisplayBuffers()
{
  // copy the content of the circular buffer into a temporary array:
  int c, k;
  int n_r, n_w; // read and write positions

  for(c=0; c<numChannels; c++)
  {
    for(n_w=0; n_w<blockSize; n_w++)
    {
      n_r = n_w + sampleCounter - blockSize;
      while( n_r < 0 )
        n_r += maxBlockSize;
      while( n_r >= maxBlockSize )
        n_r -= maxBlockSize;


      tmpBuffer[c][n_w] = windowBuffer[n_w] * inBuffer[c][n_r];
    }
    // zero padding:
    for(n_w=blockSize; n_w<fftSize; n_w++)
      tmpBuffer[c][n_w] = 0.0;
  }

  fourierTransformer.setBlockSize(fftSize);
  fourierTransformer.getRealSignalMagnitudes( &(tmpBuffer[0][0]), &(magnitudeBuffer[0][0]) );
  fourierTransformer.getRealSignalMagnitudes( &(tmpBuffer[1][0]), &(magnitudeBuffer[1][0]) );

  double normalizer = 2.0 / blockSize;
  for(c=0; c<numChannels; c++)
  {
    for(k=0; k<blockSize; k++)
    {
      magnitudeBuffer[c][k] *= normalizer;
    }
  }
}

void SpectrumAnalyzer::makeWindow()
{
  int n;
  for(n=0; n<maxBlockSize; n++)
    windowBuffer[n] = 0.0;

  RAPT::rsWindowFunction::createWindow(windowBuffer, blockSize, windowType, true);
  // todo: use flat top window - better for ampltude estimation

  /*
  // normalize the window to unit mean (maybe move this to RAPT::rsWindowFunction or RAPT::rsArrayTools):
  double sum          = 0.0;
  for(n=0; n<blockSize; n++)
    sum        += windowBuffer[n];
  double mean   = sum / (double) blockSize;
  double scaler = 1.0 / mean;
  for(n=0; n<blockSize; n++)
    windowBuffer[n] *= scaler;
    */
}