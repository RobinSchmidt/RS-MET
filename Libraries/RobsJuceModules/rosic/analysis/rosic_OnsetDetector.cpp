#include "../_third_party/bernsee_fft/BernseeFFT.cpp"

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OnsetDetector::OnsetDetector()
{
  blockSize       = 1024;
  hopSize         = 128;    // blockSize/8 is a reasonable value for the hopSize
  sampleRate      = 44100;
  length          = 0;
  signal          = NULL;
  numBlocks       = 0;
  bufferIndex     = 0;
  nextBlockEnd    = blockSize-1;
  magnitudes      = new float[maxBlockSize/2];
  magnitudesOld   = new float[maxBlockSize/2];
  weights         = new float[maxBlockSize/2];
  window          = new float[maxBlockSize];
  linearBuffer    = new float[maxBlockSize];
  complexSpectrum = new float[2*maxBlockSize];
  circularBuffer  = new float[2*maxBlockSize];
  clearBuffers();
  createWindow();
}

OnsetDetector::~OnsetDetector()
{
  delete[] magnitudes;
  delete[] magnitudesOld;
  delete[] weights;
  delete[] window;
  delete[] linearBuffer;
  delete[] complexSpectrum;
  delete[] circularBuffer;
}

//-------------------------------------------------------------------------------------------------
// setup:

void OnsetDetector::setSampleRate(int newSampleRate)
{
  sampleRate = newSampleRate;

  if( sampleRate <= 48000 )
    blockSize = 1024;
  else if( sampleRate <= 96000 )
    blockSize = 2048;
  else
    blockSize = 4096;

  hopSize = blockSize/8;
  createWindow();
  computeSpectralWeights();
}

void OnsetDetector::reset()
{


}

//-------------------------------------------------------------------------------------------------
// processing:

void OnsetDetector::processSignalAtOnce(float *sampleData, int numSamples, int sampleRate)
{
  setSampleRate(sampleRate);  // adjusts blockSize and creates the window
  signal = sampleData;
  length = numSamples;
  onsets.clear();
  if( signal != NULL )
  {
    computeSpectralFlux();
    findOnsetsFromFluxMaxima();
  }
}

void OnsetDetector::prepareForBlockProcessing(int expectedSignalLength, int sampleRate)
{
  setSampleRate(sampleRate);  // adjusts blockSize and creates the window
  numBlocks    = (expectedSignalLength-blockSize) / hopSize + 1;
  bufferIndex  = 0;
  nextBlockEnd = blockSize-1;
  clearBuffers();
  rms.clear();
  flux.clear();
  onsets.clear();
  flux.reserve(numBlocks);
  rms.reserve(numBlocks);
}

void OnsetDetector::feedSignalBlock(float *sampleData, int numSamples)
{
  for(int n=0; n<numSamples; n++)
  {
    circularBuffer[bufferIndex] = sampleData[n];

    if( bufferIndex == nextBlockEnd )
    {
      // find the start-index for copying a chunk from the circular into the linear buffer:
      int start = bufferIndex - (blockSize-1);
      if( start < 0 )
        start += 2*maxBlockSize;

      // find the end index (keep track, if wraparound is necesarry and how much):
      int end         = start + (blockSize-1);
      int wrapSamples = 0;
      if( end >= 2*maxBlockSize )
      {
        wrapSamples = end - (2*maxBlockSize-1);
        end         = 2*maxBlockSize-1;
      }

      // copy the non-wrapped part of the data into the linear buffer:
      int i = start;
      int j = 0;
      while( i <= end )
      {
        linearBuffer[j] = circularBuffer[i];
        i++;
        j++;
      }

      // possibly, copy the wrapped part into the linear buffer:
      if( wrapSamples > 0 )
      {
        i = 0;
        while( i < wrapSamples )
        {
          linearBuffer[j] = circularBuffer[i];
          i++;
          j++;
        }
      }

      // process the current linear buffer:
      createComplexBlockForTransform(linearBuffer, complexSpectrum);
      smbFft(complexSpectrum, blockSize, -1);
      computeMagnitudes(complexSpectrum, magnitudes, blockSize);
      rms.push_back(  computeBlockRms(linearBuffer) );
      flux.push_back( computeSpectralFluxValue(magnitudes, magnitudesOld, weights) );
      RAPT::rsArrayTools::copy(magnitudes, magnitudesOld, blockSize/2);

      // set the mark for the next completed block:
      nextBlockEnd += hopSize;
      if( nextBlockEnd >= 2*maxBlockSize )
        nextBlockEnd -= 2*maxBlockSize;
    }

    bufferIndex++;
    if( bufferIndex >= 2*maxBlockSize )
      bufferIndex = 0;
  }
}

void OnsetDetector::finishBlockProcessing()
{
  numBlocks = (int)flux.size();
  if( numBlocks > 0 )
    flux[0] = 0.0;   // ignore very first value
  findOnsetsFromFluxMaxima();
}

//-------------------------------------------------------------------------------------------------
// inqiury:

/*
float OnsetDetector::getProgressInPercent()
{
  return 0.f;   // preliminary
}
*/

//-------------------------------------------------------------------------------------------------
// internal functions:

float* OnsetDetector::createAllZeroFloatArray(int numValues)
{
  float *a = new float[numValues];
  for(int n=0; n<numValues; n++)
    a[n] = 0.f;
  return a;
}

void OnsetDetector::clearBuffers()
{
  int n;
  for(n=0; n<maxBlockSize/2; n++)
  {
    magnitudes[n]    = 0.f;
    magnitudesOld[n] = 0.f;
  }
  for(n=0; n<maxBlockSize; n++)
  {
    linearBuffer[n] = 0.f;
  }
  for(n=0; n<2*maxBlockSize; n++)
  {
    circularBuffer[n]  = 0.f;
    complexSpectrum[n] = 0.f;
  }
}

void OnsetDetector::createWindow()
{
  // create a Hann-window:
  for(int n=0; n<blockSize; n++)
    window[n] = (float) (0.5 * ( 1.0 - cos(2.0*PI*n / (double) (blockSize-1)) ));
}

void OnsetDetector::computeSpectralWeights()
{
  // algorithm parameters (todo - make them members and let them be accessed from client code):
  float weight       = -1.f/4.f;     // exponent for spectral weighting values < 0 weight low
                                     // frequencies stronger
  float weightCutoff = 22000;        // cutoff frequency above which the flux does not contribute
                                     // to the onset's strength

  int n;
  for(n=0; n<blockSize/2; n++)
    weights[n] = pow((float) (n+1), weight);
  int kCutoff = (int) round(weightCutoff * (double)blockSize / (double)sampleRate);
  for(n=kCutoff; n<blockSize/2; n++)
    weights[n] = 0.f;
}

void OnsetDetector::createComplexBlockForTransform(float *realSignal, float *destination)
{
  for(int n=0; n<blockSize; n++)
  {
    destination[2*n]   = realSignal[n] * window[n];
    destination[2*n+1] = 0.f;
  }
}

void OnsetDetector::computeMagnitudes(float *complexSpectrum, float *magnitudes, int numBins)
{
  float scaler = (float) (1.0 / blockSize);
  for(int k=0; k<numBins/2; k++)
  {
    magnitudes[k] = scaler * sqrt(   complexSpectrum[2*k]   * complexSpectrum[2*k]
                                   + complexSpectrum[2*k+1] * complexSpectrum[2*k+1] );
  }
}

float OnsetDetector::computeBlockRms(float *block)
{
  float accu = 0.f;
  for(int n=0; n<blockSize; n++)
    accu += block[n]*block[n];
  accu /= (float) blockSize;
  return sqrt(accu);
}

float OnsetDetector::computeSpectralFluxValue(float *magnitudes, float *oldMagnitudes,
                                              float *spectralWeights)
{
  float accu = 0.0;
  float diff = 0.0;
  for(int k=0; k<blockSize/2; k++)
  {
    diff = magnitudes[k] - oldMagnitudes[k];
    if( diff > 0.f )
      accu += diff * spectralWeights[k];
  }
  return accu;
}

void OnsetDetector::computeSpectralFlux()
{
  // initializations:
  int blockIndex = 0;
  int blockStart = 0;
  int blockEnd   = blockStart + blockSize - 1;
  numBlocks      = (length-blockSize) / hopSize + 1;
  flux.reserve(numBlocks);
  rms.reserve(numBlocks);

  // compute first block seperately:
  createComplexBlockForTransform(&signal[blockStart], complexSpectrum);
  smbFft(complexSpectrum, blockSize, -1); // Stephan Bernsee's FFT routine
  computeMagnitudes(complexSpectrum, magnitudes, blockSize);
  rms.push_back(  computeBlockRms(&signal[blockStart]) );
  flux.push_back( 0.f );
  blockIndex += 1;
  blockStart += hopSize;
  blockEnd    = blockStart + blockSize - 1;
  RAPT::rsArrayTools::copy(magnitudes, magnitudesOld, blockSize/2);

  // loop over the blocks:
  while( blockEnd < length )
  {
    // compute magnitude spectrum fo current block:
    createComplexBlockForTransform(&signal[blockStart], complexSpectrum);
    smbFft(complexSpectrum, blockSize, -1);
    computeMagnitudes(complexSpectrum, magnitudes, blockSize);

    // compute the block's RMS and it's (weighted) spectral flux with respect to the previous
    // block:
    rms.push_back(  computeBlockRms(&signal[blockStart]) );
    flux.push_back( computeSpectralFluxValue(magnitudes, magnitudesOld, weights) );

    // increments and state-updates for next iteration:
    blockIndex += 1;
    blockStart += hopSize;
    blockEnd    = blockStart + blockSize - 1;
    RAPT::rsArrayTools::copy(magnitudes, magnitudesOld, blockSize/2);
  }
}

void OnsetDetector::findOnsetsFromFluxMaxima()
{
  // algorithm parameters (todo - make them members and let them be accessed from client code):
  int w              = blockSize/64;   // window size for local maximum
  int m              = 4;              // multiplier for extending the window into the past, 1...4
                                       // smaller -> more onsets
  int k              = (m+1)*w+1;      // size of the neighborhood for the local mean

  float threshold    = 0.7f;                  // relative threshold
  float absThreshold = (float)RAPT::rsDbToAmp(-60.0);  // absolute threshold in dB fo the RMS value of a
                                              // block to be considered at all
  bool  refineOnsets = false;                  // toggles onset time/strength refinement on

  for(int blockIndex = m*w+w; blockIndex<numBlocks; blockIndex++)
  {
    int   n         = blockIndex-w;                              // block under investigation
    float value     = flux[n];                                   // flux of the block
    float localMean = RAPT::rsArrayTools::mean(    &flux[n-m*w], k);  // mean of the block's neighborhood
    float localMax  = RAPT::rsArrayTools::maxValue(&flux[n-m*w], k);  // maximum in the block's neighborhood

    // onset decision function - an onset must the local maximum of the chunk and some threshold
    // above the local mean (modified with respect to Dixon so as to use a relative threshold):
    bool onsetDetected =     value  == localMax
                         &&  rms[n] >= absThreshold
                         &&  value  >= localMean + threshold*localMean;

    if( onsetDetected == true )
    {
      Onset newOnset;
      newOnset.timeInSamples = (n+1)*hopSize + blockSize/2 - 1;
      newOnset.strength      = value;

      // if so desired, refine the time and strength of the onset by fitting a quadratic parabola
      // through the local maximum and its two neighbours and locate the maximum of the parabola:
      if( refineOnsets == true )
      {
        float t[3];
        t[1] = (float) newOnset.timeInSamples;
        t[0] = t[1] - hopSize;
        t[2] = t[1] + hopSize;
        float s[3];
        s[1] = value;
        s[0] = flux[n-1];
        s[2] = flux[n+1];
        float a, b, c;
        fitQuadratic(t, s, a, b, c);
        float tExact = -b / (2.f*a);
        float sExact = a*tExact*tExact + b*tExact + c;
        newOnset.timeInSamples = (int) round(tExact);
        newOnset.strength      = sExact;
      }

      onsets.push_back(newOnset);
    }

  }
}

void OnsetDetector::fitQuadratic(float *x, float *y, float &a, float &b, float &c)
{
  float k1 = y[1]      - y[0];
  float k2 = x[0]*x[0] - x[1]*x[1];
  float k3 = x[1]      - x[0];
  float k4 = k1/k3;
  float k5 = k2/k3;
  a = (y[2]-k4*x[2]-y[0]+k4*x[0]) / (x[2]*x[2]+k5*x[2]-x[0]*x[0]-k5*x[0]);
  b = (k1+k2*a)/k3;
  c = y[0]-a*x[0]*x[0]-b*x[0];
}
