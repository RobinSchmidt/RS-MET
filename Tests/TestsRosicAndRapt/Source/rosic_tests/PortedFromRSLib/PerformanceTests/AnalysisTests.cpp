#include "AnalysisTests.h"

typedef std::complex<double> rsComplexDbl;

void testFourierTransformer(std::string &reportString)
{
  static const int bufferSize = 4096;
  static const int numBuffers = 32;

  rsComplexDbl x[bufferSize];
  rsComplexDbl X[bufferSize];
  RAPT::rsArrayTools::fillWithRandomValues(x, bufferSize, -1.0, +1.0, 1);

  rsFourierTransformerRadix2D ft;
  ft.setBlockSize(bufferSize);
  ::PerformanceCounterTSC counter;
  counter.init();
  int b;
  for(b = 0; b < numBuffers; b++)
    ft.transformComplexBuffer(x, X);
  double cycles = (double) counter.getNumCyclesSinceInit();
  double cyclesPerSample = cycles / (bufferSize*numBuffers);
  printPerformanceTestResult("FourierTransformerRadix2", cyclesPerSample);

  counter.init();
  for(b = 0; b < numBuffers; b++)
  {
    RAPT::rsArrayTools::copy(x, X, bufferSize);
    rsFFT(X, bufferSize);
  }
  cycles = (double) counter.getNumCyclesSinceInit();
  cyclesPerSample = cycles / (bufferSize*numBuffers);
  printPerformanceTestResult("rsFFT", cyclesPerSample);

}

void testAutoCorrelationPitchDetector2(std::string &reportString)
{
  static const int bufferSize     = 2048;
  static const int updateInterval = 256;
  static const int blockSize      = 128;
  static const int numBlocks      = 16;

  double fs = 44100;   // samplerate in Hz
  double f  = 1000.0;  // frequency of the sinusoid

  // create and setup the pitch detector object:
  rsAutoCorrelationPitchDetectorD pd;
  pd.setBufferSize(bufferSize);
  pd.setUpdateInterval(updateInterval);
  pd.setSampleRate(fs);
  pd.setMaxFundamental(8000.0);

  // do the performance test:
  double w = 2*PI*f/fs;
  double fe; // frequency estimate (irrelevant in performance test)
  double block[blockSize];   
  int b, n;
  ::PerformanceCounterTSC counter;
  counter.init();
  for(b = 0; b < numBlocks; b++)
  {
    for(n = 0; n < blockSize; n++)
      block[n] = sin( w * (b*blockSize+n) );
    fe = pd.processBlock(block, blockSize);
  }
  double cycles = (double) counter.getNumCyclesSinceInit();
  double cyclesPerBlock  = cycles/numBlocks;
  double cyclesPerSample = cyclesPerBlock/blockSize;
  printPerformanceTestResult("AutoCorrelationPitchDetector", cyclesPerSample);
}
