#include "ModalTests.h"

// move to PerformanceTestTools:
template<class T>
double getCyclesPerSample(T &module, int numSamples = 1000, int numTests = 5)
{  
  ::ProcessorCycleCounter counter;

  // create noise to be used as input signal:
  double *x = new double[numSamples];
  RAPT::rsArray::fillWithRandomValues(x, numSamples, -1.0, 1.0, 0);

  // do the test numTests times, use the minimum as result:
  double y, cycles, minCycles;
  minCycles = RS_INF(double);
  for(int i = 1; i <= numTests; i++)
  {
    counter.init();
    for(int n = 0; n < numSamples; n++)
      y = module.getSample(x[n]);
    cycles = (double) counter.getNumCyclesSinceInit();
    if( cycles < minCycles )
      minCycles = cycles;
  }

  // cleanup and return result:
  delete[] x;
  return minCycles / numSamples;
}

template<class T>
double getCyclesPerSampleBlockWise(T &module, 
                                   int numSamples = 10000, 
                                   int numTests = 5, 
                                   int blockSize = 512)
{  
  ::ProcessorCycleCounter counter;

  double *x = new double[numSamples];
  double *y = new double[numSamples];
  RAPT::rsArray::fillWithRandomValues(x, numSamples, -1.0, 1.0, 0);

  double cycles, minCycles;
  minCycles = RS_INF(double);
  int numFullBlocks = numSamples/blockSize;
  int lastBlockSize = numSamples - numFullBlocks*blockSize;

  for(int i = 1; i <= numTests; i++)
  {
    counter.init();

    int blockStart = 0;
    for(int n = 0; n < numFullBlocks; n++)
    {
      module.processBlock(&x[blockStart], &y[blockStart], blockSize);
      blockStart += blockSize;
    }
    module.processBlock(&x[blockStart], &y[blockStart], lastBlockSize);

    cycles = (double) counter.getNumCyclesSinceInit();
    if( cycles < minCycles )
      minCycles = cycles;
  }

  // cleanup and return result:
  delete[] x;
  return minCycles / numSamples;
}

void testModalFilter3(std::string &reportString)
{
  double fs  = 44100;  // samplerate in Hz
  double ta  = 0.02;   // decay time in seconds
  double td  = 0.1;    // decay time constant in seconds
  double f   = 55;     // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.5;    // amplitude as raw factor
  
  int blockSize  = 512;
  int numSamples = 10000;


  rsModalFilterDD mf;
  mf.setModalParameters(f, A, td, phs, fs);

  double cyclesPerSample = getCyclesPerSample(mf, 1000);
  printPerformanceTestResult("ModalFilter", cyclesPerSample);

  mf.reset();
  double cyclesPerSampleBlockWise = getCyclesPerSampleBlockWise(mf, numSamples, 5, blockSize);
  printPerformanceTestResult("  blockwise:", cyclesPerSampleBlockWise);

  rsModalFilterWithAttackDD mfa;
  mfa.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa, 1000);
  printPerformanceTestResult("ModalFilterAttack", cyclesPerSample);


  rsModalFilterWithAttackDD mfa2;
  mfa2.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa2, 1000);
  printPerformanceTestResult("ModalFilterAttack2", cyclesPerSample);

  /*
  rsNonlinearModalFilter nmf;
  nmf.setModalParameters(f, A, td, phs, fs);
  double cyclesPerSample = getCyclesPerSample(nmf);
  std::cout << cyclesPerSample;
  */
}

void testModalFilterBank(std::string &reportString)
{
  double sampleRate = 44100.0;  // samplerate in Hz
  double decay      = 0.5;      // decay time constant in seconds
  double frequency  = 55.0;     // fundamental frequency in Hz

  int numPartials = 100;
  int numSamples  = 10000;
  int blockSize   = 512;

  /*
  rsVectorDbl f = linearRangeVector(numPartials, 1.0, numPartials);
  rsVectorDbl a = applyFunction(f, -1.0,  &pow);
  rsVectorDbl d = applyFunction(f, -0.5,  &pow);
  rsVectorDbl p = randomVector(numPartials, 0.0, 360.0, 0); 
  */

  std::vector<double> f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  std::vector<double> a = rsApplyFunction(f, -0.7,  &pow); 
  std::vector<double> d = rsModalFilterBankDD::modeDecayTimes(f, 4.0, 0.95); 
  std::vector<double> p = rsRandomVector(numPartials, 0.0, 360.0, 0);

  rsModalFilterBankDD mfb;
  mfb.setSampleRate(sampleRate);
  mfb.setModalParameters(f, a, 0.1*d, d, p); // we need an operator taking a double and a vector
  mfb.setReferenceFrequency(frequency);
  mfb.setReferenceDecay(decay);

  mfb.reset();
  double cyclesPerSample = getCyclesPerSample(mfb, numSamples) / numPartials;
  printPerformanceTestResult("ModalFilterBank", cyclesPerSample);

  mfb.reset();
  double cyclesPerSampleBlockWise = 
    getCyclesPerSampleBlockWise(mfb, numSamples, 5, blockSize) / numPartials;
  printPerformanceTestResult("  blockwise:", cyclesPerSampleBlockWise);


  int dummy = 0;
}
