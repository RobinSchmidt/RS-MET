#include "ModalTests.h"

template<class T>
double getCyclesPerSample(T &module, int numSamples = 1000, int numTests = 5)
{  
  ProcessorCycleCounter counter;

  // create noise to be used as input signal:
  double *x = new double[numSamples];
  RSLib::rsFillWithRandomValues(x, numSamples, -1.0, 1.0, 0);

  // do the test numTests times, use the minimum as result:
  double y, cycles, minCycles;
  minCycles = rsInfDouble;
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
  ProcessorCycleCounter counter;

  double *x = new double[numSamples];
  double *y = new double[numSamples];
  RSLib::rsFillWithRandomValues(x, numSamples, -1.0, 1.0, 0);

  double cycles, minCycles;
  minCycles = rsInfDouble;
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

void testModalFilter2(std::string &reportString)
{
  double fs  = 44100;  // samplerate in Hz
  double ta  = 0.02;   // decay time in seconds
  double td  = 0.1;    // decay time constant in seconds
  double f   = 55;     // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.5;    // amplitude as raw factor
  
  int blockSize  = 512;
  int numSamples = 10000;


  rsModalFilter mf;
  mf.setModalParameters(f, A, td, phs, fs);

  double cyclesPerSample = getCyclesPerSample(mf, 1000);
  appendResultToReport(reportString, "ModalFilter", cyclesPerSample);

  mf.reset();
  double cyclesPerSampleBlockWise = getCyclesPerSampleBlockWise(mf, numSamples, 5, blockSize);
  appendResultToReport(reportString, "  blockwise:", cyclesPerSampleBlockWise);

  rsModalFilterWithAttack mfa;
  mfa.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa, 1000);
  appendResultToReport(reportString, "ModalFilterAttack", cyclesPerSample);


  rsModalFilterWithAttack mfa2;
  mfa2.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa2, 1000);
  appendResultToReport(reportString, "ModalFilterAttack2", cyclesPerSample);

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

  rsVectorDbl f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  rsVectorDbl a = rsApplyFunction(f, -0.7,  &pow); 
  rsVectorDbl d = rsModalFilterBank::modeDecayTimes(f, 4.0, 0.95); 
  rsVectorDbl p = rsRandomVector(numPartials, 0.0, 360.0, 0);

  rsModalFilterBank mfb;
  mfb.setSampleRate(sampleRate);
  mfb.setModalParameters(f, a, 0.1*d, d, p);
  mfb.setReferenceFrequency(frequency);
  mfb.setReferenceDecay(decay);

  mfb.resetModalFilters();
  double cyclesPerSample = getCyclesPerSample(mfb, numSamples) / numPartials;
  appendResultToReport(reportString, "ModalFilterBank", cyclesPerSample);

  mfb.resetModalFilters();
  double cyclesPerSampleBlockWise = 
    getCyclesPerSampleBlockWise(mfb, numSamples, 5, blockSize) / numPartials;
  appendResultToReport(reportString, "  blockwise:", cyclesPerSampleBlockWise);


  int dummy = 0;
}
