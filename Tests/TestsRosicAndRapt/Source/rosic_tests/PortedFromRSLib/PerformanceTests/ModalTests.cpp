#include "ModalTests.h"

// move to PerformanceTestTools:
template<class TMod, class TSig>
double getCyclesPerSample(TMod &module, int numSamples = 1000, int numTests = 5, 
  TSig dummy = 1.0) // the dummy is to let the compiler determine the signal type
{  
  ::PerformanceCounterTSC counter;
  //::PerformanceCounterQPC counter;
  //::PerformanceCounterPMC counter;  // requires privileged instructions

  // create noise to be used as input signal:
  TSig *x = new TSig[numSamples];
  RAPT::rsArray::fillWithRandomValues(x, numSamples, -1.0, 1.0, 0);

  // do the test numTests times, use the minimum as result:
  TSig y; 
  double cycles, minCycles;
  minCycles = RS_INF(double);
  for(int i = 1; i <= numTests; i++)
  {
    counter.init();
    for(int n = 0; n < numSamples; n++)
      y = module.getSample(x[n]);
    cycles = (double) counter.getNumCyclesSinceInit();
    if( cycles < minCycles && cycles > 0 )
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
  ::PerformanceCounterTSC counter;

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
  double ta  = 0.1;    // decay time in seconds
  double td  = 0.5;    // decay time constant in seconds
  double f   = 55;     // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.5;    // amplitude as raw factor

  int blockSize  = 512;
  int numSamples = 2048;
  //int numSamples = 1000000;
  int numTests   = 15;


  double w = 2*PI*f/fs;
  double p = RAPT::rsDegreeToRadiant(phs);

  rsModalFilterDD mf;
  mf.setModalParameters(f, A, td, phs, fs);

  double cyclesPerSample = getCyclesPerSample(mf, numSamples, numTests, 1.0);
  printPerformanceTestResult("ModalFilter", cyclesPerSample);

  mf.reset();
  double cyclesPerSampleBlockWise = getCyclesPerSampleBlockWise(mf, numSamples, numTests, blockSize);
  printPerformanceTestResult("  blockwise:", cyclesPerSampleBlockWise);

  rsModalFilterWithAttackDD mfa;
  mfa.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa, numSamples, numTests, 1.0);
  printPerformanceTestResult("ModalFilterAttack", cyclesPerSample);

  rsModalFilterWithAttackDD mfa2;
  mfa2.setModalParameters(f, A, ta, td, phs, fs);
  cyclesPerSample = getCyclesPerSample(mfa2, numSamples, numTests, 1.0);
  printPerformanceTestResult("ModalFilterAttack2", cyclesPerSample);

  rsModalFilterFloatSSE2 mf4; // 4 because of the 4 sinusoids
  mf4.setParameters(w, A, p, 0.1*ta*fs, ta*fs, 0.5, 0.1*td*fs, td*fs, 0.5);
  cyclesPerSample = getCyclesPerSample(mf4, numSamples, numTests, rsFloat32x4(1.f));
  //cyclesPerSample = getCyclesPerSample(mf4, numSamples, numTests, 1.f);
  printPerformanceTestResult("ModalFilterFloatSSE2", cyclesPerSample);




  // todo: check, how just instatiating the template with vector types performs:
  //rsModalFilterWithAttack<rsFloat64x2, rsFloat64x2> mfa64x2;

  /*
  rsNonlinearModalFilter nmf;
  nmf.setModalParameters(f, A, td, phs, fs);
  double cyclesPerSample = getCyclesPerSample(nmf);
  std::cout << cyclesPerSample;
  */

  // The rsModalFilterFloatSSE2 version performs very badly when the difference equation is 
  // implemented correctly - compare that to the stereo version of EngineersFilter which takes
  // 13.75 cycles per sample and biquad - for scalar double and rsFloat64x2. Try a modal filter
  // based on rsFloat64x2 ...and/or maybe the inlining in the loop over the stages helps? Try
  // a ModalBank with rsFloat32x4. ...when removing certain simple state-update instructions from
  // the difference equation, it runs like lightning - 7 cycles per sample and mode - but adding
  // the update back let's the cycles skyrocket to several hundreds :-(  ...this is very weird!
  // maybe other implementation structures could be better (try a rotating phasor). Maybe also
  // make a test of rsBiquadCascade<rsFloat32x4, rsFloat32x4> - see if it also takes just 13.75
  // cycles per 4-vector - try it also on the other machine
  // maybe check also the output signals - could it be that the filter is unstable and the 
  // arithmetic produces exceptions/errors that depend on which coeffs are active? denormals
  // or something? yes - actually, the coeffs are even denormal - figure out, why - fix it and
  // test again - yeah - fuck - it was the denormals!
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
  double cyclesPerSample = getCyclesPerSample(mfb, numSamples, 5, 1.0) / numPartials;
  printPerformanceTestResult("ModalFilterBank", cyclesPerSample);

  mfb.reset();
  double cyclesPerSampleBlockWise = 
    getCyclesPerSampleBlockWise(mfb, numSamples, 5, blockSize) / numPartials;
  printPerformanceTestResult("  blockwise:", cyclesPerSampleBlockWise);


  int dummy = 0;
}
