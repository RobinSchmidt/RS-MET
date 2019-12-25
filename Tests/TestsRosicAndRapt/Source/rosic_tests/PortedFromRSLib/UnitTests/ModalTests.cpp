//#include "ModalTests.h"

bool testModalFilter2()
{
  bool testResult = true;

  // user parameters:
  static const int N = 1000;   // # samples to plot
  double fs  = 44100;  // samplerate in Hz
  double ta  = 0.05;   // attack time in seconds
  double td  = 0.1;    // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.5;    // amplitude as raw factor
  //double pm  = 5.0;    // phase-modulation (for the 2nd implementation)

  // create the target signal:
  double xt[N];
  double w = 2*PI*f/fs;    // normalized radian frequency
  double a = 1/(td*fs);    // normalized decay-rate
  double p = phs*PI/180.0; // start-phase in radians
  int n;
  for(n = 0; n < N; n++)
    xt[n] = A * exp(-a*n) * sin(w*n + p);

  // the impulse-response of an rsModalFilter should give the same signal:
  rsModalFilterDD mf;
  mf.setModalParameters(f, A, td, phs, fs);

  // check inquiry functions:
  testResult &= rsIsCloseTo(mf.getDecayTime(fs), td, 1.e-13);

  // check impulse-response:
  double x1[N];
  getImpulseResponse(mf, x1, N);
  double err = RAPT::rsArrayTools::maxDeviation(xt, x1, N);
  testResult &= err < 1.e-11;


  // the same procedure for an object of class rsNonlinearModalFilter:
  rsNonlinearModalFilterDD nmf;
  nmf.setModalParameters(f, A, td, phs, fs);
  double x2[N];
  getImpulseResponse(nmf, x2, N);
  err = RAPT::rsArrayTools::maxDeviation(xt, x2, N);
  testResult &= err < 1.e-13;  // rsNonlinearModalFilter has less error than rsModalFilter


  // check class rsModalFilterWithAttack:
  rsModalFilterWithAttackDD mfa;
  mfa.setModalParameters(f, A, ta, td, phs, fs);

  // check inquiry functions:
  double tmp = mfa.getLength(rsAmp2dB(1.0/EULER), fs);
  testResult &= rsIsCloseTo(tmp, ta+td, 1.e-13);

  // check impulse-response:
  double td2, scaler;
  expDiffScalerAndTau2(td, ta, &td2, &scaler);
  double a2 = 1/(td2*fs);    // normalized decay-rate 2
  for(n = 0; n < N; n++)
    xt[n] = A * scaler * (exp(-a*n)-exp(-a2*n)) * sin(w*n + p);
  getImpulseResponse(mfa, x1, N);
  err = RAPT::rsArrayTools::maxDeviation(xt, x1, N);
  testResult &= err < 1.e-11;

  // check class rsModalFilterWithAttack2:
  rsModalFilterWithAttack2DD mfa2;
  mfa2.setModalParameters(f, A, ta, td, phs, fs);
  getImpulseResponse(mfa2, x1, N);
  err = RAPT::rsArrayTools::maxDeviation(xt, x1, N);
  testResult &= err < 1.e-7;
    // 4 orders of magnitude less precise than rsModalFilterWithAttack (with GCC)

  // check rsModalFilterFloatSSE2:
  rsModalFilterFloatSSE2 mf_f32x4;
  double scaledA = A*scaler;
  mf_f32x4.setParametersTwoEnvs(w, scaledA, p, td2*fs, td2*fs, 0.0, td*fs, td*fs, 0.0);
  getImpulseResponse(mf_f32x4, x1, N);
  //plotData(N, 0, 1/fs, xt, x1); 
  err = RAPT::rsArrayTools::maxDeviation(xt, x1, N);
  testResult &= err < 0.002;

  return testResult;
}


bool testModalSynth()
{
  return true; 
  // maybe we can get rid of this "test" - it doesn't really test anything meaningful

  bool testResult = true;

  rsModalFilterBankDD ms;
  ms.setSampleRate(44100);
  ms.setReferenceFrequency(1000);
  ms.setReferenceDecay(0.2);

  static const int N = 100;
  double y[N];

  // 4 modes:
  double f4[4] = {1.0, 2.0, 3.0, 4.0};
  double a4[4] = {1.0, 0.5, 0.3, 0.4};
  double d4[4] = {1.0, 0.9, 0.8, 0.7};

  rsAssert(false); // this test is not yet fully updated bcs of this rsVectorDbl class - commented 
  // code below should be re-activated
  //rsVectorDbl f(4, f4);
  //rsVectorDbl g(4, a4);
  //rsVectorDbl d(4, d4);
  //rsVectorDbl p = rsModalFilterBankDD::randomModePhases(g);
  //ms.setModalParameters(f, g, 0.25*d, d, p);

  int n;
  y[0] = ms.getSample(1.0);
  for(n = 1; n < N; n++)
    y[n] = ms.getSample(0.0);

  testResult &= fabs(y[0]) < 1.e-17;  //y[0] should be close to zero

  return testResult;
}


