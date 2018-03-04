#include "rosic_GeneratorsTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
#include "../Shared/Plotting/rosic_Plotter.h"
#include "../Shared/Plotting/GNUPlotter.h"
using namespace rosic;

void rotes::testOscillatorStereo()
{
  // set up the WaveTable:
  static const int prototypeLength = 2048;
  double prototypeWaveFlat[prototypeLength];
  double *prototypeWave[2];
  prototypeWave[0] = &prototypeWaveFlat[0];
  prototypeWave[1] = &prototypeWaveFlat[0];
  int n;
  for(n=0; n<prototypeLength; n++)
    prototypeWaveFlat[n] = sawWave( 2*PI*n / prototypeLength);
  MipMappedWaveTableStereo waveTable;
  waveTable.setWaveform(prototypeWave, prototypeLength);

  // set up the oscillator:
  OscillatorStereo osc;
  osc.setWaveTableToUse(&waveTable);
  osc.setStereoPhaseShift(90.0);
  osc.setStartPhase(95.0);

  // obtain the data for the plot:
  static const int plotLength = 160;
  double plotIndices[plotLength];
  fillWithIndex(plotIndices, plotLength);

  double plotDataFlat1[2*plotLength];
  double *plotData1[2];
  plotData1[0] = &plotDataFlat1[0];
  plotData1[1] = &plotDataFlat1[plotLength];
  osc.getWaveformForDisplay(plotData1, plotLength);

  osc.setStartPhase(96.0);
  double plotDataFlat2[2*plotLength];
  double *plotData2[2];
  plotData2[0] = &plotDataFlat2[0];
  plotData2[1] = &plotDataFlat2[plotLength];
  osc.getWaveformForDisplay(plotData2, plotLength);


  Plotter::plotData(plotLength, plotIndices, plotData1[0], plotData2[0]);
}

void rotes::testLorentzSystem()
{
  static const int N = 2000;


  rosic::LorentzSystem lorentzSystem;
  lorentzSystem.setPseudoFrequency(1000);


  double t[N], x[N], y[N], z[N];
  lorentzSystem.getState(&x[0], &y[0], &z[0]);
  for(int n = 1; n < N; n++)
  {
    lorentzSystem.iterateState();
    lorentzSystem.getState(&x[n], &y[n], &z[n]);
  }

  rosic::fillWithIndex(t, N);
  Plotter::plotData(N, t, x, y, z);
}

bool rotes::testSnowflake()
{
  bool testResult = true;

  rosic::Snowflake sf;
  sf.setSampleRate(1.0);
  sf.setFrequency(0.35 / 4); // inc = frequency*numLines/sampleRate

  int N = 1000;  // number of samples
  int n;         // sample index
  std::vector<double> xt(N), yt(N), xf(N), yf(N); // table-based and on-the-fly generated outputs

  // create on-the-fly output:
  sf.setUseTable(false);
  sf.reset();
  for(n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xf[n], &yf[n]);

  // create table-based output:
  sf.setUseTable(true);
  sf.reset();
  for(n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xt[n], &yt[n]);

  // xt should equal xf, same for yt,yf - check this and use it as return value
  double xError = 0, yError = 0;
  for(n = 0; n < N; n++)
  {
    xError = rmax(xError, fabs(xt[n]-xf[n]));
    yError = rmax(yError, fabs(yt[n]-yf[n]));
    testResult &= xError == 0 && yError == 0;
    rsAssert(testResult);
  }

  // compare outputs with and without reset in on-the-fly mode (there are some audible buzzy 
  // artifacts when reset is off - it should be equal up to numerical roundoff - but isn't):

  std::vector<double> x0(N), y0(N), x1(N), y1(N); // reset modes 0 (never) and 1 (each cycle)

  // create a unit square:
  sf.clearRules();
  sf.setAxiom("F+F+F+F+");
  sf.setNumIterations(0);
  sf.setAngle(90);
  sf.updateAllInternals();

  sf.setUseTable(false);

  sf.reset();
  sf.setResetAfterCycles(1);  // reset after each cycle
  for(n = 0; n < N; n++) 
    sf.getSampleFrameStereo(&x1[n], &y1[n]);

  sf.reset();
  sf.setResetAfterCycles(0);  // never reset
  for(n = 0; n < N; n++) 
    sf.getSampleFrameStereo(&x0[n], &y0[n]);


  GNUPlotter plt;
  plt.setRange(-1.1, +1.1, -1.1, +1.1);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");

  //plt.addDataArrays(N, &xt[0], &yt[0]);
  //plt.addDataArrays(N, &xf[0], &yf[0]);

  plt.addDataArrays(N, &x0[0], &y0[0]);
  plt.addDataArrays(N, &x1[0], &y1[0]);
  plt.plot();


  return testResult;
}

bool rotes::testResetter()
{
  bool result = true;



  return result;
}

void rotes::testSnowflakeResetting()
{
  // We test how the apparent "modulation" frequency depends on the reset interval, number of
  // turtle lines and signal frequency

  // create a snowflake that produces a Moore curve:
  rosic::Snowflake sf;
  sf.clearRules();
  sf.setAxiom("LFL+F+LFL+F+");
  sf.addRule('L', "-RF+LFL+FR-");
  sf.addRule('R', "+LF-RFR-FL+");
  sf.setAngle(90);
  sf.setUseTable(false);

  // these are the parameters, on which this modulation frequency depends - tweak them:
  sf.setResetAfterCycles(1);
  sf.setResetAfterLines(66);
  sf.setNumIterations(2);    // 0->4, 1->16, 2->64, 3->256, 4->1024
  sf.setSampleRate(8192.0);
  sf.setFrequency(24.0);
  sf.updateAllInternals();


  // create test output
  int N = 8192;  // number of samples  
  int n;         // sample index
  std::vector<double> x(N), y(N), r(N);
  for(n = 0; n < N; n++)
  {
    sf.getSampleFrameStereo(&x[n], &y[n]);
    if(sf.getLineCount() == 0)
      r[n] = 1;
  }

  std::vector<int> lineCountResets;
  for(n = 1; n < N; n++) {
    if(r[n] == 1 && r[n-1] == 0)
      lineCountResets.push_back(n); }



  //rosic::writeToStereoWaveFile("MooreCurveResetting.wav", &x[0], &y[0], N, 8192, 16);

  // plot left and right signal aginst sample index:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays(N, &y[0]);
  plt.addDataArrays(N, &r[0]); // samples, weher lineCount == 0
  plt.plot();

  // Observations:
  // Notation: L: num lines, C: cycle count reset interval, R: line count reset interval, 
  // f: frequency, fs: sample rate, P: reset period
  // Resets due to lineCount 
  // fs=8192, L=64, C=1, f=32: R=63: P=k*256, R=64: P=k*256, R=65: P=k*260, R=66: P=k*264, k int
  // fs=8192, L=64, C=1, f=24: R=63: P=342,683,1024,1366,1707,2048,2390,2731,3072,3414
  //                           R=64: dito
  //                           R=65: P=347,694,1040,1387,1734,2080,2427,2774,3120,3467
  //                           R=66:

  /*
  -that seems to depend on the difference between (a multiple of) numLines and lineCountReset 
    if lineCountReset is close to a multiple of numLines, modulation is slow - but the multiple
      can also be 1.5 ...but also, the higher the multiple, the slower and the higher the played
      note, the faster, so 
      speed = k * (numLines-resetInterval) * noteFreq?
      that k depends in some way on the ratio resetInterval/numLines - maybe gcd/lcm is involved?
   */

  // or: speed = modFreq = fm = noteFreq * F(resetInterval R, numLines L, cycleResetInterval C), 
  // assume C = const
  // fm = f * F(R,L) - we want a formula for F(R,L), the formula has to involve max(R,L) and R-L
  // fm ?= a * f * abs(L-max(R,L))

  // ...naah - i think we need a resetIntervalInSeconds = sampleRate / (resetFreqFactor * noteFreqinHz)
  // and compute the resetInterval as numLines*resetIntervalInSeconds
  // ...or something

  // give the user a parameter p, when p=0 then R=inf, when p=1 then R=L -> R = L/p ...but that
  // doesn't involve the note-frequency ...actually we want also, when inc increases, that the 
  // difference between R and (some multiple of) L decreases, let D = R-L (or L-R), we want
  // D = const/f

  // no - we should have R = p*L and p should have keytracking like:
  // p(key) = p0 + k*(key-referenceKey)
}