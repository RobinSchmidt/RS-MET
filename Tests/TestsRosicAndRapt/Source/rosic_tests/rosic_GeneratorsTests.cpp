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
  sf.setTurnAngle(90);
  sf.updateAllInternals();

  sf.setUseTable(false);

  sf.reset();
  sf.setResetRatio(0, 1);  // reset after each cycle
  for(n = 0; n < N; n++) 
    sf.getSampleFrameStereo(&x1[n], &y1[n]);

  sf.reset();
  sf.setResetRatio(0, 0);  // never reset
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

void fillResetInstantArray(rosic::ResetCounter& rc, int* r, int N)
{
  int i = 0;      // counts iterations/ticks
  int j = 0;      // counts resets
  rc.reset();
  while(true) {
    if(rc.tick()) {
      r[j] = i;
      j++;
      if(j == N)
        break;
    }
    i++;
  }
}

bool rotes::testResetter()
{
  bool result = true;

  static const int N = 15;
  int r[N];   // reset instants
  rosic::ResetCounter rc;

  rc.setInterval(5.0);
  fillResetInstantArray(rc, r, N);

  // check, if resets occur at expected times:
  int i;
  for(i = 0; i < N; i++)
  {

  }

  rc.setInterval(5.5);
  fillResetInstantArray(rc, r, N);
  // difference between entries should alternate between 5 and 6: 5,6,5,6,5,6,...

  rc.setInterval(5.25);  
  fillResetInstantArray(rc, r, N);
  // difference between entries tould follow the pattern: 5,5,5,6,5,5,5,6,...

  return result;
}


void rotes::testTurtleReverse()
{
  std::string commands = "F+F+F+F+F+F+F-F-F-F+F+F-F+F+F-F";

  TurtleGraphics tg;
  tg.setAngle(30);
  tg.init();

  std::vector<double> x, y;

  // forward part:
  x.push_back(tg.getX());
  y.push_back(tg.getY());
  int i = 0;
  while(i < commands.size()) {
    if(tg.interpretCharacter(commands[i])) {
      x.push_back(tg.getX());
      y.push_back(tg.getY()); }
    i++; }

  // reversed part:
  tg.setReverseMode(true);
  i = (int) commands.size()-1;
  while(i >= 0) {
    if(tg.interpretCharacter(commands[i])) {
      x.push_back(tg.getX());
      y.push_back(tg.getY()); }
    i--; }

  // plot:
  int N = (int) x.size();
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}


void setupTurtleSource(rosic::TurtleSource& ts)
{
  // used to set up a TurtleSource to specific settings (the idea is to call it for two objects to
  // make sure, they have the same settings)


  //ts.setTurtleCommands("F+F+F+F+F+");  // regular pentagon
  //ts.setTurnAngle(72);

  ts.setTurtleCommands("F+F+F+F+");   // square
  ts.setTurnAngle(90);

  //ts.setTurtleCommands("F+F+F+");   // triangle
  //ts.setTurnAngle(120);

  //ts.setTurtleCommands("F+F+");   // just forward/backward
  //ts.setTurnAngle(180);

  //ts.setTurtleCommands("FF");
  //ts.setTurnAngle(0);
  //// a simple "F" doesn't work right in forward mode, "FF" works right in forward but not reverse 
  //// mode

  ts.setSampleRate(1);
  ts.setFrequency(0.05);
  ts.setResetRatio(0, 0); // no reset
}

bool runTurtleTest(rosic::TurtleSource& ts1, rosic::TurtleSource& ts2, int numSamples, 
  bool plot = false)
{
  // compares the output of two TurtleSource objects and returns true, if they are almost equal

  int N = numSamples;
  std::vector<double> x1(N), y1(N), x2(N), y2(N);
  int n;
  double err;
  double tol = 1.e-14;
  bool result = true;

  for(n = 0; n < N; n++)
  {
    ts1.getSampleFrameStereo(&x1[n], &y1[n]);
    ts2.getSampleFrameStereo(&x2[n], &y2[n]);
    err = rmax(fabs(x2[n]-x1[n]), fabs(y2[n]-y1[n]));
    //rsAssert(ts2.checkIndexConsistency());

    result &= err <= tol;
    //rsAssert(result);
  }

  if(plot)
  {
    GNUPlotter plt;
    //plt.addDataArrays(N, &x1[0], &y1[0]);
    //plt.addDataArrays(N, &x2[0], &y2[0]);
    plt.addDataArrays(N, &x1[0]);
    plt.addDataArrays(N, &y1[0]);
    plt.addDataArrays(N, &x2[0]);
    plt.addDataArrays(N, &y2[0]);
    plt.plot();
  }

  return result;
}

void rotes::testTurtleSource()
{
  int N = 100; // number of samples
  bool result = true;

  // set up two TurtleSource objects with the same settings, the only difference being that one 
  // uses the table and the other one doesn't:
  rosic::TurtleSource ts1, ts2;
  ts1.setUseTable(true);

  setupTurtleSource(ts1);
  setupTurtleSource(ts2);
  ts1.setFrequencyScaler(-1.0);
  ts2.setFrequencyScaler(-1.0);
  result &= runTurtleTest(ts1, ts2, N, true);


  /*
  std::vector<double> x1(N), y1(N), x2(N), y2(N);
  int n;
  double err;
  double tol = 1.e-14;

  // for a simpler test, starting already in reverse mode:
  ts1.setFrequencyScaler(-1.0);
  ts2.setFrequencyScaler(-1.0);  
  // it already get wrong at sample index 1 - the lineBuffers x,y are wrong in ts2

  // forward:
  for(n = 0; n < N/2; n++)
  {
    ts1.getSampleFrameStereo(&x1[n], &y1[n]);
    ts2.getSampleFrameStereo(&x2[n], &y2[n]);

    err = rmax(fabs(x2[n]-x1[n]), fabs(y2[n]-y1[n]));
    rsAssert(ts2.checkIndexConsistency());
    //rsAssert(err < tol);
  }

  // backward:
  ts1.setFrequencyScaler(-1.0);
  ts2.setFrequencyScaler(-1.0);  
  for(n = N/2; n < N; n++)
  {
    ts1.getSampleFrameStereo(&x1[n], &y1[n]);
    ts2.getSampleFrameStereo(&x2[n], &y2[n]);

    err = rmax(fabs(x2[n]-x1[n]), fabs(y2[n]-y1[n]));
    rsAssert(ts2.checkIndexConsistency());
    //rsAssert(err < tol);
  }

  // make a simpler experiment: let the turtles start in reverse mode (don't switch mode during
  // run - it's already wrong when doing this) ...perhaps it starts out with a wrong reverse-flag?

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x1[0], &y1[0]);
  plt.addDataArrays(N, &x2[0], &y2[0]);
  //plt.addDataArrays(N, &x1[0]);
  //plt.addDataArrays(N, &y1[0]);
  //plt.addDataArrays(N, &x2[0]);
  //plt.addDataArrays(N, &y2[0]);
  plt.plot();
  */
}