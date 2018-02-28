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
  rosic::Snowflake sf;

  int N = 50;  // number of samples
  int n;         // sample index

  std::vector<double> xt(N), yt(N), xf(N), yf(N); // table-based and on-the-fly generated outputs

  // create a unit square:
  //sf.clearRules();
  //sf.setAxiom("F+F+F+F");
  //sf.setNumIterations(0);
  //sf.setAngle(90);
  sf.setSampleRate(1.0);
  sf.setFrequency(0.35 / 4); // inc = frequency*numLines/sampleRate
  sf.updateAllInternals();

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
    rsAssert(xError == 0 && yError == 0);
  }

  


  GNUPlotter plt;
  plt.setRange(-1.1, +1.1, -1.1, +1.1);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");
  plt.addDataArrays(N, &xt[0], &yt[0]);
  plt.addDataArrays(N, &xf[0], &yf[0]);
  plt.plot();


  return false;
}