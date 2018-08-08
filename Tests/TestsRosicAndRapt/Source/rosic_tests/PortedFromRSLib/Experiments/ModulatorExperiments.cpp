#include "ModulatorExperiments.h"

void breakpointModulator()
{
  double fs  = 100;  // samplerate in Hz

  rsBreakpointModulatorD bm;
  bm.setSampleRate(fs);

  // generate and plot envelope:
  static const int N = 500;  // number of samples to plot
  double t[N], x[N];
  createTimeAxis(N, t, fs);
  bm.noteOn();
  for(int n = 0; n < N; n++)
  {
    if( n == 300 )
      bm.noteOff();
    x[n] = bm.getSample();
  }
  plotData(N, t, x);
}

void breakpointModulatorSmoothFadeOut()
{
  double fs  = 100;  // samplerate in Hz

  rsBreakpointModulatorD bm;
  bm.setSampleRate(fs);
  bm.initialize();

  // generate and plot envelope:
  static const int N = 500;  // number of samples to plot
  bm.setBreakpointLevel(1, 0.0);
  bm.setBreakpointTime( 1, (N-1)/fs);
  bm.setBreakpointShape(1, rsModBreakpoint<double>::SMOOTH);

  double t[N], x[N];
  createTimeAxis(N, t, fs);
  bm.noteOn();
  for(int n = 0; n < N; n++)
    x[n] = bm.getSample();
  plotData(N, t, x);
}

void triSawModulator()
{
  double fs = 1000;          // sample rate
  static const int N = 5000;  // number of samples to plot

  rosic::rsTriSawModulator tsm;
  tsm.setSampleRate(fs);
  tsm.setAttackTime(20);
  tsm.setDecayTime(200);

  double t[N], x[N];
  createTimeAxis(N, t, fs);
  tsm.reset();
  for(int n = 0; n < N; n++)
    x[n] = tsm.getSample();
  plotData(N, t, x);
}

