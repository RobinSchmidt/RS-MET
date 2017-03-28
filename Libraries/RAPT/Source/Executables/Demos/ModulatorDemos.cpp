#include "ModulatorDemos.h"

void breakpointModulatorDefault()
{
  // Creates a BreakpointModulator object and plots its default output signal - i.e. without 
  // setting up any breakpoints. This will produce an attack-decay-sustain-release (ADSR) shape.

  int N = 8000;                    // total length in samples
  int noteOffAt = 5000;            // sample at which we trigger the note-off
  vector<double> y(N);

  rsBreakpointModulatorF envGen;
  //rsBreakpointModulator<double> envGen;
  envGen.setTimeScale(0.05);       // otherwise, it would be too long for plotting
  envGen.noteOn();
  for(int n = 0; n < N; n++)
  {
    if(n == noteOffAt)
      envGen.noteOff();
    y[n] = envGen.getSample();
  }

  GNUPlotter plt;
  plt.addDataArrays(N, &y[0]);
  plt.plot();

  // todo: later plot the output against an x-axis scaled in seconds (we need the Array-fill functions 
  // there...)
}
