#include "FilterDemos.h"
#include "../../Modules/RAPT.cpp"
#include "../Shared/Plotting/Plotting.cpp"

void ladderImpulseResponse()
{
  // Demonstrates, how to use the LadderFilter class. 

  typedef float Real; // real numbers are "float" - you can change this to "double"
  typedef RAPT::LadderFilter<Real, Real> LadderFilterR;  // for convenience

  // create and set up the filter:
  LadderFilterR filter;
  filter.setSampleRate((Real)44100.0);
  filter.setCutoff(    (Real)1000.0);
  filter.setResonance( (Real)   0.8);
  filter.setMode(LadderFilterR::LP_24);

  // plot the impulse response (1000 samples long, scaled by 1):
  plotImpulseResponse(filter, 1000, (Real) 1);
}