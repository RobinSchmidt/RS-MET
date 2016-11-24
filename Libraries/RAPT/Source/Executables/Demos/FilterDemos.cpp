#include "FilterDemos.h"
#include "../Shared/Plotting/Plotting.cpp"  // templates from there will be instantiated here

void ladderImpulseResponse()
{
  // Demonstrates, how to use the LadderFilter class. 

  typedef double Real; // real numbers are "float" - you can change this to "double"
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

void svfImpulseResponse()
{
  // Demonstrates, how to use the StateVariableFilter class. 

  typedef double Real;
  typedef RAPT::StateVariableFilter SVF;  // for convenience
  //typedef RAPT::StateVariableFilter<Real, Real> SVF;  // for convenience - to be used after templatization

  // create and set up the filter:
  RAPT::StateVariableFilter flt;
  flt.setSampleRate( 44100);
  flt.setFrequency(  1000);
  flt.setMode(SVF::LOWPASS);
  flt.setBandwidth(0.25);       // irrelevant for lowpass
  flt.setGain(2.0);

  // plot the impulse response:
  plotImpulseResponse(flt, 200, (Real) 1);
}