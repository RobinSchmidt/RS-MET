#include "FilterDemos.h"
#include "../Shared/Plotting/Plotting.cpp"  // templates from there will be instantiated here

void ladderImpulseResponse()
{
  // Plots an impulse response of the LadderFilter. 

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

void phasorImpulseResponse()
{
  // Plots an impulse response of the PhasorFilter. 

  typedef RAPT::PhasorFilter<float, float> PhsFlt;  // for convenience
  PhsFlt flt;
  flt.setSampleRate(44100); 
  flt.setFrequency(441);
  flt.setDecayTime(0.01f);
  //flt.setDecayTime(RS_INF(float));

  // plot the impulse response:
  plotImpulseResponse(flt, 1600, 1.f);
}

void svfImpulseResponse()
{
  // Plots an impulse response of the StateVariableFilter. 

  typedef float Real;
  typedef RAPT::StateVariableFilter<Real, Real> SVF;  // for convenience

  // create and set up the filter:
  SVF flt;
  flt.setSampleRate(44100); 
  flt.setFrequency(1000);
  flt.setMode(SVF::LOWPASS);
  flt.setBandwidth(0.25);       // irrelevant for lowpass
  flt.setGain(2.0);

  // plot the impulse response:
  plotImpulseResponse(flt, 200, (Real) 1);
}