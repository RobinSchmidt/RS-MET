#pragma once


/** Under Construction - just a stub at the moment.

A brickwall lowpass filter class that is made from a chain of 3 filter parts: (1) a lowpass, 
(2) a notch/bandstop, (3) an allpass. The lowpass is responsible for the general lowpass nature of
the filter. The notch is responsible for reducing the ringing at the cutoff frequency by notching 
some band around that frequency out. The allpass is responsible for moving a part of the ringing 
over to the left side of the edge, if you think in terms of the step response or a square wave 
input. The settings of these partial filters have been hand tuned to strike an optimal balance 
between the desirable steepness of the filter in the frequency domain and undesirable ringing of 
the filter in the time domain. The different modes that can be set by setMode switch between 
different configurations for the 3 filters that have been found to give good results.

...TBC...  

See the brickwallAndAllpass() in FilterExperiments.cpp and the BrickwallFilter presets for
ToolChain- this class is meant to encapsulate the findings of these experiments. */

template<class TSig, class TPar>
class rsBrickwallFilter
{

public:


  enum class Mode
  {
    halp12_halp4_ap8,   // 12th order Halpern lowpass, 4th order Halpern notch, 8th order allpass
    bess6               // 6th order Bessel lowpass, ....
  };


  void setMode(Mode newMode)
  {
    mode = newMode;
    // setDirty();
  }


protected:

  // User parameters:
  Mode mode       = Mode::halpern12;
  TPar sampleRate = TPar(44100);
  TPar cutoff     = TPar(1000);

  // Embedded objects:
  RAPT::rsEngineersFilter<TSig, TPar> lowpass;
  RAPT::rsEngineersFilter<TSig, TPar> notch;
  rosic::rsFlatZapper                 allpass;  // Maybe replace by rsAllpassChain

};
