#pragma once


//=================================================================================================

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
ToolChain - this class is meant to encapsulate the findings of these experiments. */

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


//=================================================================================================

template<class T>
struct rsStateVariableFilterCoeffs
{
  T g;          // Embedded integrator gain.
  T R2pg;       // 2*R + g
  T h;          // 1 / (1 + 2*R*g + g*g), factor for zero delay feedback prediction.
  T cL;         // Coefficient for lowpass output.
  T cB;         // Coefficient for bandpass output.
  T cH;         // Coefficient for highpass output.

  // Notation:
  // -R is the damping, R2 == 2*R == 1/Q.
};

template<class T>
struct rsStateVariableFilterState
{
  T s1;         // State of first integrator
  T s2;         // State of second integrator
};

template<class TSig, class TCoef>
struct rsStateVariableFilterData
{
  rsStateVariableFilterCoeffs<TCoef> coeffs;
  rsStateVariableFilterState<TSig>   state;
};


 /** A class that implements a chain of zero delay feedback (ZDF) state variable filters (SVFs). It
 can be used to replace a rsBiquadCascade ...TBC...  */

template<class TSig, class TCoef>  // types for signal and coefficients
class rsStateVariableFilterChain
{


public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setupFrom(const RAPT::rsBiquadCascade<TSig, TCoef>& biquadChain);


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig in);



protected:

  inline TSig getStageOutput(int stage, TSig in);


  std::vector<rsStateVariableFilterData<TSig, TCoef>> data;

};



template<class TSig, class TCoef> 
TSig rsStateVariableFilterChain<TSig, TCoef>::getStageOutput(int stage, TSig in)
{
  // Shorthands for convenience:
  rsStateVariableFilterCoeffs<TCoef>& c = data[stage].coeffs;
  rsStateVariableFilterState<TSig>&   s = data[stage].state;

  // Compute the 3 outputs (LP, BP, HP):
  TSig yH = (in - c.R2pg * s.s1 - s.s2) * c.h;
  TSig yB = c.g*yH + s.s1;
  s.s1    = c.g*yH + yB; 
  TSig yL = c.g*yB + s.s2;
  s.s2    = c.g*yB + yL;

  // Combine the 3 outputs to final output:
  return c.cL*yL + c.cB*yB + c.cH*yH;

  // See comments in rsStateVariableFilter<TSig, TPar>::getSample for what's going on
}

template<class TSig, class TCoef> 
TSig rsStateVariableFilterChain<TSig, TCoef>::getSample(TSig in)
{
  TSig y = in;
  for(int i = 0; i < (int)data.size(); i++)  // ToDo: try to avoid conversion
    y = getStageOutput(i, y);
  return y;
}

