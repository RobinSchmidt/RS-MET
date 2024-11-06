#ifndef RAPT_STATEVARIABLEFILTERSIMPER_H
#define RAPT_STATEVARIABLEFILTERSIMPER_H

/** An implementation of Andrew Simper's circuit modeled state variable filter from here:

https://www.cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf

It provides similar responses as the RBJ cookbook filters but the SVF is better suited to be used
as a VCF in the context of a synthesizer because it responds nicely to modulation.

This class implements the minimal core of the filter without any convenience features such as 
setters like setSampleRate, setFrequency, setMode, etc. We don't do this here because we don't want
to keep so many member variables around. These convenience features can be implemented on top of 
the core class in a subclass or in some object that embeds such an SVF core. Here, we have only the
minimum set of member variables that is needed to make the filter work. */


template<class T>                   // ToDo: have TSig, TPar template parameters
class rsStateVariableFilterSimper
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Enumeration of the available filter modes. */
  enum Mode
  {
    Bypass,
    Lowpass,
    Highpass,
    Bandpass,        // ToDo: rename to BandpassSkirt
    //BandpassPeak,
    Notch,
    Allpass,
    Bell,            // This is what RBJ calls "peak"
    LowShelf,
    HighShelf,
    Peak,            // Not the same as the "peak" characteristic in the RBJ filters


    NumModes
  };
  // ToDo: adjust the order of the modes to be the same as in the RBJ filter...but RBJ has two 
  // bandpass variants. I think, this here is a const skirt gain bandpass. Maybe to obtain const 
  // peak gain behavior, we just need to scale by k = 1/Q? ...just a guess - figure it out!
  // The RBJ filters are also missing a "peak" filter in the sense meant here. I think, it's just
  // a resonator? If so, try to introduce it in the RBJ filters as well. Maybe rename the mode to
  // "Reson" or "Resonator".
  // 
  // rosic::CookBookFilter has the modes in that order:  BYPASS = 0, LOWPASS, HIGHPASS, 
  // BANDPASS_CONST_SKIRT,  BANDPASS_CONST_PEAK, BANDREJECT, ALLPASS, PEAK, LOW_SHELF, HIGH_SHELF

  /** Sets up the filter coefficients so as to achieve the desired mode, cutoff, Q and gain. The 
  mode must be one of the values from the Mode enum, omega = 2*pi*freq/sampleRate is the usual 
  normalized radian frequency, Q is the quality factor which determines the resonance and A is the
  linear (!) gain for bell and shelf filter modes. */
  void setup(Mode mode, T omega, T Q, T A = T(1));


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. */
  inline T getSample(T in);

  /** Resets the internal state. */
  void reset() { ic1eq = ic2eq = 0; }


protected:

  // State:
  T ic1eq = 0;  // Maybe rename to i1
  T ic2eq = 0;
  // I think these may be currents into the two capacitors?

  // Coeffs:
  T a1 = 0, a2 = 0, a3 = 0;  // Filter coeffs (ToDo: explain better)
  T m0 = 1, m1 = 0, m2 = 0;  // Mixing coeffs

};

//-------------------------------------------------------------------------------------------------
// Implementation

template<class T>
void rsStateVariableFilterSimper<T>::setup(Mode mode, T omega, T Q, T A)
{
  // Prewarping cutoff, I guess:
  T tw2 = tan(0.5*omega); 

  // Helper function to calculate the a-coefficients from the intermediate variables g and k:
  auto calcFilterCoeffs = [&](T g, T k)
  {
    a1 = 1 / (1 + g*(g + k));
    a2 = g*a1;
    a3 = g*a2;
  };

  // Filter- and mixing coefficient calculations according to desired mode:
  switch(mode)
  {

  case Mode::Bypass:
  {
    m0 = 1;
    a1 = a2 = a3 = m1 = m2 = 0;
  }
  break;

  case Mode::Lowpass:
  {
    calcFilterCoeffs(tw2, 1/Q);
    m0 = 0;
    m1 = 0;
    m2 = 1;
  }
  break;

  case Mode::Highpass:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 = -1;
  }
  break;

  case Mode::Bandpass:
  {
    calcFilterCoeffs(tw2, 1/Q);
    m0 = 0;
    m1 = 1;
    m2 = 0;
  }
  break;

  case Mode::Notch:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 =  0;
  }
  break;

  case Mode::Peak:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 = -2;
  }
  break;

  case Mode::Allpass:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -2*k;
    m2 =  0;
  }
  break;

  case Mode::Bell:
  {
    T k = 1/(Q*A);
    calcFilterCoeffs(tw2, k);
    m0 = 1;
    m1 = k*(A*A - 1);
    m2 = 0;
  }
  break;

  case Mode::LowShelf:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2 / sqrt(A), k);
    m0 = 1;
    m1 = k*(A - 1);
    m2 = (A*A - 1);
  }
  break;

  case Mode::HighShelf:
  {
    T k = 1/Q;
    calcFilterCoeffs(tw2 * sqrt(A), k);
    m0 = A*A;
    m1 = k*(1 - A)*A;
    m2 = (1 - A*A);
  }
  break;

  default:
  {
    rsError("Unknown filter type in rsStateVarFilterSimper::setup");
    a1 = a2 = a3 = m0 = m1 = m2 = 0;  // We will produce a zero output in such a case.
  };

  }

  // ToDo:
  //
  // - Figure out and document what the coefficients and intermediat variables mean. Looking at the
  //   scribble on the front page on the paper, it seem like k is the feedback factor after the 1st
  //   integrator? And the a1, a2 are the gains of the integrators? And g is affecting them?
}

template<class T>
T rsStateVariableFilterSimper<T>::getSample(T v0)
{
  // Intermediate variables (voltages?):
  T v3 = v0 - ic2eq;                     // Feedback?
  T v1 = a1*ic1eq + a2*v3;
  T v2 = a2*ic1eq + a3*v3 + ic2eq;

  // State update (capacitor currents?):
  ic1eq = 2*v1 - ic1eq;
  ic2eq = 2*v2 - ic2eq;

  // Mix final output:
  return m0*v0 + m1*v1 + m2*v2;
}


// ToDo:
// - Document the code more - explain what the variables mean, etc.
// - Maybe have two template parameters TSig, TPar as in the other filters. I think,
//   v0,v1,v2,v3,ic1eq,ic2eq must all be TSig, a1,a2,a3,m1,m2,m3 must be TPar
// - Maybe move implementation into .cpp file ...but maybe not.
// - Figure out the z-domain transfer function and implement a function 
//   getTransferFunctionAt(rsComplex<TPar> z)


#endif