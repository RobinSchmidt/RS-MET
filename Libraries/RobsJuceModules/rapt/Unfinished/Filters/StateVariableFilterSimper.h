#ifndef RAPT_STATEVARIABLEFILTERSIMPER_H
#define RAPT_STATEVARIABLEFILTERSIMPER_H

/** An implementation of Andrew Simper's circuit modeled state variable filter from here:

https://www.cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf  

*/


template<class T>             // have TSig, TPar template parameters
class rsStateVarFilterSimper
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
    Bandpass,
    Notch,
    Peak,            // Not the same as the "peak" characteristic in the RBJ filters
    Allpass,
    Bell,            // This is what RBJ calls "peak"
    LowShelf,
    HighShelf,

    NumModes
  };

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
  void reset()
  {
    ic1eq = 0;
    ic2eq = 0;
  }


protected:

  // State:
  T ic1eq = 0;
  T ic2eq = 0;
  // I think these may be currents into the two capacitors?

  // Coeffs:
  T a1 = 0, a2 = 0, a3 = 0;  // Filter coeffs (ToDo: explain better)
  T m0 = 1, m1 = 0, m2 = 0;  // Mixing coeffs

};

//-------------------------------------------------------------------------------------------------
// Implementation

template<class T>
void rsStateVarFilterSimper<T>::setup(Mode mode, T omega, T Q, T A)
{
  // Prewarping cutoff, I guess:
  T tw2 = tan(0.5*omega); 

  // Helper function to calculate the a-coefficients from g and k:
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
T rsStateVarFilterSimper<T>::getSample(T v0)
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


#endif