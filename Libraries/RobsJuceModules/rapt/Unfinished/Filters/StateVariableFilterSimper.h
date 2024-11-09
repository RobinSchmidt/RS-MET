#ifndef RAPT_STATEVARIABLEFILTERSIMPER_H
#define RAPT_STATEVARIABLEFILTERSIMPER_H

/** An implementation of Andrew Simper's circuit modeled state variable filter from here:

      https://www.cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf

It provides the same responses as the RBJ cookbook filters but the SVF is better suited to be used
as a VCF in the context of a synthesizer because it responds nicely to modulation.

This class implements the minimal core of the filter without any convenience features such as 
setters like setSampleRate, setFrequency, setMode, etc. We don't do this here because we don't want
to keep so many member variables around. These convenience features can be implemented on top of 
the core class in a subclass or in some object that embeds such an SVF core. Here, we have only the
minimum set of member variables that is needed to make the filter work. */


template<class TSig, class TPar>
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
    BandpassSkirt,   // Bandpass with constant skirt gain
    BandpassPeak,    // Bandpass with constant peak gain
    Notch,
    Allpass,
    Bell,            // This is what RBJ calls "peak"
    LowShelf,
    HighShelf,
    Peak,            // Not the same as the "peak" characteristic in the RBJ filters
                     // Maybe rename to resonator (but first figure out if it actually is one)

    NumModes
  };

  /** Sets up the filter coefficients so as to achieve the desired mode, cutoff, Q and gain. The 
  mode must be one of the values from the Mode enum, omega = 2*pi*freq/sampleRate is the usual 
  normalized radian frequency, Q is the quality factor which determines the resonance and A is the
  linear (!) gain for bell and shelf filter modes. */
  void setup(Mode mode, TPar omega, TPar Q, TPar A = TPar(1));




  /** UNDER CONSTRUCTION....Does not yet work!
  Sets up the filter coefficients to simulate a biquad filter with given coeffs. */
  void setupFromBiquad(TPar b0, TPar b1, TPar b2, TPar a1, TPar a2);

  /** Initializes all coefficients to achieve a neutral "bypass" response. */
  void initCoeffs()
  {
    a1 = 0, a2 = 0, a3 = 0;
    m0 = 1, m1 = 0, m2 = 0;
  }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. */
  inline TSig getSample(TSig in);

  /** Resets the internal state. */
  void reset() { i1 = i2 = 0; }


protected:

  /** Helper function to calculate and assign the a-coefficients from the intermediate variables 
  g and k.  ToDo: explain meaning of g and k.  */
  void calcFilterCoeffs(TPar g, TPar k);


  // State:
  TSig i1 = 0, i2 = 0;          // Capacitor equivalent(?) currents ic1eq, ic2eq in the paper.

  // Coeffs:
  TPar a1 = 0, a2 = 0, a3 = 0;  // Filter coeffs (ToDo: explain better)
  TPar m0 = 1, m1 = 0, m2 = 0;  // Mixing coeffs

};

//-------------------------------------------------------------------------------------------------
// Implementation

template<class TSig, class TPar>
void rsStateVariableFilterSimper<TSig, TPar>:: calcFilterCoeffs(TPar g, TPar k)
{
  a1 = 1 / (1 + g*(g + k));
  a2 = g*a1;
  a3 = g*a2;
}

template<class TSig, class TPar>
void rsStateVariableFilterSimper<TSig, TPar>::setup(Mode mode, TPar omega, TPar Q, TPar A)
{
  // Prewarping cutoff (I guess):
  TPar tw2 = tan(0.5*omega);

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
    TPar k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 = -1;
  }
  break;

  case Mode::BandpassSkirt:
  {
    calcFilterCoeffs(tw2, 1/Q);
    m0 = 0;
    m1 = 1;
    m2 = 0;
  }
  break;

  case Mode::BandpassPeak:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 = 0;
    m1 = k;
    m2 = 0;
  }
  break;

  case Mode::Notch:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 =  0;
  }
  break;

  case Mode::Allpass:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -2*k;
    m2 =  0;
  }
  break;

  case Mode::Bell:
  {
    TPar k = 1/(Q*A);
    calcFilterCoeffs(tw2, k);
    m0 = 1;
    m1 = k*(A*A - 1);
    m2 = 0;
  }
  break;

  case Mode::LowShelf:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2 / sqrt(A), k);
    m0 = 1;
    m1 = k*(A - 1);
    m2 = (A*A - 1);
  }
  break;

  case Mode::HighShelf:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2 * sqrt(A), k);
    m0 = A*A;
    m1 = k*(1 - A)*A;
    m2 = (1 - A*A);
  }
  break;

  case Mode::Peak:
  {
    TPar k = 1/Q;
    calcFilterCoeffs(tw2, k);
    m0 =  1;
    m1 = -k;
    m2 = -2;
  }
  break;

  default:
  {
    rsError("Unknown filter type in rsStateVarFilterSimper::setup");
    a1 = a2 = a3 = m0 = m1 = m2 = 0;  
    // We will produce an output of zero in such a case.
  };

  }

}

template<class TSig, class TPar>
TSig rsStateVariableFilterSimper<TSig, TPar>::getSample(TSig v0)
{
  // Compute node voltages:
  TSig v3 = v0 - i2;                       // Feedback (?)
  TSig v1 = a1*i1 + a2*v3;                 // Voltage at node 1, Bandpass output (?)
  TSig v2 = a2*i1 + a3*v3 + i2;            // Voltage at node 2, Lowpass output (?)

  // State update (by computing "equivalent"(?) capacitor currents):
  i1 = TPar(2)*v1 - i1;                    // Eq. 2?
  i2 = TPar(2)*v2 - i2;

  // Mix final output:
  return m0*v0 + m1*v1 + m2*v2;
}


#endif