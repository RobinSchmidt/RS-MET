#ifndef RAPT_STATEVARIABLEFILTERMYSTRAN_H
#define RAPT_STATEVARIABLEFILTERMYSTRAN_H


template<class TSig, class TPar> // signal, parameter types
class rsStateVariableFilterMystran
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
    BandpassSkirt,
    BandpassPeak,
    Bandstop,
    Allpass,
    Bell,
    LowShelf,
    HighShelf,

    NumModes
  };

  void setup(Mode mode, TPar omega, TPar Q, TPar A = TPar(1));
  // Convenience function


  // Separate setup functions for the different modes to allow to bypass the switch-statement in 
  // the general setup function

  void setupBypass();
  void setupLowpass(      TPar omega, TPar Q);
  void setupHighpass(     TPar omega, TPar Q);
  void setupBandpassSkirt(TPar omega, TPar Q);
  void setupBandpassPeak( TPar omega, TPar Q);
  void setupBandstop(     TPar omega, TPar Q);
  void setupAllpass(      TPar omega, TPar Q);
  void setupBell(         TPar omega, TPar Q, TPar A);
  void setupLowShelf(     TPar omega, TPar Q, TPar A);
  void setupHighShelf(    TPar omega, TPar Q, TPar A);


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. */
  inline TSig getSample(TSig in);

  /** Resets the internal state. */
  void reset() { z1 = z2 = 0; }


protected:

  // State:
  TSig z1 = 0, z2 = 0;

  // Coeffs:
  TPar a0 = 0, a1 = 0, a2 = 0;  // Mixing coeffs - maybe rename to aL, aB, aH
  TPar g  = 0;                  // Integrator gain (?)
  
  TPar r  = 0;          // Damping (?)


  TPar gpr = 0;
  TPar scl = 1;

};

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setup(Mode mode, TPar w, TPar Q, TPar A)
{
  switch(mode)
  {
  case Mode::Bypass:        setupBypass();               break;
  case Mode::Lowpass:       setupLowpass(      w, Q);    break;
  case Mode::Highpass:      setupHighpass(     w, Q);    break;
  case Mode::BandpassSkirt: setupBandpassSkirt(w, Q);    break;
  case Mode::BandpassPeak:  setupBandpassPeak( w, Q);    break;
  case Mode::Bandstop:      setupBandstop(     w, Q);    break;
  case Mode::Allpass:       setupAllpass(      w, Q);    break;
  case Mode::Bell:          setupBell(         w, Q, A); break;
  case Mode::LowShelf:      setupLowShelf(     w, Q, A); break;
  case Mode::HighShelf:     setupHighShelf(    w, Q, A); break;
  default:
  {
    rsError("Unknown filter type in rsStateVariableFilterMystran::setup");
    a0 = a1 = a2 = 0;
    g = r = 0;
    gpr = 0;
    scl = 0;
  };
  }
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBypass()
{
  // H(s) = 1 

  g   = 0; 
  r   = 0;
  gpr = 0;
  scl = 1;
  a0  = 0;
  a1  = 0;
  a2  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupLowpass(TPar w, TPar Q)
{
  // H(s) = 1 / (s^2 + s/Q + 1)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 1;
  a1  = 0;
  a2  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupHighpass(TPar w, TPar Q)
{
  // H(s) = s^2 / (s^2 + s/Q + 1)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 0;
  a1  = 0;
  a2  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandpassSkirt(TPar w, TPar Q)
{
  // H(s) = s / (s^2 + s/Q + 1)  (constant skirt gain, peak gain = Q)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 0;
  a1  = 1;
  a2  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandpassPeak(TPar w, TPar Q)
{
  // H(s) = (s/Q) / (s^2 + s/Q + 1)      (constant 0 dB peak gain)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 0;
  a1  = r;
  a2  = 0;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBandstop(TPar w, TPar Q)
{
  // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 1;
  a1  = 0;
  a2  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupAllpass(TPar w, TPar Q)
{
  // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)

  g   = tan(0.5*w); 
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 1;
  a1  = -r;
  a2  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupBell(TPar w, TPar Q, TPar A)
{
  // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)

  g   = tan(0.5*w); 
  r   = 1/(Q*A);                  // P = Q*A, r = 1/P
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 1;
  a1  = A*A*r;                    // A^2 / P = A/Q = A^2 * r
  a2  = 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupLowShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)

  g   = tan(0.5*w) / sqrt(A);
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = A*A;
  a1  = A*r;
  a2  = 1;                        // High-freq gain should be one for a low-shelf.
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupHighShelf(TPar w, TPar Q, TPar A)
{
  // H(s) = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)

  g   = tan(0.5*w) * sqrt(A);
  r   = 1/Q;
  gpr = g + r;
  scl = 1 / (1 + g*gpr);
  a0  = 1;                        // Low-freq gain should be one for a low-shelf.
  a1  = A*r;
  a2  = A*A;
}


template<class TSig, class TPar>
TSig rsStateVariableFilterMystran<TSig, TPar>::getSample(TSig in)
{
  // Compute outputs:
  //TSig hp = (in - (g+r)*z1 - z2) / (1 + g*(g+r));
  TSig hp = (in - gpr*z1 - z2) * scl;
  TSig bp = z1 + g*hp; 
  TSig lp = z2 + g*bp;

  // State variable update:
  z1 = 2*bp - z1;                // Equivalent to: z1 += 2*g*hp
  z2 = 2*lp - z2;                // Equivalent to: z2 += 2*g*bp

  // Mix final output:
  return a2*hp + a1*bp + a0*lp;

  // ToDo:
  //
  // - Optimize: Precompute g+r and 1/(1+g*(g+r)) in the setup functions. We may then get rid of
  //   the r member variable
}







#endif