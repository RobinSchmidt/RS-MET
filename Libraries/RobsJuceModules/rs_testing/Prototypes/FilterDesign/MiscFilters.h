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

//=================================================================================================

/** Subclass of rsStateVariableFilterMystran that extends it by a general setup() function that
takes a mode parameter and then dispatches to the different setup functions. It also adds some
inquiry functions that retrieve the design parameter from the coefficients. That's probably a
rather useless functionality - if something like that is desired, it would make more sense to
just store it in additional member variables in some wrapper subclass. */

template<class TSig, class TPar>
class rsStateVariableFilterMystran2 : public rsStateVariableFilterMystran<TSig, TPar>
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

  void setup(Mode mode, TPar w, TPar Q, TPar A = TPar(1))
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
      g = 0;
      gpr = 0;
      scl = 0;
    };
    }
  }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  bool isLowpass()       const { return a0 == 1 && a1 == 0 && a2 == 0; }
  bool isHighpass()      const { return a0 == 0 && a1 == 0 && a2 == 1; }
  bool isBandpass()      const { return a0 == 0 && a1 >  0 && a2 == 0; }
  bool isBandpassSkirt() const { return a0 == 0 && a1 == 1 && a2 == 0; }
  bool isBandpassPeak()  const { return a0 == 0 && a1 != 1 && a2 == 0; }        // a1 = r = 1/Q
  bool isBandstop()      const { return a0 == 1 && a1 == 0 && a2 == 1; }
  bool isAllpass()       const { return a0 == 1 && a1 <  0 && a2 == 1; }        // a1 = -1/Q
  bool isBell()          const { return a0 == 1 && a1 >  0 && a2 == 1; }        // a1 = A^2/Q
  bool isShelf()         const { return isLowShelf() || isHighShelf(); }
  bool isLowShelf()  const { return /*a0 >  0 &&*/ a0 != 1 && a1 >  0 && a2 == 1; } // a0=A^2, a1=A/Q
  bool isHighShelf() const { return a0 == 1 && a1 >  0 && /* a2 >  0 && */ a2 != 1; } // a2=A^2, a1=A/Q



  // ToDo: check, if we really need the a0 > 0 condition for LS and a2 > 0 condition for HS

  // I think there's an edge case of Q = 1 where constant peak and constant skirt bandpasses are
  // indistinguishable. I think, in this case, both should return true - but they currently don't
  // I think. We need to check that in any case one and only one of them returns true, i.e. that
  // the conditions are disjoint or mutually exclusive - except in edge cases maybe.


  TPar getOmega() const
  {
    if(isLowShelf())
    {
      TPar r = gpr - g;
      TPar A = a1 / r;
      return 2*atan(g*sqrt(A));        // g = tan(w/2) / sqrt(A);
    }
    else if(isHighShelf())
    {
      TPar r = gpr - g;
      TPar A = a1 / r;
      return 2*atan(g/sqrt(A));        // g = tan(w/2) * sqrt(A)
    }
    else
      return 2*atan(g);                // g = tan(w/2)
  }

  TPar getQualityFactor() const
  {
    if(isBell())
    {
      return sqrt(1/(a1*(gpr-g)));

      // We have 3 equations involving r, Q, A: (1) r = 1/(Q*A), (2) gpr = g + r, (3) a1 = A^2 * r
      // where the knowns are  gpr, g, a1  and the unknowns are r, Q, A. These equations can be 
      // grabbed directly from the code in setupBell(). Sage can solve this simple nonlinear 
      // system of equations for us with the following code:
      //
      //  var("gpr g a1 r  Q A")
      //  e1 = r   == 1/(Q*A)
      //  e2 = gpr == g + r
      //  e3 = a1  == A^2 * r
      //  solve([e1,e2,e3],[r,Q,A])
      //
      // Picking the 1st solution and manually making it prettier gives the result above.
    }
    else
      return 1 / (gpr - g);  // gpr = g + r  ->  r = gpr - g = 1/Q
  }

  TPar getBellGain() const
  {
    rsAssert(isBell(), "Calling this function only makes sense for bell filters");
    TPar Q = sqrt(1/(a1*(gpr-g)));
    return 1 / (Q*(gpr-g));
  }

  TPar getShelfGain() const
  {
    rsAssert(isShelf(), "Calling this function only makes sense for shelf filters");
    return a1 / (gpr - g);
  }

  rsComplex<TPar> getTransferFunctionAtOld(const rsComplex<TPar>& z)
  {
    // We cheat here. We know that the filter has the s-domain trasfer function 
    // H(s) = (a0 + a1*s + a2*s^2) / (1 + s/Q + s^2) and we know that we need to substitute
    // s according to the bilinear transform as s = k * (z-1)/(z+1) where the scaling factor k is
    // given by 1/g which I figured out by trial and error (ToDo: give an explanation why it is 
    // that factor)

    //TPar w = getOmega();
    TPar Q = getQualityFactor();
    rsComplex<TPar> s = (1/g) * (z-TPar(1)) / (z+TPar(1));

    if(isBell())
    {
      rsError("This does not yet work. The formula is still wrong.");
      TPar A = getBellGain();
      //TPar P = A*Q;

      //return (a0 + a1*s + a2*s*s) / (TPar(1) + s/(Q*A) + s*s); // Nope! Wrong!
      return (a0 + a1*s + a2*s*s) / (TPar(1) + s/(Q/A) + s*s); // 

      //return (a0 + a1*s*A + a2*s*s) / (TPar(1) + s/Q + s*s);  // Wrong

      // Or maybe there were not wrong - it seems that the reference may have been wrong
    }

    return (a0 + a1*s + a2*s*s) / (TPar(1) + s/Q + s*s);

    // Someday, we want to have a proper implementation that directly computes H(z) in terms of
    // our coefficients without reconstructing the design parameters. I have not yet figured that
    // out, though.
  }






  // New - needs test:
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z)
  {
    ////// New:
    //TPar b0, b1, b2, a0, a1;
    //getBiquadCoeffs(&b0, &b1, &b2, &a1, &a2);
    //rsComplex<TPar> d  = TPar(1)/z, d2 = d*d;                  // d = z^-1, d2 = z^-2
    //rsComplex<TPar> H_ = (b0 + b1*d + b2*d2) / (TPar(1) + a1*d + a2*d2);
    ////// Does not work! Why?

    // Old:
    rsComplex<TPar> H_lp = getLowpassTransferFunctionAt(z);
    rsComplex<TPar> H_bp = getBandpassTransferFunctionAt(z);
    rsComplex<TPar> H_hp = getHighpassTransferFunctionAt(z);
    rsComplex<TPar> H    = a0*H_lp + a1*H_bp + a2*H_hp;
    return H;


  }

  rsComplex<TPar> getLowpassTransferFunctionAt(const rsComplex<TPar>& z)
  {
    TPar s  =  scl;
    TPar c  =  gpr;
    TPar b0 =  g*s;
    TPar b1 =  2*g*s;
    TPar b2 =  g*s;
    TPar a0 =  1/g;
    TPar a1 = -2/g + 2*c*s + 2*g*s;
    TPar a2 =  1/g - 2*c*s + 2*g*s;
    rsComplex<TPar> d = TPar(1)/z, d2 = d*d; // d = z^-1, d2 = z^-2
    rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (a0 + a1*d + a2*d2);
    return H;
  }

  rsComplex<TPar> getHighpassTransferFunctionAt(const rsComplex<TPar>& z)
  {
    TPar s  =  scl;
    TPar c  =  gpr;
    TPar b0 =  s;
    TPar b1 = -2*s;
    TPar b2 =  s;
    TPar a0 =  1;
    TPar a1 =  2*c*g*s + 2*g*g*s - 2;
    TPar a2 = -2*c*g*s + 2*g*g*s + 1;
    rsComplex<TPar> d = TPar(1)/z, d2 = d*d;
    rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (a0 + a1*d + a2*d2);
    return H;
  }

  rsComplex<TPar> getBandpassTransferFunctionAt(const rsComplex<TPar>& z)
  {
    TPar s  =  scl;
    TPar c  =  gpr;
    TPar b0 =  g*s;
    TPar b1 =  0;
    TPar b2 = -g*s;
    TPar a0 =  1;
    TPar a1 =  2*c*g*s + 2*g*g*s - 2;
    TPar a2 = -2*c*g*s + 2*g*g*s + 1;
    rsComplex<TPar> d = TPar(1)/z, d2 = d*d;
    rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (a0 + a1*d + a2*d2);
    return H;

    // Note: The a-coeffs are the same as in the highpass case.
  }




  // New: needs tests:


  void getBiquadCoeffs(TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2)
  {
    getBiquadNumeratorCoeffs(b0, b1, b2);
    getBiquadDenominatorCoeffs(  a1, a2);
  }
  // Maybe rename to toBiquad


  void getBiquadDenominatorCoeffs(TPar* a1, TPar* a2)
  {
    TPar s = scl;
    TPar c = gpr;
    *a1 =  2*(c*g + g*g)*s - 2;
    *a2 = -2*(c*g - g*g)*s + 1;
    // Simplify: factor out g, create variable for the common subexpression
  }

  void getBiquadNumeratorCoeffs(TPar* b0, TPar* b1, TPar* b2)
  {
    TPar t0, t1, t2;  // Temporaries

    getBiquadNumeratorCoeffsLP(&t0, &t1, &t2);
    *b0 = a0*t0;
    *b1 = a0*t1;
    *b2 = a0*t2;

    getBiquadNumeratorCoeffsBP(&t0, &t1, &t2);
    *b0 += a1*t0;
    *b1 += a1*t1;
    *b2 += a1*t2;

    getBiquadNumeratorCoeffsHP(&t0, &t1, &t2);
    *b0 += a2*t0;
    *b1 += a2*t1;
    *b2 += a2*t2;

    // We really should rename a0 to aL etc.
  }

  void getBiquadNumeratorCoeffsLP(TPar* b0, TPar* b1, TPar* b2)
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =   s*g*g;
    *b1 = 2*s*g*g;
    *b2 =   s*g*g;
  }

  void getBiquadNumeratorCoeffsBP(TPar* b0, TPar* b1, TPar* b2)
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =  g*s;
    *b1 =  0;
    *b2 = -g*s;
  }

  void getBiquadNumeratorCoeffsHP(TPar* b0, TPar* b1, TPar* b2)
  {
    TPar s = scl;
    TPar c = gpr;
    *b0 =  s;
    *b1 = -2*s;
    *b2 =  s;
  }





};

