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
class rsStateVariableFilter2 : public rsStateVariableFilter<TSig, TPar>  // Rename!
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
    case Mode::Lowpass:       setupLowpass(      w, Q   ); break;
    case Mode::Highpass:      setupHighpass(     w, Q   ); break;
    case Mode::BandpassSkirt: setupBandpassSkirt(w, Q   ); break;
    case Mode::BandpassPeak:  setupBandpassPeak( w, Q   ); break;
    case Mode::Bandstop:      setupBandstop(     w, Q   ); break;
    case Mode::Allpass:       setupAllpass(      w, Q   ); break;
    case Mode::Bell:          setupBell(         w, Q, A); break;
    case Mode::LowShelf:      setupLowShelf(     w, Q, A); break;
    case Mode::HighShelf:     setupHighShelf(    w, Q, A); break;
    default:
    {
      rsError("Unknown filter type in rsStateVariableFilterMystran::setup");
      aL = aB = aH = 0;
      g = 0;
      c = 0;
      s = 0;
    };
    }
  }

  /*
  void setupFromBiquad(TPar b0, TPar b1, TPar b2, TPar a1, TPar a2)
  {
    TPar T = (a1*a1 - a2*a2 - 2*a2 - 1);
    TPar r;

    if(T < 0)
    {
      TPar S = sqrt(-1/T);

      aB =  (2*b0*S - 2*b2*S);
      g  = -(1/((a1 - a2 - 1)*S));
      r  =  (2*(a2 - 1) / (T*S));
    }
    else
    {
      rsError("Not yet implemented correctly");

      TPar S = sqrt(1/T);

      aB = +(2*b0*S - 2*b2*S);                    // May need a minus
      g  = +(1/((a1 - a2 - 1)*S));                // May need a minus
      r  = +(2*(a2 - 1) / (T*S));                 // May need a minus

      // Tested and found to be wrong: ---,--+,-+-,-++,+--,+-+,++-,+++
      // Hmm - just inverting the sign of the argument of the square root and then compensating via
      // sign flips in the formulas for aB,g,r doesn't work
    }
    

    aH = -(b0 - b1 + b2) / (a1 - a2 - 1);
    aL =  (b0 + b1 + b2) / (a1 + a2 + 1);
    c  =  g + r;
    s  =  1 / (1 + g*c);

    // Sage gave me a second solution:
    //   aB == -2*b0*S + 2*b2*S
    //   g  ==  1/((a1 - a2 - 1)*S)
    //   r  == -2*(a2 - 1) / (T * S)
    // The rest is the same. Maybe sometimes we need that solution? Maybe when the b-coeffs are
    // negative?
  }
  */
  // ToDo: check what happens with the second solution. Maybe it's needed in certain cases when
  // the filter includes a sign inversion? Maybe make unit tests with random biquad coeffs and
  // make roundtrips. Maybe also try random SVF coeffs. I think, we can distinguish between the 
  // applicability of solution 1 or 2 by looking at the rhs in the computation of r. I think, r 
  // must always be positive - or does it? It is 1/Q and also 2*R where R is the damping. Maybe
  // a negative damping coeffs would lead to an unstable filter? Try it!
  // Can it happen that the argument of the square-root is negative, i.e. T is positive? Maybe in 
  // this case, we need to take the sqrt of the absolute value of T and then invert the signs of
  // aB, g, r?
  // Optimize the common subexpressions: 2*(a2 - 1), (a1 - a2 - 1). Maybe we can also take
  // the sqrt of T itself rather than 1/T and then compute r = S * 2*(a2 - 1) / T


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  bool isLowpass()       const { return aL == 1 && aB == 0 && aH == 0; }
  bool isHighpass()      const { return aL == 0 && aB == 0 && aH == 1; }
  bool isBandpass()      const { return aL == 0 && aB >  0 && aH == 0; }
  bool isBandpassSkirt() const { return aL == 0 && aB == 1 && aH == 0; }
  bool isBandpassPeak()  const { return aL == 0 && aB != 1 && aH == 0; }        // a1 = r = 1/Q
  bool isBandstop()      const { return aL == 1 && aB == 0 && aH == 1; }
  bool isAllpass()       const { return aL == 1 && aB <  0 && aH == 1; }        // a1 = -1/Q
  bool isBell()          const { return aL == 1 && aB >  0 && aH == 1; }        // a1 = A^2/Q
  bool isShelf()         const { return isLowShelf() || isHighShelf(); }
  bool isLowShelf()  const { return /*aL >  0 &&*/ aL != 1 && aB >  0 && aH == 1; } // a0=A^2, a1=A/Q
  bool isHighShelf() const { return aL == 1 && aB >  0 && /* aH >  0 && */ aH != 1; } // a2=A^2, a1=A/Q



  void getCoeffs_g_c_s(TPar* g, TPar* c, TPar* s) const
  {
    *g = this->g;
    *c = this->c;
    *s = this->s;
  }

  void getMixCoeffs_l_b_h(TPar* l, TPar* b, TPar* h) const
  {
    *l = this->aL;
    *b = this->aB;
    *h = this->aH;
  }


  // ToDo: check, if we really need the aL > 0 condition for LS and aH > 0 condition for HS

  // I think there's an edge case of Q = 1 where constant peak and constant skirt bandpasses are
  // indistinguishable. I think, in this case, both should return true - but they currently don't
  // I think. We need to check that in any case one and only one of them returns true, i.e. that
  // the conditions are disjoint or mutually exclusive - except in edge cases maybe.

  bool hasSameCoeffsAs(const rsStateVariableFilter2<TSig, TPar>& f, TPar tol)
  {
    if(!rsIsCloseTo(aL, f.aL, tol)) return false;
    if(!rsIsCloseTo(aB, f.aB, tol)) return false;
    if(!rsIsCloseTo(aH, f.aH, tol)) return false;
    if(!rsIsCloseTo(g,  f.g,  tol)) return false;
    if(!rsIsCloseTo(c,  f.c,  tol)) return false;
    if(!rsIsCloseTo(s,  f.s,  tol)) return false;
    return true;
  }

  TPar getOmega() const
  {
    if(isLowShelf())
    {
      TPar r = c - g;
      TPar A = aB / r;
      return 2*atan(g*sqrt(A));        // g = tan(w/2) / sqrt(A);
    }
    else if(isHighShelf())
    {
      TPar r = c - g;
      TPar A = aB / r;
      return 2*atan(g/sqrt(A));        // g = tan(w/2) * sqrt(A)
    }
    else
      return 2*atan(g);                // g = tan(w/2)
  }

  TPar getQualityFactor() const
  {
    if(isBell())
    {
      return sqrt(1/(aB*(c-g)));

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
      return 1 / (c - g);  // gpr = g + r  ->  r = gpr - g = 1/Q
  }

  TPar getBellGain() const
  {
    rsAssert(isBell(), "Calling this function only makes sense for bell filters");
    TPar Q = sqrt(1/(aB*(c-g)));
    return 1 / (Q*(c-g));
  }

  TPar getShelfGain() const
  {
    rsAssert(isShelf(), "Calling this function only makes sense for shelf filters");
    return aB / (c - g);
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
      return (aL + aB*s + aH*s*s) / (TPar(1) + s/(Q/A) + s*s); // 

      //return (a0 + a1*s*A + a2*s*s) / (TPar(1) + s/Q + s*s);  // Wrong

      // Or maybe there were not wrong - it seems that the reference may have been wrong
    }

    return (aL + aB*s + aH*s*s) / (TPar(1) + s/Q + s*s);

    // Someday, we want to have a proper implementation that directly computes H(z) in terms of
    // our coefficients without reconstructing the design parameters. I have not yet figured that
    // out, though.
  }


};

//=================================================================================================

/** A class for implementing digital filters in state space form. It implements the (vector/MIMO)
difference equation:

  y[n]   = C * x[n] + D * u[n]           output generation for sample index n
  x[n+1] = A * x[n] + B * u[n]           state update to prepare for the next sample

where x is a length N state vector, u is a length p input vector, y is a length q output vector, A
is an N-by-N state transition matrix, B is an N-by-p injection matrix, C is an q-by-N output matrix 
and D is a q-by-p feedthrough matrix. (verify sizes and terminology - I made some of them up. Edit:
on page 357, in a code comment JOS calls B,C,D the input, output, feed-around matrices 
respectively)

The system has the q-by-p MIMO transfer function matrix:

  H(z) = D + C * (z*I - A)^(-1) * B     (1) Eq G.5

The H(i,j) element of this matrix gives the transfer function from the j-th input to the i-th 
output [VERIFY!].


ToDo: 

- Give expression for impulse response (see (1) page 346)

- Add conversions from/to direct forms


References:

  (1) Introduction to Digital Filters with Audio Application (Julius O. Smith)

*/


template<class T>     // ToDo: have TSig and TPar
class rsStateSpaceFilter
{

public:

  /** Sets up the shapes of our matrices to allow for the desired number of ins/outs/states but 
  doesn't initialize the contents of those matrices. If you expect to change these sizes later 
  during realtime processing, this function can also be used once with the maximum epxected sizes
  to reserve enough memory to avoid need ing re-allocations later when the size changes dynamically
  during processing. */
  void setDimensions(int numIns, int numOuts, int numStates);

  /** Sets up our A,B,C,D matrices. This may reshape our member matrices and will copy the data 
  from the given argument matrices into them. Reshaping may re-allocate memory unless our member 
  matrices already have (more than) enough space, which you can ensure using setDimensions(). */
  void setup(const rsMatrixView<T>& A, const rsMatrixView<T>& B, const rsMatrixView<T>& C,
    const rsMatrixView<T>& D);
  // ToDo:
  // -Maybe we should keep only references or pointers to these matrices here? The way we do it now
  //  requires redundant existence of these matrices in memory which is bad. But maybe the state 
  //  vector x and its temporary storage t should nevertheless stay non-reference members?
  //  ...hmmm...

  /** Computes the transfer function matrix at the given complex number z. The (i,j)-th element
  of this matrix is the point-to-point transfer function from input j to output i. ...I think...
  or maybe it's the other way around? ...Figure out!  */
  rsMatrix<rsComplex<T>> getTransferFunctionAt(rsComplex<T> z);



  /** Processes a single MIMO output frame at a time. */
  void processFrame(T* ins, T* outs)
  {
    using AT = rsArrayTools;
    using MV = rsMatrixView<T>;

    // Wrap some matrix-view objects around inputs and output to view them as p-by-1 and q-by-1 
    // column-vectors:
    MV u(p, 1, ins);                             // u = input vector
    MV y(q, 1, outs);                            // y = output vector

    // Compute output y = C*x + D*u:
    MV::matrixMultiply(          C, x, &y);      // y = C*x
    MV::matrixMultiplyAccumulate(D, u, &y);      // y = C*x + D*u

    // Update the state x = A*x + B*u:
    MV::matrixMultiply(          A, x, &t);      // t = A*x
    MV::matrixMultiplyAccumulate(B, u, &t);      // t = A*x + B*u
    AT::copy(&t(0,0), &x(0,0), N);               // x = t
  }
  // Notes:
  // -Needs more tests, if it does the right thing.
  // -We currently have no safeguards against the client passing too short arrays - should we? We 
  //  already know, that ins should have length q and outs have length p and expect the caller to
  //  know and respect that, too. After all, the client code must have set us up that way at some 
  //  previous point via calling e.g. setup(). Safeguarding here would require otherwise useless 
  //  and redundant function parameters like numIns, numOuts. Not sure, if that's a good idea.
  //  ...but maybe...we'll see....
  // -Maybe try to shorten the code. instead of AT::copy(...) we could use something like
  //  x.copyDataFrom(t). But calling AT::copy may be more efficient because it bypasses the
  //  setShape() call in copyDataFrom(). The function names matrixMultiply/Accumulate could be 
  //  shortened to something like mul/Accum ...the matrix prefix is redundant because we already 
  //  are in class rsMatrixView


  /** Resets out internal state vector to all zeros. */
  void reset() { x.setToZero(T(0)); }
  // ToDo: 
  // -Maybe have a function setState(T* newState) that lets the client explicitly set up any 
  //  desired initial state/condition.


protected:


  int N = 0;  // number of internal states
  int p = 0;  // number of inputs
  int q = 0;  // number of outputs
  // These are actually redundant but convenient. They could be inferred from certain row- and 
  // column settings in our matrices below if saving that little amount of extra space seems 
  // worthwhile. Maybe rename them into numStates, numIns, numOuts.


  rsMatrix<T> x, t, A, B, C, D;
  // Meaning of those matrices
  //   x: state vector, N-by-1 column vector
  //   t: temporary storage for x during state update
  //   A: state transition matrix, N-by-N matrix
  //   B: injection matrix, N-by-p matrix (verify!)
  //   C: output matrix, q-by-N matrix (verify!)
  //   D: passthrough matrix, q-by-p (verify!)
  // We also use the notation:
  //   u: input vector, p-by-1
  //   y: output vector, q-by-1
  // but we don't need any class members for these I/O variables. See (1) pg 345, Appendix G.

  // ToDo: 
  // -Check terminology. I made some of it up myself (feedthrough, injection). The book (1) 
  //  doesn't give them special names.
  // -Perhaps production code should use sparse matrices? I think, the state update matrices are
  //  typically sparse, right? But what about the other matrices? Are they also typically sparse?
  //  Maybe only A should be sparse but B,C,D dense? ...figure out!
  // -Implement a getTransferFunction() function that returns an rsMatrix of type
  //  rsRationalFunction

};

template<class T> 
void rsStateSpaceFilter<T>::setDimensions(int numIns, int numOuts, int numStates)
{
  p = numIns;
  q = numOuts;
  N = numStates;
  x.setShape(N, 1);
  t.setShape(N, 1);
  A.setShape(N, N);
  B.setShape(N, p);
  C.setShape(q, N);
  D.setShape(q, p);
  reset();
}

template<class T> 
void rsStateSpaceFilter<T>::setup(const rsMatrixView<T>& newA, const rsMatrixView<T>& newB,
  const rsMatrixView<T>& newC, const rsMatrixView<T>& newD)
{
  // Retrieve and set up desired dimensions:
  N = newA.getNumRows();     // number of states
  p = newB.getNumColumns();  // number of inputs
  q = newC.getNumRows();     // number of outputs
  setDimensions(p, q, N);
  // ToDo: Verify and document that no allocations take place here, when already enough memory was
  // allcoated previously. See rsMatrix::setShape - it calls resize on a std::vector which should
  // reallocate only in case of growth.

  // Perform some sanity checks on the input matrices:
  rsAssert(newA.isSquare(), "State transition matrices must be square");
  // ...more checks to come: Make sure, that all the desired relations between the shapes of the 
  // given  matrices are satisfied. Maybe factor these checks out into a function checkSanity()
  // or something.

  // Copy the new matrix data into our members:
  A.copyDataFrom(newA);
  B.copyDataFrom(newB);
  C.copyDataFrom(newC);
  D.copyDataFrom(newD);
  // Hmm...copyDataFrom also calls setShape. These calls are redundant with those in setDimensions.
  // Maybe it doesn't matter but perhaps it would be nicer to avoid it...we'll see...
  // Maybe we should just keep references to some A,B,C,D matrices owned by cleint code anyway.
  // That avoids redundancies and makes it easier to implement time-variant operation. Client code
  // could just vary the matrices. Maybe we should have a 2-level API. A lower level that avoids
  // redundancies and a higher level for convenience. Then, one can start with the high-level API
  // and optimize later.

  //rsError("Not yet implemented");
}

template<class T> 
rsMatrix<rsComplex<T>> rsStateSpaceFilter<T>::getTransferFunctionAt(rsComplex<T> z)
{
  using Comp = rsComplex<T>;
  using Mat  = rsMatrix<Comp>;

  // Complexify our matrices:
  Mat Ac = rsConvert<T, Comp>(A);
  Mat Bc = rsConvert<T, Comp>(B);
  Mat Cc = rsConvert<T, Comp>(C);
  Mat Dc = rsConvert<T, Comp>(D);

  // Form the matrix M = (z*I - A)^(-1):
  Mat M = -Ac;                                           // M =      - A
  for(int i = 0; i < M.getNumRows(); i++) M(i,i) += z;   // M = (z*I - A)
  M = rsLinearAlgebraNew::inverse(M);                    // M = (z*I - A)^(-1)

  // Evaluate the transfer function:
  Mat H = Dc + Cc*M*Bc;                                  // H(z) = D + C * (z*I - A)^(-1) * B 
  return H;

  // This can probably be optimized a lot. Firstly, we may want to have a lower level API that lets
  // the caller pass pre-allocated pointers to rsMatrixView to get around the internal memory 
  // allocations. Secondly, maybe we can get away without the inversion and formulate it as a 
  // solution to a linear system? Thirdly, not all matrices need to be complex. We just do it that
  // way because the operators +,* need matrices of the same element type so we just complexify all
  // our matrices.
}


/** Converts a state variable filter with coefficients g,c,s into a state space form with one input
and 3 outputs for the highpass, bandpass and lowpass part. The matrices for the state space form 
are given by:

  A = [1-2gsc     -2gs ]     B = [2gs ]     C = [-sc        -s]     D = [  s]
      [2g-2ggsc  1-2ggs]         [2ggs]         [1-gsc     -gs]         [ gs]
                                                [g-ggsc  1-ggs]         [ggs]  

See Notes/StateVariableFilter.txt for derivation of these formulas. The caller needs to pass 
pre-allocated objects of type rsMatrixView of the right shapes. The required shapes are: 
A: 2x2, B: 2x1, C: 3x2, D: 3x1. */
template<class T>
void stateVariableToStateSpace(T g, T c, T s, 
  rsMatrixView<T>* A, rsMatrixView<T>* B, rsMatrixView<T>* C, rsMatrixView<T>* D)
{
  T gs  = g*s;
  T ggs = g*gs;

  // 2x2 state transition matrix:
  rsAssert(A->hasShape(2,2));
  (*A)(0,0) =  1  -2*gs*c;
  (*A)(0,1) =     -2*gs;
  (*A)(1,0) =  2*g-2*ggs*c;
  (*A)(1,1) =  1  -2*ggs;

  // 2x1 input matrix:
  rsAssert(B->hasShape(2,1));
  (*B)(0,0) =  2*gs;
  (*B)(1,0) =  2*ggs;

  // 3x2 output matrix:
  rsAssert(C->hasShape(3,2));
  (*C)(0,0) = -s*c;
  (*C)(0,1) = -s;
  (*C)(1,0) = 1-gs*c;
  (*C)(1,1) =  -gs;
  (*C)(2,0) = g-ggs*c;
  (*C)(2,1) = 1-ggs;

  // 3x1 feedaround matrix:
  rsAssert(D->hasShape(3,1));
  (*D)(0,0) =   s;
  (*D)(1,0) =  gs;
  (*D)(2,0) = ggs;
}
// Maybe make this a static member function of rsStateSpaceFilter. Rename it to 
// fromStateVariableFilter or something. ...or maybe the name is ok. Dunno. Or maybe it should
// go into class rsFilterCoefficientConverter


/** Converts a direct form filter into state space form. The formulas are taken from Julius Smith's
book about filters, pages 351-352. ...TBC... */
template<class T>
void directFormToStateSpace(std::vector<T> b, std::vector<T> a,
  rsMatrix<T>* A, rsMatrix<T>* B, rsMatrix<T>* C, rsMatrix<T>* D)
{
  using Vec = std::vector<T>;

  // Normalize to a[0] = 1:
  if(a[0] != T(1))
  {
    T s = T(1) / a[0];
    rsScale(b, s);
    rsScale(a, s);
  }

  // Zero-pad the shorter array, if needed:
  rsPadToSameSize(b, a, T(0));

  // Compute some intermediate variables:
  int N = b.size() - 1;     // Filter order
  T   b0 = b[0];
  Vec beta(b.size());
  beta[0] = 0;              // Not used
  for(int k = 1; k < beta.size(); k++)
    beta[k] = b[k] - b0*a[k];

  // Compute SSF matrices:
  A->setShape(N,N); A->setToZero(b0);
  B->setShape(N,1); B->setToZero(b0);
  C->setShape(1,N); C->setToZero(b0);
  D->setShape(1,1); D->setToZero(b0);
  for(int i = 1; i <= N; i++)
    (*A)(0,i-1) = -a[i];
  for(int i = 1; i < N; i++)
    (*A)(i,i-1) = 1;
  (*B)(0,0) = 1;
  for(int i = 1; i <= N; i++)
    (*C)(0,i-1) = beta[i];
  (*D)(0,0) = b0;


  // Notes and ToDo:
  //
  //
  // - We deliberately pass b,a by value rather than by const reference because we may modify 
  //   them here.
  //
  // - Maybe try to avoid creation of the beta array. We can use the formula directly in the
  //   loop that assigns the C-matrix.
  //
  // - Create unit tests that also test some edge cases like empty a and/or b, a[0] != 1, etc.
}


//=================================================================================================


/** A class for representing (univariate) monomials, i.e. expressions of the form  c * x^p  for
some coefficient c and integer power (or exponent) p.

See: https://en.wikipedia.org/wiki/Monomial   */

template<class T> 
class rsMonomial
{

public:

  explicit rsMonomial(T newCoeff = T(0), int newPower = 0) : coeff(newCoeff), power(newPower) { }
  // Marked as explicit because we want to avoid hidden automatic conversions from type T


  void setup(T newCoeff, int newPower)
  {
    coeff = newCoeff;
    power = newPower;
  }

  void setCoeff(T newCoeff)
  {
    coeff = newCoeff;
  }

  void setPower(int newPower)
  {
    power = newPower;
  }




  T getCoeff() const { return coeff; }

  int getPower() const { return power; }



  T evaluateAt(T x) const { return coeff * rsPow(x, T(power)); }
  // Preliminary. We may want to use rsPowInt for integer exponents. That may be more efficient.
  // Here, we explicitly first convert the exponent to type T and then call rsPow(T x, T y).
  // But currently rsPowInt is not suitably defined. It expects two unsigned ints.

  //T evaluateAt(T x) const { return coeff * rsPowInt(x, power); }
  // rsPowInt is defined for x and power being integers. We really need a function where the
  // base is an arbitrary type and the epxonent is an integer

protected:

  T   coeff = T(0);
  int power = T(0);

};

///** Function to compare two monomials for a less-than relation that is defined by comparing
//the powers only. Monomials with the same power but with different coefficients are considered
//equivalent by this relation, i.e. if neither  lhs < rhs  nor  rhs < lhs  via inspecting the 
//powers only, the terms are considered equivalent. This kind of less-than relation is needed to sort
//the terms in rsSparsePolynomial to bring it into a canonical representation. */
//template<class T>
//bool rsLessByPower(const rsMonomial<T>& lhs, const rsMonomial<T>& rhs)
//{
//  if(lhs.getPower() < rhs.getPower())
//    return true;
//  return false;
//}
//// Needs tests



//=================================================================================================

/** A class for representing sparse polynomials, i.e. polynomials that have many zero coefficients.
We represent such sparse polynomials basically as a std::vector of monomials. */

template<class T>
class rsSparsePolynomial
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Default constructor. Constructs an empty sparse polynomial. This represents, by definition,
  the zero polynomial. */
  rsSparsePolynomial() {}

  /** Creates a polynomial from an initializer list for the terms. */
  rsSparsePolynomial(std::initializer_list<rsMonomial<T>> initList) : terms(initList) {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Reserves memory for the given number of terms. Can be called before calling functions like 
  addTerm() to pre-allocate the desired amount of memory beforehand when multiple terms are being
  added in a sequence. */
  void reserve(size_t numTerms) { terms.reserve(numTerms); }

  /** Clears the array of terms. */
  void clear() { terms.clear(); }

  /** Sets up the polynomial from a dense arrays of polynomial coeffs. When a coefficient in the 
  dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<T>& newCoeffs, T tol);


  /** Appends a term with given coeff and power to the end of our terms array. Beware that this 
  may decanonicalize the representation. */
  void appendTerm(T coeff, int power) { terms.emplace_back(rsMonomial<T>(coeff, power)); } 
  

  //void addTerm(const rsMonomial<T>& newTerm, T tol)
  //{ addTerm(newTerm.getCoeff(), newTerm.getPower(), tol); }


  void addTerm(T coeff, int power, T tol);


  // ToDo: write a function addTerm that also adds the term but maintains a canonical 
  // representation by scanning through the existing coeffs to try to find a term with same 
  // exponent. If one is found, add the coeff. If none is found, insert the coeff/power pair at 
  // the right position.




  /** Sets the number of terms. If the new number is less than the current number, it will just 
  cut off terms from the end. If the new number is greater than the current number, it will just
  extend our vector of terms and the added terms at the end are uninitialized, i.e. may contain 
  garbage. This function should only be used if you intend to set up the new terms via e.g. 
  setTerm() after calling setNumTerms(). So, it's a function that needs a lot of care to be used
  properly. */
  void setNumTerms(int newNumTerms) { terms.resize(newNumTerms); }
  // This may put the terms array into a non-canonical (or even invalid) state! Maybe it shouldn't
  // be used. We'll see....


  void setTerm(int index, T coeff, int power) 
  { rsAssert(isValidIndex(index));  terms[index].setup(coeff, power); }

  void setPower(int index, int newPower)
  { rsAssert(isValidIndex(index)); terms[index].setPower(newPower); }

  void setCoeff(int index, T newCoeff)
  { rsAssert(isValidIndex(index)); terms[index].setCoeff(newCoeff); }

  void scaleCoeff(int index, T scaler) { setCoeff(index, scaler * getCoeff(index)); }

  void scale(T scaler)
  {
    for(int i = 0; i < getNumTerms(); i++)
      scaleCoeff(i, scaler);
  }
  // rename to scaleCoeffs for consistency


  void shiftCoeff(int index, T amount) { setCoeff(index, amount + getCoeff(index)); }

  void shiftCoeffs(T amount)
  {
    for(int i = 0; i < getNumTerms(); i++)
      shiftCoeff(i, amount);
  }


  void shiftPower(int index, int amount) { setPower(index, amount + getPower(index)); }

  void shiftPowers(int amount)
  {
    for(int i = 0; i < getNumTerms(); i++)
      shiftPower(i, amount);
  }






  /** Reverses the array of terms. */
  void reverse() { rsReverse(terms); }

  /** Turns the representation of the polynomial into a canonical one. A canonical representation 
  has the following properties: (1) The powers are strictly increasing as function of index. 
  (2) No power appears more than once. (3) No zero coefficient appear. We achieve this by 
  first sorting the terms, then consolidating multiple terms with equal exponents into single
  terms and finally deleting all terms that have a coefficient zero (up to the given tolerance). */
  void canonicalize(T tol);

  void copyDataFrom(const rsSparsePolynomial<T>& other)
  {
    setNumTerms(other.getNumTerms());
    for(int i = 0; i < getNumTerms(); i++)
      setTerm(i, other.getCoeff(i), other.getPower(i));
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, iff this sparse polynomial is empty, i.e. has no terms. */
  bool isEmpty() const { return terms.empty(); }

  bool isZero(T tol) const
  {
    for(int i = 0; i < getNumTerms(); i++)
      if( rsAbs(getCoeff(i)) > tol )
        return false;
    return true;
  }

  /** Return true, iff the given index is valid, i.e. the object has a term with given index. */
  bool isValidIndex(int i) const { return i >= 0 && i < getNumTerms(); }

  /** Returns the number of terms in this polynomial. The i-th term is a monomial of the form 
  ci * x^pi with a coefficient ci and a power/exponent pi. */
  int getNumTerms() const { return (int) terms.size(); }

  /** Returns the minimum power that occurs in this polynomial. */
  int getMinPower() const;

  /** Returns the maximum power that occurs in this polynomial. In mathematical jargon, the 
  highest power in a polynomial is also known as the degree or order of the polynomial. */
  int getMaxPower() const;

  /** Returns the index of the maximum power. */
  int getMaxPowerIndex() const;

  /** Returns the degree of the polynomial. This is mathematical term for the term with the
  highest power/exponent that has a nonzero coefficient. */
  int getDegree() const { return getMaxPower(); }
  // This is basically an alias name for getMaxPower(). I'm not sure, if it's a good idea to have 
  // two functions that do the exact same thing. Maybe get rid of it. But on the other hand, it's 
  // nice to have to be consistent with the API of class rsPolynomial. 

  /** Returns the leading coefficient, i.e. the coefficient that multiplies the highest power of
  the input variable x. */
  T getLeadingCoeff() const;

  /** Returns the term (i.e. the monomial) at the given index. */
  rsMonomial<T> getTerm(int index) const { rsAssert(isValidIndex(index)); return terms[index]; }

  /** Returns the coefficient of the term with given index. */
  T getCoeff(int index) const {  rsAssert(isValidIndex(index)); return terms[index].getCoeff(); }

  /** Returns the power of the term with given index. */
  int getPower(int index) const { rsAssert(isValidIndex(index)); return terms[index].getPower(); }

  /** Checks if this sparse polynomial is in canonical representation. A representation is 
  canonical if it has no zero coefficients (up to a given tolerance) and if the powers are strictly
  increasing (as function of term-index). The empty polynomial is also accepted as a canonical 
  representation. It represents the zero polynomial. */
  bool isCanonical(T tol = T(0)) const;


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Evaluates the polynomial at the given x and returns the result. */
  T evaluateAt(T x) const;


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the polynomial at the given input x. */
  T operator()(T x) const { return evaluateAt(x); }

  /** Adds two polynomials. */
  rsSparsePolynomial<T> operator+(const rsSparsePolynomial<T>& q) const 
  { rsSparsePolynomial<T> r; add(*this, q, &r, T(0)); return r; }

  /** Subtracts two polynomials. */
  rsSparsePolynomial<T> operator-(const rsSparsePolynomial<T>& q) const 
  { rsSparsePolynomial<T> r; subtract(*this, q, &r, T(0)); return r; }

  /** Multiplies two polynomials. */
  rsSparsePolynomial<T> operator*(const rsSparsePolynomial<T>& q) const 
  { rsSparsePolynomial<T> r; multiply(*this, q, &r, T(0)); return r; }


  //-----------------------------------------------------------------------------------------------
  /** \name Low Level API. These functions operate on pre-allocated output parameters (passed by 
  pointer) which potentially avoids heap allocations. */

  static void add(
    const rsSparsePolynomial<T>& p,
    const rsSparsePolynomial<T>& q,
    rsSparsePolynomial<T>* r, T tol);

  static void subtract(
    const rsSparsePolynomial<T>& p,
    const rsSparsePolynomial<T>& q,
    rsSparsePolynomial<T>* r, T tol);

  static void weightedSum(
    const rsSparsePolynomial<T>& p, T wp,
    const rsSparsePolynomial<T>& q, T wq,
    rsSparsePolynomial<T>* r, T tol);

  static void multiply(
    const rsSparsePolynomial<T>& p,
    const rsSparsePolynomial<T>& q,
    rsSparsePolynomial<T>* r, T tol);


protected:

  std::vector<rsMonomial<T>> terms;

};


template<class T>
void rsSparsePolynomial<T>::setupFromDenseCoeffs(const std::vector<T>& newCoeffs, T tol)
{
  terms.clear();
  terms.reserve(newCoeffs.size());
  for(int i = 0; i < (int) newCoeffs.size(); i++)
    if(rsAbs(newCoeffs[i]) > tol)
      terms.emplace_back(rsMonomial<T>(newCoeffs[i], i));

  //canonicalize(); // Not sure, if we should do this automatically...maybe not

  // It's really important to use  >  rather than  >=  in the conditional. Consider tol = 0. If we
  // would use  >=  then  >= 0  would return true when the coeff is zero, so zero coeffs would get 
  // accepted which is not what we want.
}

template<class T>
void rsSparsePolynomial<T>::addTerm(T coeff, int power, T tol)
{
  rsAssert(isCanonical());

  int i = 0;
  while(i < getNumTerms())
  {
    if(getPower(i) == power)
    {
      shiftCoeff(i, coeff);
      if(rsAbs(getCoeff(i)) <= tol)
        rsRemove(terms, (size_t) i);
      return;
    }
    else if(getPower(i) < power)
    {
      i++;
    }
    else
    {
      break;
    }
  }
  rsInsert(terms, rsMonomial<T>(coeff, power), (size_t) i);
}

template<class T>
void rsSparsePolynomial<T>::canonicalize(T tol)
{
  // In the empty case, we have nothing to and we really *need* to return early in order to not 
  // get an access violation in the code below (in the  int p = getPower(0);  line):
  if(isEmpty())
    return;

  // Sort the terms by power/exponent:
  using Mon = rsMonomial<T>;
  std::sort(terms.begin(), terms.end(), 
            [](const Mon& lhs, const Mon& rhs){ return lhs.getPower() < rhs.getPower(); });

  // Consolidate multiple terms with equal power/exponent into single term: 
  int numTerms = getNumTerms();
  int p = getPower(0);              // Current power
  int r = 1;                        // Read index
  int w = 0;                        // Write index
  while(r < numTerms) 
  {
    if(getPower(r) == p)
      shiftCoeff(w, getCoeff(r));
    else 
    {
      w++;
      setTerm(w, getCoeff(r), getPower(r));
      p = getPower(r);
    }
    r++;
  }
  setNumTerms(w+1);
  // This algorithm works only when the terms are sorted by exponent so it doesn't really make 
  // sense to factor it out into a function in its own right. Doing so could invite calling it on 
  // unsorted term arrays in which case we would have a bug.

  // Remove terms with coefficient zero:
  rsRemoveIf(terms, [&tol](const Mon& term){ return rsAbs(term.getCoeff()) <= tol; });

  // Sanity check in debug mode:
  rsAssert(isCanonical());


  // ToDo:
  //
  // - Maybe try using  rsHeapSort(&terms[0], (int) terms.size(), &rsLessByPower);  instead of
  //   std::sort(..)
}



template<class T>
int rsSparsePolynomial<T>::getMinPower() const
{
  if(isEmpty())
    return 0;
  int minPower = std::numeric_limits<int>::max();
  for(auto& term : terms)
    minPower = rsMin(minPower, term.getPower());
  return minPower;
}

template<class T>
int rsSparsePolynomial<T>::getMaxPower() const
{
  if(isEmpty())
    return 0;
  int maxPower = std::numeric_limits<int>::min();
  for(auto& term : terms)
    maxPower = rsMax(maxPower, term.getPower());
  return maxPower;

  // The implementation is written in such a way that it should still work reasonably when the
  // client code sets up terms with negative powers. The empty polynomial will still have a max
  // power (aka degree) of zero. ...TBC...
}


template<class T>
int rsSparsePolynomial<T>::getMaxPowerIndex() const
{
  rsAssert(isCanonical());
  // The output of this function is not well defined when there are multiple terms with the highest
  // power, so this function should really only be used on canonical representations.

  if(isEmpty())
    return -1;

  int maxIndex = 0;
  int maxPower = getPower(0);
  for(int i = 1; i < getNumTerms(); i++)
  {
    if(getPower(i) > maxPower)
    {
      maxPower = getPower(i);
      maxIndex = i;
    }
  }

  return maxIndex;
}

template<class T>
T rsSparsePolynomial<T>::getLeadingCoeff() const 
{ 
  int i = getMaxPowerIndex();
  if(i != -1)
    return getCoeff(i);
  else
    return 0;
}
// Needs test.

template<class T>
bool rsSparsePolynomial<T>::isCanonical(T tol = T(0)) const
{
  // An empty polynomial is the canonical representation of the zero polynomial:
  if(isEmpty())
    return true;

  // Check that 0-th coeff is nonzero:
  if(rsAbs(getCoeff(0)) <= tol)
    return false;

  // Check that all other coeffs are also nonzero and that the powers are strictly increasing:
  int prevPow = getPower(0);                     // Previous power
  for(int i = 1; i < getNumTerms(); i++)
  {
    // Coeffs should be nonzero:
    if(rsAbs(getCoeff(i)) <= tol)
      return false;

    // Powers should be strictly increasing:
    int curPow = getPower(i);                    // Current power..
    if(curPow <= prevPow)
      return false;
    prevPow = curPow;                            // ..becomes previous power for next iteration.
  }

  return true;
}

template<class T>
T rsSparsePolynomial<T>::evaluateAt(T x) const 
{ 
  T y = 0;
  for(auto& term : terms)
    y += term.evaluateAt(x);
  return y;
}



template<class T>
void rsSparsePolynomial<T>::add(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();      // Number of terms in left operand p
  int Nq = q.getNumTerms();      // Number of terms in right operand q
  int Nr = Np + Nq;              // Number of terms in result r (before canonicalization)

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::subtract(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, -q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::weightedSum(
  const rsSparsePolynomial<T>& p, T wp,
  const rsSparsePolynomial<T>& q, T wq,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, wp * p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, wq * q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::multiply(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np * Nq;

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    for(int j = 0; j < Nq; j++)
      r->setTerm(i*Nq+j, p.getCoeff(i) * q.getCoeff(j), p.getPower(i) + q.getPower(j));

  r->canonicalize(tol);
}

// Rename to rsDivMod:
template<class T>
void rsDivide(
  const rsSparsePolynomial<T>& num,
  const rsSparsePolynomial<T>& den,
  rsSparsePolynomial<T>* quot,
  rsSparsePolynomial<T>* rem,
  T tol)
{
  rsAssert(!den.isZero(tol));

  rsSparsePolynomial<T> tmp1, tmp2;    // ToDo: Let the caller pass pre-allocated objects
  quot->clear();                       // quot = 0
  rem->copyDataFrom(num);              // rem  = num

  while(!rem->isZero(tol) && rem->getDegree() >= den.getDegree())
  {
    // t = lead(r) / lead(d):
    int iRem = rem->getMaxPowerIndex();
    int iDen = den. getMaxPowerIndex();
    T   cRem = rem->getCoeff(iRem);
    T   cDen = den. getCoeff(iDen);
    int pRem = rem->getPower(iRem);
    int pDen = den. getPower(iDen);
    T cT = cRem / cDen;
    T pT = pRem - pDen;

    // q = q + t:
    quot->addTerm(cT, pT, tol);

    // r = r - t * d:
    tmp1.copyDataFrom(*rem);
    tmp2.copyDataFrom( den);
    tmp2.scale(-cT);
    tmp2.shiftPowers(pT);
    rsSparsePolynomial<T>::add(tmp1, tmp2, rem, tol);
    // Maybe instead of using the two temp polynomials, use rem->addTerm in a loop over the terms
    // of den. But I'm not sure, if that's really better. The addTerm calls may trigger a lot of
    // data movement, too. Maybe try both variants and do benchmarks.
  }

  // The algorithm has been adapted from: 
  //
  //   https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
  //
  // In the following pseudocode, all variables (n,d,q,r,t) are polynomials (t is actually a 
  // monomial, though). Wikipedia says:
  //
  // Inputs:    n: numerator, d: denominator
  // Outputs:   q: quotient,  r: remainder
  // Require:   d != 0
  // Invariant: n = d * q + r                  # This holds at each step
  //
  // q = 0                                     # Init quotient to zero
  // r = n                                     # Init remainder to numerator
  // while( r != 0 and deg(r) >= deg(d) )
  // {
  //    t = lead(r) / lead(d)                  # t is a monomial
  //    q = q + t
  //    r = r - t * d
  // }
  // return (q, r)
}
// Needs more tests! It has already passed its first test, though.


template<class T>
bool rsIsCloseTo(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q, T tol)
{
  if(p.getNumTerms() != q.getNumTerms())
    return false;

  for(int i = 0; i < p.getNumTerms(); i++)
  {
    if(p.getPower(i) != q.getPower(i))
      return false;
    if(rsAbs(p.getCoeff(i) - q.getCoeff(i)) > tol)
      return false;
  }

  return true;
}





// ToDo:
//
// - Implement division with remainder. See:
//   https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
//   https://de.wikipedia.org/wiki/Polynomdivision#Algorithmus
//
// - Figure out what happens if client code uses negative powers. Currently, there's nothing that
//   prevents this and maybe it could even make sense to allow it. But then the notion of degree
//   gets murky. Maybe then there is indeed a difference between the degree and the max power in
//   the case of an empty polynomial? Maybe, for the time being, we should trap attempts to set up
//   terms with negative powers. This can later be relaxed, if needed.
//
// - Maybe keep the invariant that the polynomial is in canonical representation. Implementing 
//   algorithms for both cases is a mess. Maybe prepend a __ to those member functions that could
//   destroy the canonical representation to signal to the caller that they are now doing something
//   low level and potentially dangerous - like __shiftPower(int index, int amount). The regular
//   shiftPower function can still be present. It would just call __shiftPower() and then
//   canonicalize()
//
//
// Notes:
//
// - It might be tempting to write a constructor and/or setup function that takes a dense 
//   polynomial, i.e. an object of type rsPolynomial<T>. But I think, that's not a good idea 
//   because it would introduce unnecessary coupling.


//=================================================================================================

/** A class for representing sparse filters in direct form. We represent them using two sparse
polynomials. One for the numerator and one for the denominator of the transer function. The filter
is implemented in direct form 2 using a single delayline. */

template<class TSig, class TPar>
class rsSparseFilter
{

public:



  /** Sets up the filter from dense arrays of numerator and denominator coeffs. When a coefficient
  in the dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<TPar>& numCoeffs, 
    const std::vector<TPar>& denCoeffs, TPar tol)
  {
    num.setupFromDenseCoeffs(numCoeffs, tol);
    den.setupFromDenseCoeffs(denCoeffs, tol);
    updateDelayLineLength();
  }
  // This may allocate!
  // ToDo: use num.setupFromDenseCoeffs(&numCoeffs[0], (int) numCoeffs.size(), tol)


  void setNumNumeratorTerms(int newNumTerms)
  {
    num.setNumTerms(newNumTerms);
  }
  // This may allocate!

  void setNumDenominatorTerms(int newNumTerms)
  {
    den.setNumTerms(newNumTerms);
  }
  // This may allocate!

  void setNumeratorTerm(int index, TPar coeff, int delay)
  {
    num.setTerm(index, coeff, delay);
  }

  void setDenominatorTerm(int index, TPar coeff, int delay)
  {
    den.setTerm(index, coeff, delay);
  }
  // Actually, we really should call updateDelayLineLength() after setting a term because it 
  // potentially requires a change of the length. But: updateDelayLineLength() is expensive and 
  // setting terms is an operation that might be called in a loop or sequence in which case only
  // one update after the sequence of calls should be done. Maybe when calling it in a sequence,
  // we should use other functions like setNumeratorTermSuppressDelayUpdate. Or maybe give the 
  // function a boolean parameter updateDelayLength which defaults tor true


  /** Updates the length of the delayline according to the maximum power of z^-1 that occurs in
  numerator and denominator polynomial. */
  void updateDelayLineLength()
  {
    int maxDegree = rsMax(num.getDegree(), den.getDegree());
    delayLine.setMaximumDelayInSamples(maxDegree);
    delayLine.setDelayInSamples(maxDegree);
  }




  /** Applies a scaling factor to the filter. This basically means to scale all numerator coeffs by
  that factor. */
  void scale(TPar scaler)
  {
    num.scale(scaler);
  }

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given 
  amount. */
  void addPreDelay(int amountInSamples)
  {
    if(amountInSamples < 0)
    {
      rsError("Negative predelay is not allowed");
      return;
      // We could allow it if we already have some predelay and the given amount would just 
      // reduce it. Maybe we can relax the restriction to amountInSamples >= -num.getMinPower() or
      // something. If the minimum power is 3, we could allow a predelay amount of -3. This would 
      // then just reduce the predelay to zero.
    }

    num.shiftPowers(amountInSamples);
    updateDelayLineLength();
  }

  void removePreDelay()
  {
    int preDelay = num.getPower(0);
    // We assume here the the 0-th term is the one with the lowest power! This invariant should be
    // checked in isFilterValid().


    num.shiftPowers(-preDelay);
    updateDelayLineLength();
  }

  /** Turns the filter into its inverse. This basically amounts to swapping numerator and 
  denominator and possibly applying some scaling of the coefficients if b0 != 1. */
  void invert()
  {
    rsAssert(isFilterValid());

    removePreDelay();

    rsAssert(num.getPower(0) == 0);
    rsAssert(num.getCoeff(0) != 0);
    // We assume here that the 0-th num coefficient is the one that multiplies z^0 in the transfer
    // function and scales x[n] in the difference equation. Maybe we can relax that later to allow
    // the z^0 coeff to appear at a different position in the array

    TPar s = TPar(1) / num.getCoeff(0);
    scale(s);
    rsSwap(num, den);
    scale(s);

    // Figure out if rsSwap causes memory allocations when swapping the underlying std::vectors. 
    // Actually, swapping two vectors would only require pointer adjustments under the hood. Maybe
    // std::swap is clever enough to implement it that way?

    // When the filter has an initial delay, i.e. the first power in the numerator is not equal to 
    // zero, then we actually cannot invert the filter. A delay cannot be undone (at least not in
    // realtime). In this case, the best we can do is to invert the filter up to a delay. I think, 
    // we can do this by first figuring out the minimum exponent of the numerator and the 
    // subtracting that from all the numerator exponents. After that, we can invert as usual.
  }

  void reflectZeros()
  {
    int deg = num.getDegree();
    for(int i = 0; i < num.getNumTerms(); i++)
      num.setPower(i, deg - num.getPower(i));

    num.reverse();  // To make the array ordered by ascending powers
  }
  // not yet tested

  // ToDo: reflectPoles/reflectZeros - shoulv reverse the coeff arrays. i.e. the powers should remain
  // the same but the coeffs should be reversed. Wait! No! We need to modify the powers from p
  // to deg-p. Then the term array will be sorted in reverse order so we should reverse it

  void copySettingsFrom(const rsSparseFilter<TSig, TPar>& other)
  {
    num = other.num;
    den = other.den;
    updateDelayLineLength();

    // ToDo: Use num.copyFrom(other.num), den.copyFrom(other.num). maybe have a boolean parameter
    // copyState - if true, also copy the contents of the delayline from the other object. Or maybe
    // split it into two functions: copyCoeffsFrom(), copyStateFrom()
  }



  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isFilterValid() const
  {
    bool ok = true;

    // Numerator and denominator polynomials should not be empty:
    ok &= num.getNumTerms() > 0 && den.getNumTerms() > 0;

    // Filter should satisfy the a0 == 1 normalization property:
    ok &= den.getPower(0) == 0  && den.getCoeff(0) == TPar(1);

    // Length of delayline should match the maximum of the degrees of numerator and denominator:
    int maxDelay = rsMax(num.getDegree(), den.getDegree());
    ok &= delayLine.getDelayInSamples() == maxDelay;

    // Do more checks: like, the minimum power being >= zero in num and den, powers don't appear
    // twice, powers are ordered (but maybe they don't have to be - not sure yet) etc.

    // ok &= num.getMinPowerIndex() == 0

    return ok;
  }


  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    // Compute numerator N(z):
    rsComplex<TPar> N(0);
    for(int i = 0; i < num.getNumTerms(); i++)
      N += num.getCoeff(i) * rsPow(z, rsComplex<TPar>(-num.getPower(i)));

    // Compute denominator D(z):
    rsComplex<TPar> D(0);
    for(int i = 0; i < den.getNumTerms(); i++)  
      D += den.getCoeff(i) * rsPow(z, rsComplex<TPar>(-den.getPower(i)));

    // Compute transfer function H(z) = N(z) / D(z):
    return N / D;
  }



  /** Computes one output sample at a time using a direct form 2 implementation. */
  TSig getSample(TSig in)
  {
    rsAssert(isFilterValid());

    // Apply denominator as feedback part:
    TSig tmp = in;
    for(int i = 1; i < den.getNumTerms(); i++)
      tmp -= den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply numerator as feedforward path:
    tmp = 0;
    for(int i = 0; i < num.getNumTerms(); i++)
      tmp += num.getCoeff(i) * delayLine.readOutputAt(num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }

  /** Computes a sample at a time of the inverse filter. We apply the desired transformation to the
  filter into its inverse on the fly. */
  TSig getSampleInverse(TSig in)
  {
    rsAssert(isFilterValid());
    rsAssert(num.getPower(0) == 0);
    rsAssert(num.getCoeff(0) != 0);

    // Apply scaled numerator as feedback part:
    TPar s = TPar(1) / num.getCoeff(0);
    TSig tmp = in;
    for(int i = 1; i < num.getNumTerms(); i++)
      tmp -= s * num.getCoeff(i) * delayLine.readOutputAt(num.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply scaled denominator as feedforward path:
    tmp = 0;
    for(int i = 0; i < den.getNumTerms(); i++)
      tmp += s * den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }


  /** Computes a sample at a time of a filter that has the numerator transformed from min-phase to
  max-phase or vice versa. For mixed phase filters, it inverts the mix. We reflect the zeros in
  the unit circle. */
  TSig getSamplePhased(TSig in)
  {
    rsAssert(isFilterValid());

    // Apply denominator as feedback part:
    TSig tmp = in;
    for(int i = 1; i < den.getNumTerms(); i++)
      tmp -= den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply reversed numerator as feedforward path:
    int deg = num.getDegree();
    tmp = 0;
    for(int i = 0; i < num.getNumTerms(); i++)
      tmp += num.getCoeff(i) * delayLine.readOutputAt(deg - num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }
  // Needs test! This is perhaps not great for realtime use because the num.getDegree() call must
  // iterate through the whole numerator. Maybe that value coudl be cached. Not sure.







  /** Resets the filter's state. This clears the delayline. */
  void reset()
  {
    delayLine.reset();
  }



protected:

  rsSparsePolynomial<TPar> num, den;    // Numerator and denominator of transfer function
  rsBasicDelayLine<TSig> delayLine;     // Delayline used for the direct form 2 implementation

};

