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







//=================================================================================================

// Under construction

/** Implements a sparse rational function. ...TBC... */

template<class T>
class rsSparseRationalFunction
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  /** Numerator and denominator polynomials. These data members are public because it's really more
  convenient that way. We could do some sort of facade pattern and delegation but it would 
  literally just be boilerplate - here and in client code - and a lot of it. We would have to 
  implement functions like:

    void setNumeratorCoeff(int index, T newCoeff) { num.setCoeff(index, newCoeff); }

  and then client code would call things like:

    r.setNumeratorCoeff(...);

  instead of:

    r.num.setCoeff(...);

  We would have to do this basically for all setters and getters of rsSparsePolynomial - twice. 
  And there are a lot of setters and getters. Nope. Just nope! Let's make num and den public 
  instead. Yes, I'm fully aware that it goes against OOP encapsulation practices. I know the rules 
  and break them deliberately here. We do not really have to maintain any class invariants or 
  anything like that so it's ok to let client code directly access and manipulate the numerator and 
  denominator. The only donwside may be that the variable names "num" and "den" now become part of
  the public API of the class and can't be changed later. I can live with that. */
  rsSparsePolynomial<T> num, den;
  // ...well...wait: There actually is a class invariant that (maybe) should be maintained: The 
  // denominator should be nonzero...hmmm...well...or maybe we just take the position that the onus 
  // is on the client to avoid divisions by zero. That's actually also how it works for int and 
  // float. Such variables (and also rsFraction) also do not nanny the programmer that way. So why
  // should we? Or maybe I'm just too lazy to write the boilerplate and trying to rationalize it? 
  // But it's not just about writing the boilerplate. It's also about readability and bloat - not 
  // on the binary code side (the delegations would be inlined) but on the source code side. 
  //
  // Hmm...but maybe the assumption that we really want expose all the setters and getters for the
  // two polynomials is wrong? Maybe we actually want to deal with a higher level interface here?
  // If really access to the full functionality of rsSparsPolynomial is needed, we could provide
  // getters liek getNumerator/DenominatorReference() for that.
  //
  // We'll see.....



  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */


  rsSparseRationalFunction() 
  {
    den.appendTerm(T(1), 0);
  }


  rsSparseRationalFunction(
    const rsSparsePolynomial<T>& numerator, const rsSparsePolynomial<T>& denominator) 
    : num(numerator), den(denominator)  {}

  // Maybe implement it for const rsSparsePolynomial<T>&&, too. Or maybe that one is enough?


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void clear()
  {
    num.clear();
    den.clear();
    den.appendTerm(T(1), 0);
  }

  void setNumTerms(int newNumNumeratorTerms, int newNumDenominatorTerms)
  {
    num.setNumTerms(newNumNumeratorTerms);
    den.setNumTerms(newNumDenominatorTerms);
  }

  void setupFromDenseCoeffs(
    const std::vector<T>& newNumeratorCoeffs,
    const std::vector<T>& newDenominatorCoeffs,
    T tol)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs,   tol);
    den.setupFromDenseCoeffs(newDenominatorCoeffs, tol);
  }

  void setupFromDenseCoeffs(
    const T* newNumeratorCoeffs,   int newNumNumeratorTerms, 
    const T* newDenominatorCoeffs, int newNumDenominatorTerms,
    T tol)
  {
    num.setupFromDenseCoeffs(newNumeratorCoeffs,   newNumNumeratorTerms,   tol);
    den.setupFromDenseCoeffs(newDenominatorCoeffs, newNumDenominatorTerms, tol);
  }


  void copyDataFrom(const rsSparseRationalFunction<T>& q)
  {
    num.copyDataFrom(q.num);
    den.copyDataFrom(q.den);
  }



  /** Applies a scaling factor to this rational function. This basically means to scale all 
  numerator coeffs by that factor. */
  void scale(T scaler) { num.scale(scaler); }


  void multiplyBy(rsMonomial<T> factor) { num.multiplyBy(factor); }


  void multiplyBy(const rsSparseRationalFunction<T>& factor, T tol) 
  { 
    num.multiplyBy(factor.num, tol);
    den.multiplyBy(factor.den, tol);

    // I think, num and den may now have a common factor, so we potentially need to divide that
    // out:
    //canonicalize();
    // But maybe we should not automatically reduce the result by default because doing so may 
    // require memory allocations (because the GCD algo needs temporaries) and we really need this
    // function to be realtime safe (it's used in rsDampedCombAllpass::getCombTransferFunction(), 
    // for example - and that is used in rsDampedMultiCombAllpass::updateFilters() which could 
    // potentially be called on an audio thread). Maybe we should have a boolean parameter 
    // "reduce"? that deafults to true but that we set to false in a realtime context?)

    // Consider 14/15 * 3/4 = 42/60 = 7/10. Although both factors are in lowest terms, their 
    // product is not. The same thing could happen with rational functions.
  }
  



  void addConstant(T constant, T tol)
  {
    num.addScaledPolynomial(den, constant, tol);
  }


  //void canonicalize();
  // Should: (1) Divide out the GCD of num and den. (2) Canonicalize num and den. 
  // (3) Divide num and den by the leading coeff of den (i.e. make den monic)
  //
  // But maybe it could also make sense to define a canonical representation as one with monic
  // numerator? But no! This can't represent the zero function.


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  bool isCloseTo(const rsSparseRationalFunction<T>& q, T tol) const
  {
    return q.num.isCloseTo(num, tol) && q.den.isCloseTo(den, tol);
  }


  // isCanonical()
  // A canonical representation has canonical numerator and denominator with no common factors
  // and the denominator is monic


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Evaluates the function at the given input x. */
  T operator()(T x) const { return num(x) / den(x); }

  /** Evaluates the function at the given input z whose type may be different from the 
  coefficient type, for example, for evaluating functions with real coeffs at complex arguments.
  WARNING: the same considerations as for @see rsPolynomial::operator(TArg) apply. */
  template<class TArg>
  TArg operator()(TArg z) const { return num(z) / den(z); }

  /** Adds two rational functions. */
  rsSparseRationalFunction<T> operator+(const rsSparseRationalFunction<T>& q) const 
  { rsSparseRationalFunction<T> r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  /** Subtracts two rational functions. */
  rsSparseRationalFunction<T> operator-(const rsSparseRationalFunction<T>& q) const 
  { rsSparseRationalFunction<T> r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  /** Multiplies two rational functions. */
  rsSparseRationalFunction<T> operator*(const rsSparseRationalFunction<T>& q) const 
  { return rsSparseRationalFunction(num * q.num, den * q.den); }

  /** Divides two rational functions. */
  rsSparseRationalFunction<T> operator/(const rsSparseRationalFunction<T>& q) const 
  { return rsSparseRationalFunction(num * q.den, den * q.num); }


  //-----------------------------------------------------------------------------------------------
  /** \name Low level API.  */


  static void weightedSum(
    const rsSparseRationalFunction<T>& p, T wp,
    const rsSparseRationalFunction<T>& q, T wq,
    rsSparseRationalFunction<T>* r, T tol);



  static void weightedSumDestructive(
    rsSparseRationalFunction<T>* p, T wp,
    rsSparseRationalFunction<T>* q, T wq,
    rsSparseRationalFunction<T>* r, T tol);
  // The first parameter p may alias to the result r. 





};

/** Multiplies a coefficient and a sparse rational function. */
template<class T>
inline rsSparseRationalFunction<T> operator*(const T& s, const rsSparseRationalFunction<T>& p)
{
  rsSparseRationalFunction<T> r(p);
  r.scale(s);
  return r;
}

template<class T>
void rsSparseRationalFunction<T>::weightedSum(
  const rsSparseRationalFunction<T>& p, T wp,
  const rsSparseRationalFunction<T>& q, T wq,
  rsSparseRationalFunction<T>* r, T tol)
{
  r->den = p.den * q.den;
  r->num = wp * p.num * q.den  +  wq * q.num * p.den;

  // This can probably be optimized with respect to avoid unnecessary temporary objects and heap
  // allocations. We may also use the gcd instead of just cross-mutiplying the denominators.
}

template<class T>
void rsSparseRationalFunction<T>::weightedSumDestructive(
  rsSparseRationalFunction<T>* p, T wp,
  rsSparseRationalFunction<T>* q, T wq,
  rsSparseRationalFunction<T>* r, T tol)

{
  rsAssert(rsAreAddressesDistinct(*p, *q));
  //rsAssert(rsAreAddressesDistinct(*r, *p));  // We may actually allow this!
  rsAssert(rsAreAddressesDistinct(*r, *q));
  // Maybe we can relax this? It would be really nice if r could be equal to at least one of p or
  // q. Requiring p and q to be distinct is not such a big problem. I think, we can allow this, if
  // SP::weightedSum can work in place. But at the moment. I think, it can't. But maybe it can be
  // made so. Looking at the code, it seems like it could work when the result aliases to the 1st 
  // argument. Test and document this! A test indicates that this may indeed work out. Investigate
  // this further and document! We could perhaps make it work to also allow r == q by swapping p 
  // and q (and wp and wq) in this case. But what if r == p == q? ...well...in that case, we could 
  // leave the denominator of r (and p and q) alone and just multiply the numerator by wp+wq, I 
  // think. I think, the p == q case could possibly also be handled

  using SP = rsSparsePolynomial<T>;

  //              arg1        arg2        result
  SP::multiply(   p->num,     q->den,     &p->num, tol);  // Replace p->num by p->num * q->den
  SP::multiply(   q->num,     p->den,     &q->num, tol);  // Replace q->num by q->num * p->den
  SP::weightedSum(p->num, wp, q->num, wq, &r->num, tol);  // Establish r->num
  SP::multiply(   p->den,     q->den,     &r->den, tol);  // Establish r->den

  // ToDo:
  //
  // - Document exactly, how it can be used with respect to which pointers must be distinct.
}
// Needs tests



//=================================================================================================


/** A subclass of rsSparseRationalFunction that is meant to deal specifically with transfer 
functions of digital filters. A general transfer function for a digital filter looks like:

          b0 + b1 * z^-1 + b2 * z^-2 + ... + bN * z^-N
  H(z) = ----------------------------------------------
          1  + a1 * z^-1 + a2 * z^-2 + ... + aM * z^-M

Note that this is actually a rational function not in z itself but in z^-1. We override the 
function evaluation operator () to take care of this reciprocation of z. We also implement some 
additional functionality on top of the baseclass that is specific to such transfer functions. For 
example, digital filter transfer functions are usually normalized to a0 = 1, as seen above. We 
implement a check for that condition (and a few others) isCanonical(). We also provide functions to
invert the transfer function (basically, swapping numerator and denominator but maintaining that 
the a0 = 1 still holds after the swap), reflecting the zeros about the unit circle (turning 
minimum phase filters into maximum phase ones), etc. ...TBC... */


template<class T>
class rsSparseDigitalTransferFunction : public rsSparseRationalFunction<T>
{

public:

  using Base = rsSparseRationalFunction<T>;    // For convenience
  using Base::Base;                            // Inherit constructors



  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

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


    // Maybe do num.shiftPowers(rsMax(amountInSamples, -getPreDelay()) ); and write unit tests for
    // this
  }

  /** Removes the predelay from this filter, if any is present. This makes sure that the lowest
  exponent of z^-1 in the numerator is zero. see getPreDelay(), addPreDelay()  */
  void removePreDelay() { num.shiftPowers(-getPreDelay()); }

  /** Turns the filter into its inverse. This basically amounts to swapping numerator and
  denominator and possibly applying some scaling of the coefficients if b0 != 1. */
  void invert()
  {
    rsAssert(isCanonical());

    // A filter with predelay cannot be inverted in realtime. The best thing we can do in this 
    // case is to invert the filter up to the predelay. Removing the predelay ensures that the
    // 0-th term in the numerator has power of 0, i.e. it's a  b0 * z^-0  term and not some crazy
    // b7 * z^-7  term:
    removePreDelay();
    rsAssert(num.getPower(0) == 0); // Numerator is already asserted to be non-empty in isCanoncial
    rsAssert(num.getCoeff(0) != 0); // so we can access the 0-th element without risk here

    // Swap numerator and denominator while maintaining the a0 = 0 normalization condition:
    T s = T(1) / num.getCoeff(0);
    scale(s);
    std::swap(num, den);
    scale(s);

    // I'm pretty sure it doesnt' allocate. Verify and document.
  }

  /** Reflects the zeros of the filter about the unit circle. This will turn a minimum phase
  filter into a maximum phase one and vice versa. For mixed phase filters, it inverts the mix.
  It doesn't affect stability or filter order. */
  void reflectZeros()
  {
    int deg = num.getDegree();
    for(int i = 0; i < num.getNumTerms(); i++)
      num.setPower(i, deg - num.getPower(i));
    num.reverse();                               // Order array by ascending powers again

    // How about a reflectPoles() function? But that would turn stable filters into unstable ones,
    // so it's usefulness is questionable. For the time being, we can do without.
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the predelay introduced by this filter. A predelay is characterized by the fact that
  the lowest exponent of z^-1 in the numerator is not zero. That means the numerator does not look
  like b0 + b1*z-^1 + b2*z^-2 + ... but rather something like b7*z^-7 + b8*z^-8 + b9*z^-9 + ...
  for a predelay of 7 samples, for example. */
  int getPreDelay() const { return num.getPower(0); }

  /** Returns the order of the filter. This is the maximum exponent of z^-1 that occurs in the
  transfer function. */
  int getFilterOrder() const { return rsMax(num.getDegree(), den.getDegree()); }

  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isCanonical() const
  {
    bool ok = true;

    // Numerator and denominator polynomials should not be empty:
    ok &= num.getNumTerms() > 0 && den.getNumTerms() > 0;

    // We assume the filter polynomials to be in canonical shape:
    ok &= num.isCanonical();
    ok &= den.isCanonical();

    // Filter should satisfy the a0 == 1 normalization property:
    ok &= den.getPower(0) == 0 && den.getCoeff(0) == T(1);

    return ok;
  }


  /** Computes the density of the numerator defined as the number of actual nonzero coeffs divided
  by the number of potentially nonzero coeffs given the degree of the numerator. */
  double getNumeratorDensity() const
  {
    return double(num.getNumTerms()) / double(num.getDegree()+1);
  }

  /** Computes the density of the denominator defined as the number of actual nonzero coeffs 
  divided by the number of potentially nonzero coeffs given the degree of the denominator. */
  double getDenominatorDensity() const
  {
    return double(den.getNumTerms()-1) / double(den.getDegree());
  }

  /** Returns the "combined density" defined as the number of actual nonzero coeffs of the filter 
  divided by the number of potential nonzero coeffs for the given filter order. The a0 coeff 
  doesn't count because it's always 1. */
  double getCombinedDensity() const
  {
    int numPossibleCoeffs = 2*getFilterOrder() + 1;
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }

  /** Returns the "separated density" defined as the as number of actual nonzero coeffs in 
  numerator and denominator divided by the number of potential nozero coeffs for the given orders
  of numerator and denominator. */
  double getSeparatedDensity() const
  {
    int numPossibleCoeffs = num.getDegree()+1 + den.getDegree();
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }
  // ToDo: Explain this better! We consider numerator and denominator as separate filters that can
  // have their own orders and therefore the computation of the number of possibly nonzero coeffs
  // is different. For example, in an Nth order allpole filter, we have a 0th order numerator. In 
  // getCombinedDensity, we would assume that it could potentially have N+1 coeffs. Here, we assume
  // that it can have only one because the order of the numerator is zero.

  // Maybe this could be moved into the baseclass. But I'm not sure, if the +1 fo the num and
  // no +1 for the den also applies there...well...I think, it does when we assume a canonical
  // representation with monic denomionator. Here we normalize the denominator to a0=1 - but
  // whatever way we normalize the function, the denominator has always one degree of freedom
  // less than the numerator (for the same degree). A degree 1 polynomial has 2 coeffs and a degree
  // N polynomial has N+1 coeffs. That number applies to the numerator as is. But the denominator
  // is normalized so we lose one degree of freedom and subtract 1 again.

  // ToDo: Maybe return the densities as rsFraction<int>. They are rational numbers so maybe we 
  // should treat them as such.





  T operator()(T x) const 
  { 
    T xr = T(1) / x;
    return num(xr) / den(xr); 
  }


  template<class TArg>
  TArg operator()(TArg z) const 
  { 
    TArg zr = TArg(1) / z;
    return num(zr) / den(zr);
  }
  // Reciprocation of z needed because we store the coeffs of H(z^-1)


  // This boilerplate is needed to have the desired arithmetic operators available also for the 
  // derived class. They can not be inherited from the baseclass because their parameter and
  // return types are different. It may work with pointer-types but not with value-types (I guess):


  rsSparseDigitalTransferFunction<T> operator+(const rsSparseDigitalTransferFunction<T>& q) const 
  { rsSparseDigitalTransferFunction<T> r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }


  rsSparseDigitalTransferFunction<T> operator-(const rsSparseDigitalTransferFunction<T>& q) const 
  { rsSparseDigitalTransferFunction<T> r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }


  rsSparseDigitalTransferFunction<T> operator*(const rsSparseDigitalTransferFunction<T>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.num, den * q.den); }

  rsSparseDigitalTransferFunction<T> operator/(const rsSparseDigitalTransferFunction<T>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.den, den * q.num); }





};

/** Multiplies a coefficient and a sparse digital transfer function. */
template<class T>
inline rsSparseDigitalTransferFunction<T> operator*(
  const T& s, const rsSparseDigitalTransferFunction<T>& p)
{
  rsSparseDigitalTransferFunction<T> r(p);
  r.scale(s);
  return r;
}




//=================================================================================================

/** A class for representing sparse filters in direct form. We represent them using an object of
type rsSparseRationalFunction to store the coefficients of the transfer function H(z). Well, we
actually store the coeffs of H(z-^1) there because that's what's needed for implementation of the
difference equation. The filter is implemented in direct form 2 using a single delayline. */

template<class TSig, class TPar>
class rsSparseFilter
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the filter from dense arrays of numerator and denominator coeffs. When a coefficient
  in the dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<TPar>& numCoeffs, 
    const std::vector<TPar>& denCoeffs, TPar tol)
  {
    H.num.setupFromDenseCoeffs(numCoeffs, tol);
    H.den.setupFromDenseCoeffs(denCoeffs, tol);
    updateDelayLineLength();
  }
  // This may allocate!
  // ToDo: use num.setupFromDenseCoeffs(&numCoeffs[0], (int) numCoeffs.size(), tol)
  // use  H.setupFromDenseCoeffs(numCoeffs, denCoeffs, tol)


  void setup(const rsSparseRationalFunction<TPar>& newTransferFunction)
  { H.copyDataFrom(newTransferFunction); updateDelayLineLength(); }
  // The parameter should really be of type rsSparseDigitalTransferFunction


  void setNumNumeratorTerms(int newNumTerms) { H.num.setNumTerms(newNumTerms); }
  // This may allocate!

  void setNumDenominatorTerms(int newNumTerms) { H.den.setNumTerms(newNumTerms); }
  // This may allocate!

  void setNumeratorTerm(int index, TPar coeff, int delay) { H.num.setTerm(index, coeff, delay); }

  void setDenominatorTerm(int index, TPar coeff, int delay) { H.den.setTerm(index, coeff, delay); }
  // Actually, we really should call updateDelayLineLength() after setting a term because it 
  // potentially requires a change of the length. But: updateDelayLineLength() is expensive and 
  // setting terms is an operation that might be called in a loop or sequence in which case only
  // one update after the sequence of calls should be done. Maybe when calling it in a sequence,
  // we should use other functions like setNumeratorTermSuppressDelayUpdate. Or maybe give the 
  // function a boolean parameter updateDelayLength which defaults tor true

  // Maybe call them setNumeratorTermNoUpdate(). ...or maybe get rid of these function entirely and
  // instead give the user functions like
  // rsSparsePolynomial<TPar>& getNumeratorRef() for low level access. That breaks encapsulation
  // though. Maybe move them into some extra section "Low level API"



  /** Ensures that the delayline has enough memory allocated to support the desired transfer
  function. This may re-allocate heap memory in cases where the delayline does not already have
  enough capacity, so you don't want to call it on a realtime thread. */
  void ensureEnoughDelayMemory()
  { setMaxDelayInSamples(rsMax(getMaxDelayInSamples(), getFilterOrder())); }

  void setMaxDelayInSamples(int newMaxDelay)
  { delayLine.setMaxDelayInSamples(newMaxDelay); }


  /** Updates the length of the delayline according to the maximum power of z^-1 that occurs in
  numerator and denominator polynomial. */
  void updateDelayLineLength()
  { 
    rsAssert(hasEnoughDelayMemory(), "Not enough delay memory!");
    // When this happens, it means that you are trying to request a transfer function from the
    // filter that it can't support because you didn't pre-allocate enough delay memory. At some 
    // point, you need to call setMaxDelayInSamples() where you set up the maximum possible delay
    // and therefore the maximum possible filter order. If you don't allocate enough and request a
    // too high filter order later, this assertion will trigger. We do not just re-allocate here 
    // because this function is designed to be realtime safe. In a realtime context, you really
    // need to pre-allocate enough on construction or in some prepareToPlay() function or something
    // like that.

    delayLine.setDelayInSamples(getFilterOrder()); 
  }
  // Maybe we actually should re-allocate if necessary but still leave the assertion in. Maybe 
  // that's the most benign way to recover from the error condition in a release build? But nah!
  //
  // ToDo: Figure out and document what happens, when we ignore this error. I think, the delay
  // time will be wrapped around / bitmasked by the actual maxDelay in rsDelay member.
  // I say "actual" because the vaue you set up via setMaxDelayInSamples() may be "rounded up"
  // to the next power of two minus one...or something. See implementation of 
  // rsDelay::readOutputAt(). So if the actual maxDelay is 15, all delays will be 
  // interpreted modulo 16, so readOutputAt(20) would actually amount to a delay of 
  // 20 % 16 = 4 rather than 20. ...verify this!




  /** Applies a scaling factor to the filter. This basically means to scale all numerator coeffs by
  that factor. */
  void scale(TPar scaler) { H.num.scale(scaler); }

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given 
  amount. */
  void addPreDelay(int delayInSamples) { H.addPreDelay(delayInSamples); updateDelayLineLength(); }

  /** Removes any predelay that may be present in the filter. */
  void removePreDelay() { H.removePreDelay(); updateDelayLineLength(); }

  /** Turns the filter into its inverse. */
  void invert() { H.invert(); updateDelayLineLength(); }

  /** Reflects the zeros of the filter about the unit circle. */
  void reflectZeros() { H.reflectZeros(); }
  // Doesn't change the required delayline length so we don't need to call updateDelayLineLength


  /** Copies the settings (i.e. the coefficients) from the given other filter into this object and
  possibly adjusts the delayline length, if necessary. */
  void copySettingsFrom(const rsSparseFilter<TSig, TPar>& other)
  { H.copyDataFrom(other.H); updateDelayLineLength(); }

  // Maybe also provide a copyStateFrom() method which would also copy the content of the 
  // delayline. We could then also have a function copyDataFrom that calls both. Maybe the 
  // delayline class should also have a function copyStateFrom (or maybe copyDataFrom)



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the order of the filter. This is the maximum amount of delay needed to implement 
  the filter. */
  int getFilterOrder() const 
  { 
    return H.getFilterOrder();
    //return rsMax(H.num.getDegree(), H.den.getDegree()); 
  }

  int getMaxDelayInSamples() const { return delayLine.getMaxDelayInSamples(); }
  // Maybe rename to getMaxFilterOrder


  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isFilterValid() const { return H.isCanonical() && areDelaysConsistent(); }
  // ToDo: Elaborate documentation. Give some details about what it checks.





  /** Computes the transfer function H(z) of this filter at the given complex value z. */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const { return H(z); }
  // Yes - that's right! "return H(z)" is the whole implementation. Isn't that elegant? :-D


  // Reciprocation of z needed because H actually stores the coeffs of H(z^-1)
  // ToDo: factor the reciprocation out into the () operator of
  // rsSparseDigitalTransferFunction ...done!


  /** Returns a const reference to our transfer function object H(z). */
  const rsSparseDigitalTransferFunction<TPar>& getTransferFunction() const { return H; }




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  /** Computes one output sample at a time using a direct form 2 implementation. */
  TSig getSample(TSig in)
  {
    rsAssert(isFilterValid());

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < H.den.getNumTerms(); i++)
      tmp -= H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply numerator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < H.num.getNumTerms(); i++)
      tmp += H.num.getCoeff(i) * delayLine.readOutputAt(H.num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }

  /** Computes a sample at a time of the inverse filter. We apply the desired transformation to the
  filter into its inverse on the fly. */
  TSig getSampleInverse(TSig in)
  {
    rsAssert(isFilterValid());
    rsAssert(H.num.getPower(0) == 0);
    rsAssert(H.num.getCoeff(0) != 0);
    // Maybe we can relax this? If we do not expect the power of the 0-th coeff to be 0, we will 
    // just produce an inverted filter up to delay?

    // Apply scaled numerator of H as feedback part:
    TPar s = TPar(1) / H.num.getCoeff(0);
    TSig tmp = in;
    for(int i = 1; i < H.num.getNumTerms(); i++)
      tmp -= s * H.num.getCoeff(i) * delayLine.readOutputAt(H.num.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply scaled denominator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < H.den.getNumTerms(); i++)
      tmp += s * H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));

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

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < H.den.getNumTerms(); i++)
      tmp -= H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply reversed numerator of H as feedforward path:
    int deg = H.num.getDegree();
    tmp = 0;
    for(int i = 0; i < H.num.getNumTerms(); i++)
      tmp += H.num.getCoeff(i) * delayLine.readOutputAt(deg - H.num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }
  // Needs test! This is perhaps not great for realtime use because the H.num.getDegree() call must
  // iterate through the whole numerator. Maybe that value could be cached. Not sure. Although,
  // If we assume H.num to be in canonical representation (which it is, I think), then getDegree()
  // can be replaced by  H.getCoeff(H.getNumTerms()-1)  which avoids the iteration. Maybe such a 
  // call could even be encapsulated into something like H.getLastPower(). Maybe add an
  // H.isCanonical() check to isFilterValid().
  //
  // Maybe factor out the getSample() functions into free functions like 
  //
  // TSig rsGetSample(TSig in,
  //        const rsSparseDigitalTransferFunction<TPar>& H, rsDelay<TSig>* delay);
  //
  // This facilitates memory optimizations in situations where we have multiple filters with the
  // same set of coeffs but independent states, i.e. independent delaylines. We can store the
  // coeffs of H once and use them with different delaylines. As is currently is, we would have to
  // create a sparse filter object for each of the filters and therfore store the coeffs 
  // redundantly. This might become a general pattern for implementing filters: separate coeffs and
  // state and provide free functions that take a const ref to the coeffs and (mutable) a pointer
  // to the state.







  /** Resets the filter's state. This clears the delayline. */
  void reset() { delayLine.reset(); }



protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Self test */

  /** Checks if the maximum delay required by the transfer function (i.e. the highest power of 
  z^-1 that occurrs) is consistent with the length of the delayline. This is meant for internal 
  sanity checks. */
  bool areDelaysConsistent() const { return delayLine.getDelayInSamples() == getFilterOrder(); }
  // Maybe rename to areDelayAndOrderConsistent, doesDelayMatchOrder

  /** Checks, if the delayline has enough memory allocated to support the transfer function H(z).
  The maximum possible delay must be greater or equal to the order of the filter. */
  bool hasEnoughDelayMemory() const 
  { return delayLine.getMaxDelayInSamples() >= getFilterOrder(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  rsDelay<TSig> delayLine;         // Delayline for the direct form 2 implementation.
  rsSparseDigitalTransferFunction<TPar> H;  // Transfer function H(z). Contains filter coeffs.

};

