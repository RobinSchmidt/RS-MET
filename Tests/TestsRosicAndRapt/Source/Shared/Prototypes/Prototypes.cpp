#include "Prototypes.h"

#include "ClassicFilterDesign/PoleZeroPrototype.cpp"
#include "ClassicFilterDesign/PoleZeroMapper.cpp"
#include "ClassicFilterDesign/PoleZeroDesignerAnalog.cpp"
#include "ClassicFilterDesign/PoleZeroDesignerDigital.cpp"

#include "Projection3Dto2D.cpp"

using namespace RAPT;

double signalValueViaSincAt(double *x, int N, double t, double sincLength, double stretch,
  //FunctionPointer3DoublesToDouble windowFunction, 
  double (*windowFunction)(double,double,double),
  double windowParameter, bool normalizeDC)
{
  int L  = (int) floor(sincLength);
  int ti = (int) floor(t);         // integer part of t
  if( ti < 0 || ti >= N )
    return 0.0;
  double tf = t - ti;              // fractional part of t    
  double y  = 0.0;                 // output value
  double b;                        // weight for a sample
  double bs = 0.0;                 // sum of the weights
  int mMin  = -rsMin(L/2-1, ti);   // minimum shift
  int mMax  = +rsMin(L/2, N-ti-1); // maximum shift
  int    nr;                       // sample read position
  double ts;                       // time instant with shift
  for(int m = mMin; m <= mMax; m++)
  {
    ts = m-tf;
    b  = rsNormalizedSinc(ts/stretch) * windowFunction(ts, sincLength, windowParameter);
    bs += b;
    nr = ti+m;
    y += b * x[nr]; 
  }
  if( normalizeDC == true )
    return y/bs;
  else
    return y/stretch;
}

void halpernU(double *a, int K)
{
  RAPT::rsArray::fillWithZeros(a, K+1);
  int k, m;
  int f1, f2, f3, f4; // factorials
  if( rsIsEven(K) )
  {
    k = K/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*m] = (-1)^(k-m) * (m+k)! / ( (k-m)! * (m!)^2 ):
      f1 = rsFactorial(m+k);
      f2 = rsFactorial(k-m);
      f3 = rsFactorial(m);
      a[2*m] = pow(-1.0, k-m) * f1 / (f2*f3*f3);
    }
  }
  else
  {
    k = (K-1)/2;
    for(m = 0; m <= k; m++)
    {
      // a[2*(k-m)+1] = (-1)^m * (2*k+1-m)! / (m! * (k+1-m)! * (k-m)!):
      f1 = rsFactorial(2*k+1-m);
      f2 = rsFactorial(m);
      f3 = rsFactorial(k+1-m);
      f4 = rsFactorial(k-m);
      a[2*(k-m)+1] = pow(-1.0, m) * f1 / (f2*f3*f4);
    }
  }

  // we leave the normalization factors out, in order to keep all numbers integers. the squared
  // factors (which are also integers) will be applied, after the coefficient array was convolved
  // with itself
}
void halpernT2(double *c, int N)
{
  double u[30];  
  double s[30];
  double t[30];

  // create the U-polynomial:
  halpernU(u, N-1);

  // square it:
  RAPT::rsArray::convolve(u, N, u, N, s);

  // apply squared scale factors:
  int k;
  if( rsIsEven(N) )
  {
    k = (N-2)/2;
    RAPT::rsArray::scale(s, 2*N-1, 4*(k+1));
  }
  else
  {
    k = (N-1)/2;
    RAPT::rsArray::scale(s, 2*N-1, 4*k+2);
  }

  // multiply it by x:
  RAPT::rsArray::rightShift(s, 2*N, 1);

  // integrate it from 0 to w:
  double a[1] = { 0 };    // lower integration limit (as polynomial)
  double b[2] = { 0, 1 }; // upper integration limit (as polynomial)
  rsPolynomial<double>::integratePolynomialWithPolynomialLimits(s, 2*N-1, a, 0, b, 1, t);
    // this is actually the same as just taking the antiderivative of s

  // copy to output:
  RAPT::rsArray::copyBuffer(t, c, 2*N+1);
}

void papoulisL2(double *v, int N)
{
  double P[20];  // array for Legendre polynomials
  int k;
  if( rsIsOdd(N) )
  {
    k = (N-1)/2;
    RAPT::rsArray::fillWithZeros(v, k+1);
    for(int r = 0; r <= k; r++)
    {
      // add the weighted Legendre polynomial of order r to our v polynomial:
      rsPolynomial<double>::legendrePolynomial(P, r);
      rsPolynomial<double>::weightedSumOfPolynomials(v, r, 1.0, P, r, 2*r+1.0, v);
    }

    // square it:
    RAPT::rsArray::convolve(v, k+1, v, k+1, v);

    // account for leaving out the division by (k+1)*sqrt(2) in the weights in the accumulation:
    RAPT::rsArray::scale(v, 2*k+1, 1.0/(2*(k+1)*(k+1)));

    // integrate from -1 to 2*w^2-1
    double a[1] = { -1 };
    double b[3] = { -1, 0, 2};
    rsPolynomial<double>::integratePolynomialWithPolynomialLimits(v, 2*k, a, 0, b, 2, v);
  }
  else
  {
    k = (N-2)/2;

    // generate Legendre Polynomial of order k+1, store in P:
    rsPolynomial<double>::legendrePolynomial(P, k+1);

    // take the derivative, store in v:
    rsPolynomial<double>::polyDerivative(P, v, k+1);  // v has order k, length k+1

    // square it:
    RAPT::rsArray::convolve(v, k+1, v, k+1, v); // has order 2*k, length 2*k+1

    // multiply by (x+1):
    double xp1[2] = { 1, 1};
    RAPT::rsArray::convolve(v, 2*k+1, xp1, 2, v); // has order 2*k+1, length 2*k+2=N

    // integrate from -1 to 2*w^2-1
    double a[1] = { -1 };
    double b[3] = { -1, 0, 2};
    rsPolynomial<double>::integratePolynomialWithPolynomialLimits(v, 2*k+1, a, 0, b, 2, v); // order: 2*N, length 2*N+1

    // scale, such that L^2(1) = 1:
    double K = 1.0 / rsPolynomial<double>::evaluatePolynomialAt(1.0, v, 2*N);
    RAPT::rsArray::scale(v, 2*N+1, K); 
  }
}

//=================================================================================================

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setPoles(CRPar p1re, CRPar p1im, CRPar p2re, CRPar p2im)
{
  xx = p1re;
  yx = p1im;
  yy = p2re;
  xy = p2im;
  makePolesDistinct();  // avoid singular case of equal poles by fudging

  // The state update matrix will have one of these two general forms:
  // |p1 0 |     or:     r * |cos(w)  -sin(w)| 
  // |0  p2|                 |sin(w)   cos(w)|
  // where in the first case, p1 and p2 are the two real poles and the x and y states decay 
  // exponentially and independently from each other when the input is switched off. In the second 
  // case, the numbers r*cos(w), r*sin(w) are the real and imaginary parts of a pair of complex 
  // conjugate poles and we will see a spiraling/decaying rotation of the states when there's no 
  // input (anymore).
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setImpulseResponseStart(TPar h[3])
{
  TPar A[2][2];
  TPar c[2];
  A[0][0] = xx + xy;
  A[0][1] = yx + yy;
  A[1][0] = xx*A[0][0] + xy*A[0][1];
  A[1][1] = yx*A[0][0] + yy*A[0][1];
  rsLinearAlgebra::rsSolveLinearSystem2x2(A, c, &h[1]);
  cx = c[0];
  cy = c[1];
  ci = h[0] - cx - cy;

  // The first 3 impulse response samples of this filter can be calculated as:
  // h[0] = cx + cy + ci
  // h[1] = cx*(xx + xy) + cy*(yx + yy)
  // h[2] = cx*(xx*(xx+xy) + xy*(yx+yy)) + cy*(yx*(xx+xy) + yy*(yx+yy))
  // and this must match what is given to this function. From this, we can set up a system of
  // 3 linear equations for the 3 mixing coefficients and solve it. It can actually be solved as a 
  // 2x2 system and the 3rd equation is then trivial.

  // Problem: The A-matrix may become singular. This happens when we have two equal real poles as
  // happens in a cookbook lowpass with q = 0.5.
  // ...then, the filter can't be expressed as a parallel connection of two first order
  // filters - what to do in this case? ...
  // maybe work with a fudged determinant, like D_f = sign(D) * max(abs(D), fudgeFactor) where
  // fudgefactor is some small nonzero number
  // rewrite the equation system as:
  // |a b| * |x1| = |y1|
  // |c d|   |x2|   |y2|
  // and use D = a*d - b*c as preliminary determinant and then call D = fudgeDet(D) or something
  // ...very hacky/dirty/ugly but it may work (the filter may be slightly misadjusted in these
  // problematic cases but otherwise work well). don't call rsSolveLinearSystem2x2, write out the
  // equations here

  // or maybe use a fudgePoles function called at the end of setPoles to make sure, they are
  // somewhat distinct
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::setupFromBiquad(
  CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2)
{
  // compute and set up poles from a1, a2:
  TPar d = TPar(0.25)*a1*a1 - a2;  // discriminant, term inside sqrt
  TPar p = -a1*0.5;                // -p/2 in p-q-formula, term before +/- sqrt
  if(d >= 0) {                     // d >= 0: we have 2 real poles
    TPar sq = sqrt(d);
    setPoles(p+sq, 0, p-sq, 0);
  }
  else {                           // d < 0: we have complex conjugate poles
    TPar sq = sqrt(-d);
    setPoles(p, sq, p, -sq);
  }

  // compute first 3 samples of biquad impulse response and set up mixing coeffs:
  TPar h[3]; 
  h[0] = b0;
  h[1] = b1 - a1*h[0];
  h[2] = b2 - a1*h[1] - a2*h[0];
  setImpulseResponseStart(h);
}

template<class TSig, class TPar>
void rsStateVectorFilter<TSig, TPar>::makePolesDistinct()
{
  TPar minDelta = 0.001; // actually  2 * minimum pole-distance
    // with 0.001, impulse response of a cookbook lowpass (1kHz@44.1kHz, q=0.5) is visually 
    // indistiguishable from the correct, desired one and the mixing coeffs are around 9.
    // With 0.0001, mixing coeffs raise to about 90. Maybe make more formal tests, 
    // how  the upper bound of the mixing coeffs behaves as function of minDelta and choose 
    // something reasonable. The c-coeffs should ideally stay of the same order of magnitude as in
    // the naturally nonsingular cases (or at least not be vastly higher)

  TPar delta;
  if(xy == 0) { // two real poles
    delta = xx - yy;
    if(abs(delta) < minDelta) {
      TPar avg = 0.5*(xx+yy);
      if(xx >= yy) {
        xx = avg + minDelta;
        yy = avg - minDelta;
      }
      else {
        xx = avg - minDelta;
        yy = avg + minDelta;
      }
    }
  }
  else {
    delta = xy + yx; // + because they have opposite signs
    if(abs(delta) < minDelta) {
      if(xy > 0) {
        xy = +minDelta;
        yx = -minDelta;
      }
      else {
        xy = -minDelta;
        yx = +minDelta;
      }
    }
  }

  // todo: optimize: avoid the makePolesDistinct when it can be predicted that they are
  // distinct anyway (after computing sqrt(d) or sqrt(-d)) ...but not in prototype
}


template class rsStateVectorFilter<double, double>; // explicit instantiation

// can the mixing coeffs be computed more simply - from the b-coeffs of the biquad?

// todo: make functions setupFromPolesAndZeros - i think, the update matrix is (re,-im;im,re)
// setPoles(p1re, p1im, p2re, p2im)
// setImpRespStart(h[3])
// maybe refactor setupFromBiquad into pole (update matrix) and zero (mixing coeffs) computation
// in order to bypass the complicated pole computations when possible
// maybe also factor out a setImpulseResponse function which takes the desied first 3 samples
// of the impulse response - this may be useful when we already know the poles but want to match
// a filter that is not necessarily a biquad (for example, a state-variable filter)

// setupFromParameters - for various parametrizations
// like setupFromFrqDecAmpPhs (for the "modal" filter - take care to take into account the 
// injection into both parts - shifts phase by 45° and increases amplitude by sqrt(2))
// setFreqAndReso(freq, reso) 
// freq controls pole angle, reso is pole radius - maybe allow values > 1 and use saturation
// circle-saturation/distortion or (soft)clipping

// is it possible to find simpler formulas? it looks all very complicated - maybe from equating the
// transfer functions of biquad and svf? i think, the TF of this filter is:
// X(z) = 1 / (1 - xx * z^-1 * X(z) - xy * z^-1 * Y(z))
// Y(z) = 1 / (1 - yx * z^-1 * X(z) - yy * z^-1 * Y(z))
// H(z) = ci + cx*X(z) + cy*Y(z)

// todo: make unit tests that compare outputs of biquads for various settings all combinations 
// of real/complex poles/zeros, also with unused higher order coeffs (like b2=0, etc)
// against the outputs of the state vector (and also state variable) implementation
// ->also good to have in place when we later change a-coeff sign convention

// maybe try a nonlinearity: multiply x and y by 1/(1 + x^2 + y^2) after state update
// saturates(with some foldover)/contracts state vector without changing its angle
// maybe apply this factor also to "in" because it would be weird to pass the input through
// undistorted ...but might be interesting to explore
