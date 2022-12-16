template<class T>
void invertBiquad(T &b0, T &b1, T &b2, T &a1, T &a2)
{
  // normalize numerator:
  T g = b0;    // superfluous -> get rid
  T s = 1/g;
  b0  = 1;
  b1 *= s;
  b2 *= s;

  // swap numerator against denominator:
  rsSwap(a1, b1);
  rsSwap(a2, b2);

  // re-apply (inverted) overall gain:
  b0 *= s;
  b1 *= s;
  b2 *= s;
}

template<class T>
void twoOnePolesToBiquad(T b0[2], T b1[2], T a1[2],
  T &B0, T &B1, T &B2, T &A1, T &A2)
{
  B0 = b0[0] * b0[1];
  B1 = b0[0] * b1[1] + b0[1] * b1[0];
  B2 = b1[0] * b1[1];
  A1 = a1[1] + a1[0];
  A2 = a1[0] * a1[1];
}

template<class T>
void ensureOnePoleRootInsideUnitCircle(T &c0, T &c1)
{
  T root = -c1/c0;
  if(fabs(root) > 1.0)
  {
    // root is outside the unit circle
    root = 1.0/root;  // reflect into unit circle
    c1   = -root*c0;  // recompute coefficient for delayed sample
  }
}

// given a pair of intemediate coeffs for a one-pole filter C0, C1 (which can either be A0, A1 or
// B0, B1), this function computes the corresponding actual filter coefficients c0, c1 (either 
// a0, a1 or b0, b1) - used only internally by magnitudeMatchedOnePoleCoeffs:
template<class T>
void intermediatesToOnePoleCoeffs(T C0, T C1, T &c0, T &c1)
{
  // use p-q formula and backsubstitution:
  T x, x1, x2, tmp;
  tmp = rsSqrt(C0*C0/4 - C1*C1);
  x1  = C0/2 + tmp;
  x2  = C0/2 - tmp;
  x   = rsMax(x1, x2);
  c0  = rsSqrt(x);
  c1  = C1/c0;
}

template<class T>
void magnitudeMatchedOnePoleCoeffs(T &b0, T &b1, T &a1, T w[3], T m[3])
{
  int i;

  // powers (magnitudes squared) and cosine factors:
  T p[3], u[3];
  for(i = 0; i < 3; i++)
  {
    p[i] = m[i]*m[i];
    u[i] = 2*cos(w[i]);
  }

  // establish the matrix A;
  T A[3][3];
  for(i = 0; i < 3; i++)
  {
    A[i][0] = 1.0;
    A[i][1] = u[i];
    A[i][2] = -p[i]*u[i];
  }

  // compute the intermediate coefficients:
  T c[3];
  rsLinearAlgebra::rsSolveLinearSystem3x3(A, c, p);
  T B0 = c[0];
  T B1 = c[1];
  T A1 = c[2];

  // check, if it is possible to satisfy the constraints with a 1st order filter - if not, return
  // identity filter coefficients:
  if(B0*B0/4 - B1*B1 < 0.0 ||  1.0/4.0 - A1*A1 < 0.0)
  {
    b0 = 1.0;
    b1 = 0.0;
    a1 = 0.0;
    return;
  }

  // compute actual filter coefficients from the intermediate values: 
  T a0;
  intermediatesToOnePoleCoeffs(B0,  B1, b0, b1);
  intermediatesToOnePoleCoeffs(1.0, A1, a0, a1);

  // normalize all coefficients such that a0 == 0:
  T s = 1.0/a0;
  b0 *= s;
  b1 *= s;
  a1 *= s;
  a0  = 1.0;

  // evaluate power at w = PI/2 (needed later for renormalization): we choose PI/2, because the 
  // magnitude is easy to evaluate there and a 1st order filter cannot have zero gain there:
  T pOld = (b0*b0+b1*b1)/(1+a1*a1);

  // disambiguate filters with alike magnitude responses but possibly poles/zeros reflected in 
  // the unit circle
  ensureOnePoleRootInsideUnitCircle(a0, a1);  // selects stable filter
  ensureOnePoleRootInsideUnitCircle(b0, b1);  // selects minimum phase filter

  // this potential reflection of poles and zeros may have changed the overall gain of the 
  // filter, so we re-evaluate the power at PI/2 and compare it to the old power at PI/2 (before 
  // any reflections) and scale the feedforward coefficients appropriately:
  T pNew = (b0*b0+b1*b1)/(1+a1*a1);
  s = rsSqrt(pOld/pNew);
  b0 *= s;
  b1 *= s;

  // todo: couldn't we do the gain adjustment inside the ensureOnePoleRootInsideUnitCircle for 
  // poles and zeros separately?

  // disambiguate sign of the output by requiring b0 to be positive:
  if(b0 < 0.0)
  {
    b0 = -b0;
    b1 = -b1;
  }
}

/*

ToDo:

-Implement a function similar to 
   magnitudeMatchedOnePoleCoeffs(T &b0, T &b1, T &a1, T w[3], T m[3]); with API
   magnitudeMatchedBiquadCoeffs(T &b0, T &b1, T b2, T &a1, T &a2, T w[5], T m[5]);
 maybe list the input parameters w, m first to be more consistent with other parts of the library.

-See: http://www.rs-met.com/documents/dsp/DigitalBiquadDesignByMagnitudeRequirements.pdf
 -Try to solve the 5x5 system symbolically using Sage or Mathematica. Maybe a symbolic solution 
  could be implemented more efficiently than using a general numerical solver. The matrix actually
  has only 20 variables because the 1st column is all ones. 
 -In the 5x5 matrix, the v_n can be expressed using the double angle formula for cosines. Dunno, if
  that helps to simplify the system further. -> Try it

-Maybe the 5-point magnitude match biquads should match at the 5 following frequencies:
   0, w0/8,   w0/4,   w0/2, w0       or
   0, w0/4,   w0/2, 3*w0/4, w0       or
   0, w0/2, 3*w0/4,   w0,   5*w0/4
 in the 1st case, the required cosine factors can be computed by repeated application of the 
 double-angle formula. In the 2nd and 3rd case, we have an arithmetic progression and can use the
 cos(a+b) = ... addtion theorem. The 3rd needs to compute cos(w0/4) first and it can be used only 
 iff 5*w0/4 <= fs/2. If this is not the case, maybe switch to using the 2nd formula. The idea is to
 optimize the design to save on sin/cos calls.

-Perhaps it's possible to do a 5-point magnitude match even when some of the frequencies are above 
 fs/2? But we perhaps have to take care that we don't encounter a pair like fs/2 - df, fs/2 + df 
 which would both alias into the same freq at fs/2 - df and thereby make the matrix singular?
 -> Test that!




See:
https://www.kvraudio.com/forum/viewtopic.php?t=569113 EQ filters
https://www.kvraudio.com/forum/viewtopic.php?t=466312 Biquad coefficients using Magnitude (or Phase) Invariance Mapping Method


*/