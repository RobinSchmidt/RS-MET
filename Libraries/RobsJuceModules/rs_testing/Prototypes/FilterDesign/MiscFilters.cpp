// Helper function for conversion from biquad coefficients to state variable filter coeffcients:
template<class T>
void biquadToStateVariableFilter(
  T  b0, T  b1, T  b2, T  a1, T  a2,        // Inputs
  T*  g, T* R2, T* cL, T* cB, T* cH)        // Outputs
{
    // Compute intermediate values. The square roots could be imaginary but when we form their 
  // quotient and product, it will become real:
  using Complex = std::complex<T>;  
  T       u1 = -T(1) - a1 - a2;             // could be negative
  T       u2 = -T(1) + a1 - a2;             // ...dito
  // rsAssert(u1*u2 >= 0);                  // triggers in one of the unit tests
  Complex s1 = sqrt(Complex(u1));           // could be imaginary
  Complex s2 = sqrt(Complex(u2));           // ...dito  
  T       p  = real(s1 * s2);               // but their product should be real
  T       s  = T(1) / p;                    // we actually need the product's reciprocal

  // Compute coeffs:
  *g  = real(s1 / s2);                      // 16a, the quotient should also be real
  *R2 = s * T(2) * (a2 - T(1));             // 16b
  *cH = (b0 - b1 + b2) / (T(1) - a1 + a2);  // 16c, == -(b0-b1+b2) / s1 before taking the sqrt?
  *cB = s * T(2) * (b2 - b0);               // 16d, but with a factor of -1 (why?)
  *cL = (b0 + b1 + b2) / (T(1) + a1 + a2);  // 16e
  //h  = 1 / (1 + R2*g + g*g);               // factor for feedback precomputation

  // The formulas are taken from (Eq 16 a-e) here:
  // http://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf
  //
  // ToDo:
  // -Figure out why we need the factor -1 for the cB coeff with respect to the formula 16d in the 
  //  paper. Different conventions?
  // -Try to avoid complex numbers: I think, if one of the values under the sqrt gets negative, the
  //  other one must be negative, too and at the end of the day, this just results in a sign-flip
  //  in some intermediate variable. Maybe keep the original formulas in a comment for reference.
  //  But maybe one being positive and the other negative could occur for unstable filters and 
  //  maybe we want to be able to match them, too? Sometimes, they are useful. It's rare but it 
  //  happens so we'd better be prepared for it. ...OK yes - in one of the unit tests, we actually
  //  have such a case. Adding an rsAssert(u1*u2 >= 0) would trigger in this unit test. ...hmm - 
  //  but in that case p = 0  ->  s = inf, so that's not really an unstable filter but some even
  //  more drastic error condition - so maybe we can indeed assume that either u1,u2 are both 
  //  positive or both negative? Figure that out!

}
// ToDo:
// This may eventually go into RAPT::rsFilterCoefficientConverter. Maybe use rsComplex instead of
// std::complex to allow for a greater variety of types for T.
//
// Maybe make a version that computes R2+g instead of R2 and additionaly compute 
// h = 1 / (1 + R2*g + g*g); That's how the values used in the processing so that's an 
// optimization for the per sample calculations.
//
// See also rsStateVariableFilter<TSig, TPar>::setupFromBiquad. It duplicates the code and comments 
// here. It should call it from some central place instead. rsFilterCoefficientConverter seems most 
// suitable for that purpose indeed.



template<class TSig, class TCoef> 
void rsStateVariableFilterChain<TSig, TCoef>::setupFrom(
  const RAPT::rsBiquadCascade<TSig, TCoef>& bqc)
{ 
  int numStages = bqc.getNumStages();
  data.resize(numStages);
  TCoef b0, b1, b2, a1, a2;
  TCoef g,  R2, cL, cB, cH;
  for(int i = 0; i < numStages; i++)
  {
    bqc.getCoeffs(i, &b0, &b1, &b2, &a1, &a2);
    biquadToStateVariableFilter(b0, b1, b2, a1, a2, &g, &R2, &cL, &cB, &cH);
    data[i].coeffs.g    = g;
    data[i].coeffs.R2pg = R2 + g;
    data[i].coeffs.h    = 1 / (1 + R2*g + g*g); // == 1 / (1 + (R2+g)*g)
    data[i].coeffs.cL   = cL;
    data[i].coeffs.cB   = cB;
    data[i].coeffs.cH   = cH;
  }
}