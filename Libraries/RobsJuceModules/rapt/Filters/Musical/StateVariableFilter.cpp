
template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupFromBiquad(
  TPar b0, TPar b1, TPar b2, TPar a1, TPar a2)
{
  // Compute intermediates:
  TPar A  = a2 + 1;
  TPar B  = b2 + b0;
  TPar T1 = a1 - A;                                        // == a1 - a2 - 1
  TPar T2 = a1 + A;                                        // == a1 + a2 + 1
  TPar T  = T1 * T2;
  TPar S  = sqrt(-1 / T);
  TPar r  = 2*(a2 - 1) / (T*S);
  if(T >= 0)
  {
    rsError("The formulas work only for T < 0.");          // For T > 0, S and r are NaN
    setupMuted();                                          // We produce muted output in this case.
    return;
  }

  // Compute final coefficients:
  aH =  (b1 - B ) / T1;
  aB =  (b0 - b2) * 2*S;
  aL =  (b1 + B ) / T2;
  g  = -1 / (T1*S);
  c  =  g + r;
  s  =  1 / (1 + g*c);

  // ToDo:
  //
  // - Check what happens in the limit as T -> 0 from below. The T in the denominator of the 
  //   formula for r approaches 0 as well, but S approaches infinity but more slowly due to the 
  //   sqrt. So, I guess, overall r should approach infinity but sublinearly, namely as sqrt. 
  //   Verify that and give an interpretation for what that means. Perhaps rather than going into
  //   muted mode, we should go into bypass mode? And what if T > 0? But we also have a2 in the 
  //   numerator and if the approaches 1, the numerator approaches 0.
  //
  // - Check if we need some tolerance, i.e. if  T < 0  is not good enough but we rather need 
  //   something like  T < -tol  where tol is some small positive number like the machine espilon
  //   or some multiple of it or its sqrt.
  //
  // - We may be able to reduce the number of divisions by defining T1 = 1/(a1-A); T2 = 1/(a1+A).
  //   T = T1*T2 stays the same; S = sqrt(-T); r = 2*(a2 - 1) * (T*S); aH =  (b1 - B ) * T1;
  //   (b0 - b2) * 2/S; (b1 + B ) * T2; g = -1 * (T1*S); ...I think. That would be 3 divisions
  //   instead of 5 (not counting the one in s = ..., because that's unaffected). Hmm - I tried but
  //   it doesn't seem to work -> check the math! Wait: I think, the computation of S still needs
  //   the reciprocation - but then we would save only one division, I think.
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::convertToBiquad(
  TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2) const
{
  getBiquadNumeratorCoeffs(b0, b1, b2);
  getBiquadDenominatorCoeffs(  a1, a2);
}

template<class TSig, class TPar>
TPar rsStateVariableFilter<TSig, TPar>::getMagnitudeAt(TPar w) const
{
  rsComplex<TPar> j(0, 1);                                 // Imaginary unit
  rsComplex<TPar> z = rsExp(j*w);                          // Evaluation point in z-plane
  rsComplex<TPar> H = getTransferFunctionAt(z);            // Complex frequency response at w
  return rsAbs(H);                                         // Absolute value of H is magnitude
}

template<class TSig, class TPar>
rsComplex<TPar> rsStateVariableFilter<TSig, TPar>::getTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  TPar b0, b1, b2, a1, a2;
  convertToBiquad(&b0, &b1, &b2, &a1, &a2);
  rsComplex<TPar> d = TPar(1)/z, d2 = d*d;                 // d = z^-1, d2 = z^-2
  rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (TPar(1) + a1*d + a2*d2);
  return H;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::getBiquadDenominatorCoeffs(
  TPar* a1, TPar* a2) const
{
  *a1 =  2*(c*g + g*g)*s - 2;
  *a2 = -2*(c*g - g*g)*s + 1;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::getBiquadNumeratorCoeffs(
  TPar* b0, TPar* b1, TPar* b2) const
{
  TPar t0, t1, t2;  // Temporaries
  getBiquadNumeratorCoeffsLP(&t0, &t1, &t2); *b0  = aL*t0; *b1  = aL*t1; *b2  = aL*t2;
  getBiquadNumeratorCoeffsBP(&t0, &t1, &t2); *b0 += aB*t0; *b1 += aB*t1; *b2 += aB*t2;
  getBiquadNumeratorCoeffsHP(&t0, &t1, &t2); *b0 += aH*t0; *b1 += aH*t1; *b2 += aH*t2;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::getBiquadNumeratorCoeffsLP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =   s*g*g;
  *b1 = 2*s*g*g;
  *b2 =   s*g*g;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::getBiquadNumeratorCoeffsBP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =  g*s;
  *b1 =  0;
  *b2 = -g*s;
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::getBiquadNumeratorCoeffsHP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =  s;
  *b1 = -2*s;
  *b2 =  s;
}

//=================================================================================================
/*

ToDo:

- Try to achieve more general responses. Maybe also add responses of 1st order LP, HP, LS, HS, AP
  types. I think, these make only use of the first filter/integrator stage. Having these modes 
  available in a class for a 2nd order filter can be convenient when a multimode filter should 
  also provide these 1st order modes but one doesn't want to dispatch to a different filter object 
  for getting them.

- Figure out how to morph between LP/BP/HP, LP/AP/HP, LS/PK/HS, ...

- Maybe add inquiry functions such as getIntegratorGain() = g, getOmega() = 2*atan(g), 
  getQualityFactor() = 1 / (gpr - g). But the Q formula is wrong for bell filters and the omega
  formula is wrong for shelf filters. But maybe we can infer in which mode we are and then dispatch
  to the appropriate formula. The mode could be figured out by looking at the pattern of the mixing
  the mixing coeffs. I think, we have  LP: 1,0,0  HP: 0,0,1  BPS: 0,1,0  BPP: 0,+,0  BS: 1,0,1  
  AP: 1,-,1  PK: 1,+,1  LS: +,+,1  HS: 1,+,+. In the prototype folder in MiscFilters.h, there is 
  some subclass  rsStateVariableFilter2 that extends this class by some add-on functionality that
  includes things like that. Maybe someday, some of it should be dragged over. 

- Figure out if there is a more direct way to evaluate the transfer function, i.e. one that 
  doesn't go through a conversion to a direct form biquad.

- Add an experiment that looks at the DC-response when switching the cutoff freq. The Wishnick 
  paper says that this is a good test for modulation response.

- Implement functions like getMagnitudeResponse(const TPar* omegas, TPar* magnitudes, int N) that
  computes the magnitudes at a whole array of frequencies. Rationale: If one wants to compute the
  magnitude response for an array of frequencies using the existing getMagnitudeAt() function for 
  each of the frequencies, there will be a lot of redundant calculations because the SVF -> DF
  conversion will be done for each frequency anew even though the resulting coeffs will always be 
  the same. A similar function could be done for the phase response. Or maybe make a function that
  computes the complex frequency response

- Maybe try to use two independent integrator gains g1, g2 like this:

  template<class TSig, class TPar>
  inline void rsStateVariableFilter<TSig, TPar>::getPartialOutputs(
    TSig in, TSig* yL, TSig* yB, TSig* yH)  
  {
    // Compute outputs:
    *yH = (in - c*z1 - z2) * s;            // == (in - (g+r)*z1 - z2) / (1 + g*(g+r))
    *yB = z1 + g1 * *yH; 
    *yL = z2 + g2 * *yB;

    // State variable update:
    z1 += 2 * g1 * *yH;
    z2 += 2 * g2 * *yB;
  }

  and see, if this increases the space of realizable biquad transfer functions, i.e. solves the
  T >= 0 problem in setupFromBiquad().

- Maybe rename the occurences of variable r to R2 = 2*R to make the naming consistent with Vadim's 
  book.

*/
