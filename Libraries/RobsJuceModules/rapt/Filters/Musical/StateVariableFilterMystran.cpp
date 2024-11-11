
template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::convertToBiquad(
  TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2)
{
  getBiquadNumeratorCoeffs(b0, b1, b2);
  getBiquadDenominatorCoeffs(  a1, a2);
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffs(
  TPar* b0, TPar* b1, TPar* b2)
{
  TPar t0, t1, t2;  // Temporaries

  getBiquadNumeratorCoeffsLP(&t0, &t1, &t2);
  *b0 = aL*t0;
  *b1 = aL*t1;
  *b2 = aL*t2;

  getBiquadNumeratorCoeffsBP(&t0, &t1, &t2);
  *b0 += aB*t0;
  *b1 += aB*t1;
  *b2 += aB*t2;

  getBiquadNumeratorCoeffsHP(&t0, &t1, &t2);
  *b0 += aH*t0;
  *b1 += aH*t1;
  *b2 += aH*t2;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsLP(
  TPar* b0, TPar* b1, TPar* b2)
{
  TPar s = scl;
  TPar c = gpr;
  *b0 =   s*g*g;
  *b1 = 2*s*g*g;
  *b2 =   s*g*g;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsBP(
  TPar* b0, TPar* b1, TPar* b2)
{
  TPar s = scl;
  TPar c = gpr;
  *b0 =  g*s;
  *b1 =  0;
  *b2 = -g*s;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsHP(
  TPar* b0, TPar* b1, TPar* b2)
{
  TPar s = scl;
  TPar c = gpr;
  *b0 =  s;
  *b1 = -2*s;
  *b2 =  s;
}

template<class TSig, class TPar>
rsComplex<TPar> rsStateVariableFilterMystran<TSig, TPar>::getTransferFunctionAt(
  const rsComplex<TPar>& z)
{
  TPar b0, b1, b2, a1, a2;
  convertToBiquad(&b0, &b1, &b2, &a1, &a2);
  rsComplex<TPar> d = TPar(1)/z, d2 = d*d;                              // d = z^-1, d2 = z^-2
  rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (TPar(1) + a1*d + a2*d2);
  return H;
}

//=================================================================================================
/*

ToDo:

- Maybe have a "Muted" mode before "Bypass"

- Add a getTransferFunctionAt(std::complex<TPar> z).

- Add a setupFromBiquad(TPar b0, ...) function

- Try to achieve more general responses

- Provide a getOutputs(in, outLP, outBP, outHP) function such that the user can obtain all 3 
  outputs and mix them by themselves. 

- Figure out how to morph between LP/BP/HP, LP/AP/HP, LS/PK/HS, ...

- Maybe add inquiry functions such as getIntegratorGain() = g, getOmega() = 2*atan(g), 
  getQualityFactor() = 1 / (gpr - g). But the Q formula is wrong for bell filters and the omega
  formula is wrong for shelf filters. But maybe we can infer in which mode we are and then dispatch
  to the appropriate formula. The mode could be figured out by looking at the pattern of the mixing
  the mixing coeffs. I think, we have  LP: 1,0,0  HP: 0,0,1  BPS: 0,1,0  BPP: 0,+,0  BS: 1,0,1  
  AP: 1,-,1  PK: 1,+,1  LS: +,+,1  HS: 1,+,+.

- In the prototype folder in MiscFilters.h, there is some subclass  rsStateVariableFilterMystran2
  that extends this class by some add-on functionality. Maybe someday, some of it should be
  dragged over.

---------------------------------------------------------------------------------------------------
Algorithm

As mystran explains, the analog prototype response of this SVF is:

          a0 + a1 s + a2 s^2
  H(s) = --------------------
          1  + s/Q  + s^2

so the a-coefficients are the polynomial coefficients of the numerator of the s-domain 
transfer function. If we can manage to bring a given s-domain transfer function into this form,
then we can directly read off our mixing coeffs from the transfer function. From the RBJ
cookbook filters descirbed here:

  https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/OtherAuthors/Audio-EQ-Cookbook.txt

the lowpass, highpass, bandpass, bandstop and allpass transfer functions are indeed of this form, 
so we can directly read off our a-coeffs from these. 


For the peak/bell filter, the RBJ prototype response is of the form:

          1 + s*(A/Q) + s^2
  H(s) = -------------------
          1 + s/(A*Q) + s^2

which not exactly of the right form because in the denominator, we see an  s/(A*Q)  term 
instead of the desired  s/Q  term. But by letting  P = A*Q  we can replace  A*Q  by  P in
the denominator and in the numerator replace  Q  by  P/A  to get:
 
          1 + s*A^2/P + s^2
  H(s) = --------------------
          1 +   s/P   + s^2

This substitution can be automated using the following Sage code:

  var("s Q A P")
  H = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)
  G = H.subs(Q == P/A)
  G

which produces the output:  (A^2*s/P + s^2 + 1)/(s^2 + s/P + 1))  where G is in the desired 
form but with P instead of Q. We can now just use P in place of Q and get our mixing coeffs as
a0 = a2 = 1, a1 = A^2/P = A/Q. 


For the low shelving filter, the RBJ prototype transfer function is:

              s^2  + (sqrt(A)/Q)*s + A      A^2 + (sqrt(A)/Q)*s + s^2
  H(s) = A * --------------------------- = -----------------------------
              A*s^2 + (sqrt(A)/Q)*s + 1      1  + (sqrt(A)/Q)*s + A*s^2

Now we have two problems: the factor for s as well the one for s^2 is wrong. Instead of
1/Q and 1 as coeffs for s and s^2, we see sqrt(A)/Q and A. Both problems can be solved by
substituting  t = s*sqrt(A), i.e. s = t/sqrt(A). The following Sage snippet solves does this:

  var("s Q A t")
  H = A*(s^2+sqrt(A)/Q*s+A)/(A*s^2+sqrt(A)/Q*s+1)
  G = H.subs(s == t/sqrt(A))
  G

which produces:  (A + t^2/A + t/Q)*A/(t^2 + t/Q + 1). Again, we can read off the a0,a1,a2 
coeffs from G as a0 = A, a1 = 1/Q, a2 = 1/A. Our  s <-> t  substitution means that we now must
scale the frequencies because that's the effect of multiplying s by a factor. But those coeffs
will give a response that is off from the desired one by a scaling factor. The high frequency 
gain is supposed to be unity and it is given by coeff in front of t^2, so it would be 1/A. To get
it back to unity, we need to scale all coeffs by A such that:  a0 = A^2, a1 = A/Q, a2 = 1. 
ToDo: Figure out what went wrong to require this additional scaling!


For the high shelf the prototype transfer function is:

              A*s^2 + (sqrt(A)/Q)*s + 1     1 + (sqrt(A)/Q)*s + A*s^2
  H(s) = A * --------------------------- = ---------------------------
              s^2 + (sqrt(A)/Q)*s + A       A + (sqrt(A)/Q)*s + s^2

With this Sage code:

  var("s Q A t")
  H = A * (A*s^2 + (sqrt(A)/Q)*s + 1)/(s^2 + (sqrt(A)/Q)*s + A)
  G = H.subs(s == t*sqrt(A))
  G

We get:  (A^2*t^2 + A*t/Q + 1)*A/(A*t^2 + A + A*t/Q). Apparently, Sage didn't fully simplify the
expression. We can cancel the A to get: (A^2*t^2 + A*t/Q + 1)/(t^2 + 1 + t/Q) so we read off:
a0 = 1, a1 = A/Q, a2 = A^2. This time, we don't need to scale anything. The a0 coeff, i.e. the 
lowpass gain, already came out as 1 as it should for high-shelving filter.

*/
