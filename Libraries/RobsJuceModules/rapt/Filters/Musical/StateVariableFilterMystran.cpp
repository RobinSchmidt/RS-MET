
template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::setupFromBiquad(
  TPar b0, TPar b1, TPar b2, TPar a1, TPar a2)
{
  // Compute intermediates:

  TPar T1 = a1 - a2 - 1;
  TPar T2 = a1 + a2 + 1;
  TPar T  = T1 * T2;
  //TPar T = (a1*a1 - a2*a2 - 2*a2 - 1);                     // == (a1 - a2 - 1) * (a1 + a2 + 1)
  rsAssert(T < 0, "The formulas work only for T < 0.");
  TPar S = sqrt(-1/T);
  TPar r = (2*(a2 - 1) / (T*S));

  // Compute final coefficients:
  aH = -(b0 - b1 + b2) / (a1 - a2 - 1);
  aB =  2*S*(b0 - b2);
  aL =  (b0 + b1 + b2) / (a1 + a2 + 1);
  g  = -(1/((a1 - a2 - 1)*S));
  c  =  g + r;
  s  =  1 / (1 + g*c);

  // ToDo:
  //
  // - Figure out what the condition T >= 0 means. Can we deal with it somehow? For T -> 0, we have
  //   r -> inf, I think. 
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::convertToBiquad(
  TPar* b0, TPar* b1, TPar* b2, TPar* a1, TPar* a2) const
{
  getBiquadNumeratorCoeffs(b0, b1, b2);
  getBiquadDenominatorCoeffs(  a1, a2);
}

template<class TSig, class TPar>
TPar rsStateVariableFilterMystran<TSig, TPar>::getMagnitudeAt(TPar w) const
{
  rsComplex<TPar> j(0, 1);                                 // Imaginary unit
  rsComplex<TPar> z = rsExp(j*w);                          // Evaluation point in z-plane
  rsComplex<TPar> H = getTransferFunctionAt(z);            // Complex frequency response at w
  return rsAbs(H);                                         // Absolute value of H is magnitude
}

template<class TSig, class TPar>
rsComplex<TPar> rsStateVariableFilterMystran<TSig, TPar>::getTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  TPar b0, b1, b2, a1, a2;
  convertToBiquad(&b0, &b1, &b2, &a1, &a2);
  rsComplex<TPar> d = TPar(1)/z, d2 = d*d;                 // d = z^-1, d2 = z^-2
  rsComplex<TPar> H = (b0 + b1*d + b2*d2) / (TPar(1) + a1*d + a2*d2);
  return H;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadDenominatorCoeffs(
  TPar* a1, TPar* a2) const
{
  *a1 =  2*(c*g + g*g)*s - 2;
  *a2 = -2*(c*g - g*g)*s + 1;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffs(
  TPar* b0, TPar* b1, TPar* b2) const
{
  TPar t0, t1, t2;  // Temporaries
  getBiquadNumeratorCoeffsLP(&t0, &t1, &t2); *b0  = aL*t0; *b1  = aL*t1; *b2  = aL*t2;
  getBiquadNumeratorCoeffsBP(&t0, &t1, &t2); *b0 += aB*t0; *b1 += aB*t1; *b2 += aB*t2;
  getBiquadNumeratorCoeffsHP(&t0, &t1, &t2); *b0 += aH*t0; *b1 += aH*t1; *b2 += aH*t2;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsLP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =   s*g*g;
  *b1 = 2*s*g*g;
  *b2 =   s*g*g;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsBP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =  g*s;
  *b1 =  0;
  *b2 = -g*s;
}

template<class TSig, class TPar>
void rsStateVariableFilterMystran<TSig, TPar>::getBiquadNumeratorCoeffsHP(
  TPar* b0, TPar* b1, TPar* b2) const
{
  *b0 =  s;
  *b1 = -2*s;
  *b2 =  s;
}

//=================================================================================================
/*

ToDo:

- Add a setupFromBiquad(TPar b0, ...) function. See older implementation and Wishnick paper in the
  references. When done, make unit tests that test roundtrips for various settings. Maybe use also 
  random biquad coeffs in these tests (maybe with some stability constraints).

- Try to achieve more general responses. Maybe also add responses of 1st order LP, HP, LS, HS, AP
  types. I think, these make only use of the first filter/integrator stage. Having these modes 
  available in a class for a 2nd order filter can be convenient when a multimode filter should 
  also provide these 1st order modes but one doesn't want to dispatch to a different filter object 
  for getting them.

- Maybe implement a getPartialMagnitudesAt(complex z, TPar* low, TPar* band, TPar* high)
  function that produces the 3 magnitude responses for lowpass, bandpass and highpass

- Figure out how to morph between LP/BP/HP, LP/AP/HP, LS/PK/HS, ...

- Maybe add inquiry functions such as getIntegratorGain() = g, getOmega() = 2*atan(g), 
  getQualityFactor() = 1 / (gpr - g). But the Q formula is wrong for bell filters and the omega
  formula is wrong for shelf filters. But maybe we can infer in which mode we are and then dispatch
  to the appropriate formula. The mode could be figured out by looking at the pattern of the mixing
  the mixing coeffs. I think, we have  LP: 1,0,0  HP: 0,0,1  BPS: 0,1,0  BPP: 0,+,0  BS: 1,0,1  
  AP: 1,-,1  PK: 1,+,1  LS: +,+,1  HS: 1,+,+.

- In the prototype folder in MiscFilters.h, there is some subclass  rsStateVariableFilterMystran2
  that extends this class by some add-on functionality. Maybe someday, some of it should be
  dragged over. It has also a preliminary setupFromBiquad() function. But it doesn't always work.
  I think, it works only when T = (a1*a1 - a2*a2 - 2*a2 - 1) < 0.

- Figure out if there is a more direct way to evaluate the transfer function, i.e. one that 
  doesn't go through a conversion to a direct form biquad.

- Add an experiment that looks at the DC-response when switching the cutoff freq. The Wishnick 
  paper says that this is a good test for modulation response.

- Maybe rename the old implementation to rsStateVariableFilterOld and this class to 
  rsStateVariableFilter. The old code is kinda rubbish and should be deprecated.

- Maybe move the desription of the algorithm below into a separate text file.

---------------------------------------------------------------------------------------------------
Algorithm for computing the mixing coefficients aL, aB, aH

As mystran explains, the analog prototype response of this SVF is:

          aL + aB s + aH s^2
  H(s) = --------------------
          1  + s/Q  + s^2

so the a-coefficients are the polynomial coefficients of the numerator of the s-domain 
transfer function. If we can manage to bring a given s-domain transfer function into this form,
then we can directly read off our mixing coeffs from the transfer function. Among the RBJ
cookbook filters, the lowpass, highpass, bandpass, bandstop and allpass transfer functions are 
indeed of this form, so we can directly read off our a-coeffs from these. For the peak/bell 
filter, the RBJ prototype response is of the form:

          1 + s*(A/Q) + s^2
  H(s) = -------------------
          1 + s/(A*Q) + s^2

which is not exactly of the right form because in the denominator, we see an  s/(A*Q)  term 
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

---------------------------------------------------------------------------------------------------
References:

- https://www.kvraudio.com/forum/viewtopic.php?p=8992653#p8992653  
  Teemu Voipio (aka mystran) explains how to set up the mixing coefficients to achieve the well 
  known RBJ cookbook transfer functions. This discussion was what prompted me to re-implement the 
  SVF using these ideas to compute the mixing coefficients. 

- The Art of Virtual Analog Filter Design
  Vadim Zavalishin's excellent book has a chapter about ZDF-SVF filters (and much more good stuff).
  This was what my earlier implementation was based on.

- https://www.cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf
  Andrew Simper describes a filter implementation that is very similar. But this filter mixes its 
  final output not from highpass, bandpass and lowpass parts but rather from input, bandpass and 
  lowpass. These are two variations of the same filter.

- https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/FilterTransferFunctions.txt
  My (Robin Schmidt's) derivations for the formulas to convert from our coeffs here to direct form 
  biquad coeffs.

- https://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf
  Aaron Wishnick's paper has formulas (equation 16 a-c) for converting from direct form biquad 
  coefficients to SVF coeffs that can be used here. (This is not yet implemented but may be added 
  later.)

- https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/OtherAuthors/Audio-EQ-Cookbook.txt
  Robert Bristow Johnson's classic biquad filter cookbook. The response types realized here are
  precisely those described there.

*/
