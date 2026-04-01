
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

  // Check sanity:
  if(T >= 0)
  {
    rsError("The formulas work only for T < 0.");
    setupMuted();
    return;
    // For T > 0, S and r are NaN because we are trying to take a square root of a negative number. 
    // If this happens, it means your biquad coeffs were unstable. We produce muted output in this 
    // case. It doesn't happen for all unstable biquads, though. For some unstable biquads, the 
    // formulas still work. In these cases, we'll just use them anyway and you'll get a 
    // corresponding unstable SVF. I'm not sure, if that behavior is best, though. 
    // 
    // ToDo: Document what happens in the T == 0 case. We include it here in the if-statement such 
    // when T == 0, we also mute the output but the comment only talks about the T > 0 case. We 
    // have T = T1*T2 such that T == 0 implies T1 == 0 or T2 == 0 but below, we divide by T1,T2,
    // so it does indeed seem to be important to include the T == 0 case in the if condition to
    // avoid division by zero. Such a division by zero actually occurs above in the compuation of S
    // and r, so they should be inf (or maybe nan) in this case. It doesn't really matter, though 
    // because we then don't use these values for any further computations. ...TBC...
  }

  // Compute final coefficients:
  aH =  (b1 - B ) / T1;
  aB =  (b0 - b2) * 2*S;
  aL =  (b1 + B ) / T2;
  g  = -1 / (T1*S);
  c  =  g + r;
  s  =  1 / (1 + g*c);
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
  for getting them. Try also to realize a 2-pole resonator. Maybe that's what the "peak" mode in 
  Andy Simper's SVF is?

- Figure out how to morph between LP/BP/HP, LP/AP/HP, LS/PK/HS, ... Maybe try this by manually 
  setting the mixing coeffs for LP,BP,HP.

- Add tilt mode. See: https://github.com/zalthyrexor/QuasarEQ/blob/main/Source/zlth_dsp_filter.h

- Try to figure out how to translate the mixing coefficients between the SVF variation that mixes
  input, bandpass and lowpass with the variation that mixes highpass, bandpass and lowpass (which 
  is the variant implemented here). I think, we may need this to adapt the tilt mode above to this
  filter. Maybe try first to get the input back from HP,BP,LP.

- Maybe add inquiry functions such as getIntegratorGain() = g, getOmega() = 2*atan(g), 
  getQualityFactor() = 1 / (gpr - g). But the Q formula is wrong for bell filters and the omega
  formula is wrong for shelf filters. But maybe we can infer in which mode we are and then dispatch
  to the appropriate formula. The mode could be figured out by looking at the pattern of the mixing
  coeffs. I think, we have  LP: 1,0,0  HP: 0,0,1  BPS: 0,1,0  BPP: 0,+,0  BS: 1,0,1  AP: 1,-,1  
  PK: 1,+,1  LS: +,+,1  HS: 1,+,+. In the prototype folder in MiscFilters.h, there is some subclass
  rsStateVariableFilter2 that extends this class by some add-on functionality that includes things
  like that. Maybe someday, some of it should be dragged over. 

- Figure out if there is a more direct way to evaluate the transfer function, i.e. one that 
  doesn't go through a conversion to a direct form biquad. Somewhere is a text file where I convert
  between SVF and state-space filter coeffs. Maybe that could be useful for evaluating the transfer
  function, too? See here:
  https://github.com/zalthyrexor/QuasarEQ/blob/main/Source/zlth_dsp_filter.h
  But this filter is based on mixing input with bandpass and lowpass rather than highpass with 
  bandpass and lowpass. But for our structure here, a similar formula should be possible. Go back 
  to the derivations here:
  https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/StateVariableFilter.txt
  to figure out the formula for H(z) directly in terms of our SVF coeffs. The conversion from/to
  biquad should be kept as a nice feature anyway but we shouldn't use it to compute the transfer
  function anymore.

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

- How could we adapt this filter to non-uniformly sampled data? Could we just scale the g-coeffs in 
  the update equations? Or maybe we would need to call setup...() before each sample with 
  w = 2*pi*freq*dt where dt is the time increment at that sample? But should that be t[n]-t[n-1]
  or t[n+1]-t[n] or maybe the (weighted?) average of both? I assume that the t-array gives 
  timestamps for the corresponding samples x[n].

- Try to implement a function setupFromAnalogBiquad. The goal is to realize a filter with s-domain
  transfer function H(s) = (B0 + B1*s + B2*s^2) / (A0 + A1*s + A2*s^2) maybe normalized to A0 = 1.
  It may use BLT or MZT - or maybe we should have two separate functions for these purposes.

- Add classes for chains of state variable filters - with equal and with different coeffs per 
  stage.

- Maybe use inline or RS_INLINE also for the setup... functions. In the context of a synthesizer,
  they will typically be called at sample-rate due to envelope and LFO on the cutoff. However, in 
  other contexts (like an equalizer), the settings may be static - so I'm not sure if we really 
  want to always inline them. It may bloat the code (although: verify if the produced assembly code
  is actually bigger - the function call overhead might not be negligible in this case). It would 
  generally be really nice if we could control inlining at the call site. Figure out, if that is 
  possible with "modern" C++. If so, maybe use it.

- Add a setCoeffs() function where the user can directly set the parameters g,r,aL,aB,aH. Maybe
  in certain contexts, this may open up opportunities for optimizations. For example, when 
  designing shelving filters from a specification where the desired gain is given in dB, we could 
  compute sA = sqrt(A) = rsDbToAmp(0.5*dB) and from that we can produce A = sA*sA thereby saving
  the call to sqrt(). But with an API that expects A such as the current one, that is not possible.
  There is some tension between API convenience and efficiency.

- In setupFromBiquad(..):

  - Check what happens in the limit as T -> 0 from below. The T in the denominator of the 
    formula for r approaches 0 as well, but S approaches infinity but more slowly due to the 
    sqrt. So, I guess, overall r should approach infinity but sublinearly, namely as sqrt. 
    Verify that and give an interpretation for what that means. Perhaps rather than going into
    muted mode, we should go into bypass mode? And what if T > 0? But we also have a2 in the 
    numerator and if the approaches 1, the numerator approaches 0.
  
  - Check if we need some tolerance, i.e. if  T < 0  is not good enough but we rather need 
    something like  T < -tol  where tol is some small positive number like the machine espilon
    or some multiple of it or its sqrt.
  
  - We may be able to reduce the number of divisions by defining T1 = 1/(a1-A); T2 = 1/(a1+A).
    T = T1*T2 stays the same; S = sqrt(-T); r = 2*(a2 - 1) * (T*S); aH =  (b1 - B ) * T1;
    (b0 - b2) * 2/S; (b1 + B ) * T2; g = -1 * (T1*S); ...I think. That would be 3 divisions
    instead of 5 (not counting the one in s = ..., because that's unaffected). Hmm - I tried but
    it doesn't seem to work -> check the math! Wait: I think, the computation of S still needs
    the reciprocation - but then we would save only one division, I think.

  - The intermediate variable r is only used once, so maybe get rid of it. Maybe do the
    if(T >= 0)... test immediately after computing T. That would require to drag the S = ...
    computation into the lower part which would kind of invalidate the comments that say
    "intermediate variables" and "final coeffs", though.
  
  - Maybe wrap constants like 2,1,-1 into TPar().

  - What about those biquads that have poles exactly on the unit circle? I guess, those are the
    ones with T = 0? Can we realize them, too? Filters with poles on the unit circle can be 
    useful as sinusoidal oscillators. Maybe set up some tests with bandpasses with very high Q.
    Check what happens to the coefficients. Using setupBandpassSkirt(TPar w, TPar Q) with 
    infinite Q should lead to r = 0; c = g; s = 1/(1+g^2); That looks reasonable. Try it! Maybe
    try also lowpass and highpass with infinite Q.

*/
