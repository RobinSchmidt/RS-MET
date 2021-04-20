template<class TSig, class TPar>
rsLadderFilter<TSig, TPar>::rsLadderFilter()
{
  reset();
  setMode(Mode::LP_24);
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::setCutoff(CRPar newCutoff)
{
  cutoff = newCutoff;
  updateCoefficients();
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::setSampleRate(CRPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateCoefficients();
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::setResonance(CRPar newResonance)
{
  resonance = newResonance;
  updateCoefficients();
}

template<class TSig, class TPar>
void  rsLadderFilter<TSig, TPar>::setMixingCoefficients(
  CRPar c0, CRPar c1, CRPar c2, CRPar c3, CRPar c4)
{
  c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3; c[4] = c4;
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::setMode(int newMode)
{
  // Shorthands for convenience:
  using T = TPar;
  //using M = Mode;
  auto set = [&](T c0, T c1, T c2, T c3, T c4, T s)
  {
    setMixingCoefficients(c0, c1, c2, c3, c4); 
    this->s = s;
  };

  // The actual setup:
  if( newMode != mode )
  {
    mode = newMode;
    switch(mode)
    {                                                              // Prototype transfer function
    case FLAT:     { set(1,  0,   0,   0,  0, T(0.125)); } break;  // 1

      // lowpasses:
    case LP_6:     { set(0,  1,   0,   0,  0, T(0.25) ); } break;  // (1/(1+s))^1
    case LP_12:    { set(0,  0,   1,   0,  0, T(0.5)  ); } break;  // (1/(1+s))^2
    case LP_18:    { set(0,  0,   0,   1,  0, T(0.75) ); } break;  // (1/(1+s))^3
    case LP_24:    { set(0,  0,   0,   0,  1, T(1)    ); } break;  // (1/(1+s))^4
    // What about 1/(1+s^4) etc., i.e. Butterworth style? I think, that's not possible because
    // we are constrained by the rule that all poles should be equal

    // highpasses:
    case HP_6:     { set(1, -1,   0,   0,  0, T(0.25) ); } break;  // (s/(1+s))^1
    case HP_12:    { set(1, -2,   1,   0,  0, T(0.5)  ); } break;  // (s/(1+s))^2
    case HP_18:    { set(1, -3,   3,  -1,  0, T(0.75) ); } break;  // (s/(1+s))^3
    case HP_24:    { set(1, -4,   6,  -4,  1, T(1)    ); } break;  // (s/(1+s))^4

    // bandpasses:
    case BP_6_18:  { set(0,  0,   0,   4, -4, T(0.125)); } break;  // (s/(1+s))^1 * (1/(1+s))^3
    case BP_12_12: { set(0,  0,   4,  -8,  4, T(0.125)); } break;  // (s/(1+s))^2 * (1/(1+s))^2
    case BP_18_6:  { set(0,  4, -12,  12, -4, T(0.125)); } break;  // (s/(1+s))^3 * (1/(1+s))^1
    case BP_6_12:  { set(0,  0,   3,  -3,  0, T(0.125)); } break;  // (s/(1+s))^1 * (1/(1+s))^2
    case BP_12_6:  { set(0,  3,  -6,   3,  0, T(0.125)); } break;  // (s/(1+s))^2 * (1/(1+s))^1
    case BP_6_6:   { set(0,  2,  -2,   0,  0, T(0.125)); } break;  // (s/(1+s))^1 * (1/(1+s))^1

    //// these additional modes were found in a thread on KVR: 
    // http://www.kvraudio.com/forum/viewtopic.php?&t=466588 
    //case KVR_NF2:  setMixingCoefficients(1, -2,  2,  0,  0);  break;
    //case KVR_NF4:  setMixingCoefficients(1, -4,  8, -8,  4);  break;
    //case KVR_PF2:  setMixingCoefficients(1, -1,  1,  0,  0);  break;
    //case KVR_PF4:  setMixingCoefficients(1, -2,  3, -2,  1);  break;
    // the 1st two are notches, the 2nd two allpasses

    default:       setMixingCoefficients(0,  0,  0,  0,  0);  break; // out of range -> silence
    }
    updateCoefficients();
  }

  // The values for s have been found by trial and error to make the peak gain roughly equal 
  // between the various modes.
}

// inquiry:

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::getState(TSig *state)
{
  for(int i = 0; i < 5; i++) // later use: ArrayFunctions::copy(y, state, 5);
    state[i] = y[i];
}

template<class TSig, class TPar>
std::complex<TPar> rsLadderFilter<TSig, TPar>::getTransferFunctionAt(const std::complex<TPar>& z)
{
  using T = TPar;
  std::complex<T> G1, G2, G3, G4; // transfer functions of n-th stage output, n = 1..4
  std::complex<T> H;              // transfer function with resonance   
  std::complex<T> one(1, 0);

  T b; // = TPar(1) + a;
  if(bilinear)
  {
    b  = T(0.5) * (T(1) + a);   // factor out int getB(T a, bool bilinear)
    G1 = b*(one + one/z) / (one + a/z);
  }
  else
  {
    b  = T(1) + a;
    G1 = b / (one + a/z);
  }
  G2 = G1*G1;
  G3 = G2*G1;
  G4 = G3*G1;
  H  = g * (c[0] + c[1]*G1 + c[2]*G2 + c[3]*G3 + c[4]*G4) / (one + k * G4 / z);

  return H;
}

template<class TSig, class TPar>
TPar rsLadderFilter<TSig, TPar>::getMagnitudeResponseAt(CRPar frequency)
{
  TPar w = 2 * TPar(PI) * frequency/sampleRate;
  std::complex<TPar> j(0, 1);                      // imaginary unit
  std::complex<TPar> z = rsExp(j*w);               // location in the z-plane
  std::complex<TPar> H = getTransferFunctionAt(z); // H(z) at our z
  H *= conj(H);                                    // magnitude-squared
  return rsSqrt(H.real());                         // imaginary part should be zero anyway

  // Why not just return rsAbs(H)?

  // I wonder, if a simpler formula is possible which avoids going through the complex transfer 
  // function -> computer algebra
}

template<class TSig, class TPar>
rsRationalFunction<TPar> rsLadderFilter<TSig, TPar>::getTransferFunction()
{
  // Compute some intermediate variables:
  using T = TPar;
  T b  = T(1)+a;
  T b2 = b*b;        // b^2
  T b4 = b2*b2;      // b^4
  T a2 = a*a;        // a^2
  T d0 = c[0];       // c0 * b^0
  T d1 = c[1]*b;     // c1 * b^1
  T d2 = c[2]*b2;    // c2 * b^2
  T d3 = c[3]*b2*b;  // c3 * b^3
  T d4 = c[4]*b4;    // c4 * b^4

  // Compute coefficients of the numerator (note the binomial coeffs in the columns:
  // (1),(1,1),(1,2,1),(1,3,3,1),(1,4,6,4,1)):
  std::vector<T> N(5);
  N[4] = g * (d4 + d3 +      d2 +      d1 +      d0);
  N[3] = g * (     d3 + T(2)*d2 + T(3)*d1 + T(4)*d0) * a;
  N[2] = g * (               d2 + T(3)*d1 + T(6)*d0) * a2;
  N[1] = g * (                         d1 + T(4)*d0) * a2*a;
  N[0] = g * (                                   d0) * a2*a2;

  // Create rational function object and return it (note again the binomial coeffs 1,4,6,4,1) with 
  // the additional b4*k added in, so the denominator seems to be (a+1)^4 + b4*k*z^3
  using RF = RAPT::rsRationalFunction<T>;
  RF H(N, {a2*a2, T(4)*a*a2, T(6)*a2, T(4)*a + b4*k, T(1) });
  return H;
}
// ToDo:
// Maybe make a function getCoeffsBA(TPar* b, TPar* a) which just fills the arrays (of length 5).
// Beware that the order must be reversed for a polynomial in z^-1 (maybe make it optional to be
// in this form)
// Check, if the patterns with the binomial coeffs extend to an N-stage ladder. If so, write down
// the general form of the transfer function in canonical form and implement a function that 
// computes numerator and denominator coeffs for the general case.

template<class TSig, class TPar>
rsRationalFunction<TPar> rsLadderFilter<TSig, TPar>::getTransferFunctionOld()
{
  using RF = RAPT::rsRationalFunction<TPar>;
  TPar tol = 1024 * RS_EPS(TPar);  // ad hoc
  TPar b = 1+a;
  RF G1( { 0, b }, { a, 1 }, tol); // G1(z) = b / (1 + a/z) = (0 + b*z) / (a + 1*z)
  RF one({ 1    }, { 1    }, tol); // 1 = 1 / 1
  RF z  ({ 0, 1 }, { 1    }, tol); // z = (0 + 1*z) / 1
  RF G2 = G1*G1;                   // G1^2
  RF G3 = G1*G2;                   // G1^3
  RF G4 = G2*G2;                   // G1^4
  RF H  = g * (c[0]*one + c[1]*G1 + c[2]*G2 + c[3]*G3 + c[4]*G4) / (one + k * G4 / z); // H(z)
  return H;
}

// audio processing:

template<class TSig, class TPar>
inline TSig rsLadderFilter<TSig, TPar>::getSampleNoGain(CRSig in)
//inline TSig rsLadderFilter<TSig, TPar>::getSampleNoGain(TSig in)
{
  // Apply feedback and saturation:
  //y[4] /= 1 + y[4]*y[4];   // (ad hoc) nonlinearity applied to the feedback signal
  // try: y[4] /= 1 / (y[4]*y[4])^n, see https://www.desmos.com/calculator/rjlnhzqecs
  // (x+x^3)/(1+x^4)  https://www.desmos.com/calculator/wb39hmukry
  // (x+0.25*x^3)/(1+x^4)  ...very straight in the middle, goes through (0.5,0.5)
  // (x+0.25*x^3)/(1+x^2)

  y[0]  = in - k*y[4];        // linear
  y[0]  = rsClip(y[0], TSig(-1), TSig(+1));  // cheapest
  //y[0] /= TSig(1) + y[0]*y[0]; // nonlinearity applied to input plus feedback signal (division could be interesting with complex signals)
  //y[0]  = rsNormalizedSigmoids<TSig>::softClipHexic(y[0]);  // most expensive
  // ToDo: let the user select the saturationMode: hardClip, tanh, softClip, etc. maybe with an 
  // optional DC (add before the saturation, subtract after) ..do that in subclasses - here, we just do the 
  // cheap hardclipping

  // Apply the 1st order stages:
  y[1] = y[0] + a * (y[0] - y[1]);
  y[2] = y[1] + a * (y[1] - y[2]);
  y[3] = y[2] + a * (y[2] - y[3]);
  y[4] = y[3] + a * (y[3] - y[4]);

  // Form linear combination of the stage outputs:
  return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];

  // We should experiment with placing saturation at different points - and, of course, try more 
  // realistic functions. The y = x / (1 + x^2) that is currently used is not even a proper
  // sigmoid shape (it goes down to zero again) ...but maybe that's not necessarily a bad thing.
  // Maybe try y = b + (x-b) / (1 + a*x^2) with adjustable parameters a and b. "a" adjusts the 
  // amount of signal squasheing whereas "b" introduces asymmetry.

  // maybe try a different update euqation of the form:
  //   y[i] = y[i-1] + a * (y[i-1] - y[i]);
  // as is used in the rsSmoothingFilter. Elan says, this responds better to modulation. OK, i 
  // tried it with quite harsh cutoff modulation (by a square wave) and both formulas produce the 
  // same result. The difference of the produced files is zero. ...so - maybe make performance 
  // tests which one is more efficient and use that. maybe investigate the difference in higher 
  // precision - maybe one formula is numerically more precise than the other? Maybe test 
  // reso-modulation, too. In this case, it's particluarly interesting, if applying the 
  // compensation gain pre- or post-filter is the better way. I guess, pre-filter will smooth out
  // the discontinuities that these switches cause. But pre/post gain should also be checked with
  // a nonlinearity: pre-gain will drive the distortion much harder at high resonance
}

template<class TSig, class TPar>
inline TSig rsLadderFilter<TSig, TPar>::getSample(CRSig in)
{
  return getSampleNoGain(g * in);     // apply gain at input
  //return g * getSampleNoGain(in);   // apply gain at output
  
  // \todo Make the amount of gain compensation available as user parameter - we then compute
  // g = 1 + k * compensationAmount; instead of g = 1 + k -> filter becomes continuously adjustable 
  // between no compensation and full compensation. Then, we also don't really need the factored 
  // out getSampleNoGain fucntion anymore

  // it would perhaps make more sense to apply the compensation gain at the input side rather
  // than the output because it depends on the feedback gain k which scales y[4], so it may be the
  // case that gain*in - k*y[4] makes both terms fit together better. make tests with both variants
  // maybe use a unit impulse as input and switch between reso = 0 and reso = 1 (and back) in the 
  // middle of the signal and see which variant produces a smoother output. but on the other hand, 
  // the clipping behavior may be les desirable - the gain would be boosted before the clipper. But 
  // maybe that sounds better? Maybe we could also distribute the gain between pre and post via 
  // another parameter: preGain = gain^(gainDistribution), postGain = gain^(1-gainDistribution)
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::process(TSig *in, TSig *out, int length)
{
  for(int n = 0; n < length; n++)
    out[n] = getSample(in[n]);
}

// misc:

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::reset()
{
  for(int i = 0; i < 5; i++) // later use: ArrayFunctions::clear(y, 5);
    y[i] = TSig(0);
}

// coefficient computations:

template<class TSig, class TPar>
TPar rsLadderFilter<TSig, TPar>::computeFeedbackFactor(CRPar fb, CRPar cosWc, CRPar a)
{
  TPar b  = TPar(1) + a;
  TPar g2 = b*b / (1 + a*a + 2*a*cosWc);
  return fb / (g2*g2);
}
// for the bilinear case, we get g2 = (2*b+2*b^2*c) / (1+a^2+2*a*c) where c = cos(wc) = cosWc
// ...but what is b? i don't think, it's just 1+a. calculate the mag-resp for b=1 at w=0, then use
// the reciprocal of that for b

template<class TSig, class TPar>
TPar rsLadderFilter<TSig, TPar>::resonanceDecayToFeedbackGain(CRPar decay, CRPar cutoff)
{
  return rsExp(-1/(decay*cutoff)); // does this return 0 for decay == 0? -> test

  //if(decay > 0.0) // doesn't work with SIMD
  //  return exp(-1/(decay*cutoff));
  //else
  //  return 0.0;

  // The time tr for a sinusoid at the cutoff frequency fc to complete a single roundtrip around 
  // the filter loop is given by tr = 1/fc. The amplitude of the sinusoid as function of t is given
  // by a(t) = r^(t/tr) = r^(t*fc). The decaytime is defined as the time instant where a(t) = 1/e, 
  // so we need to solve 1/e = r^(t*fc) leading to r = (1/e)^(1/(t*fc)).
  // can be expressed as exp(-1/(decay*cutoff)) - avoid expensive pow
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::computeCoeffs(CRPar wc, CRPar fb, CRPar s, TPar *a, 
  TPar *k, TPar *g, bool bilinear)
{
  computeCoeffs(wc, fb, a, k, bilinear);
  *g = 1 + s * *k;

  //// old formula - has the same problem:
  //// evaluate the magnitude response at DC:
  //TPar b0_4   = *b * *b * *b * *b;
  //TPar dcGain = b0_4 / ((((*a + 4.0) * *a +6.0) * *a + 4.0) * *a + *k * b0_4 + 1.0);
  //*g = 1 / dcGain;

  // The simpler formula can be derived from the more complicated one by observing that the 
  // denominator is given by:
  //   d = a^4 + 4*a^3 + 6*a^2 + 4*a + 1 + k*(1+a)^4 = (1+a)^4 + k*(1+a)^4 = (1+k) * (1+a)^4
  // and the numerator is:
  //   b0^4 = (1+a)^4 
  // so we have
  //   dcGain = b0^4 / d = (1+a)^4 / ((1+k) * (1+a)^4) = 1 / (1+k)
  // and the reciprocal is just
  //   g = 1+k
  // update the pdf to include this and when done, this comment can be deleted
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::computeCoeffs(CRPar wc, CRPar fb, TPar *a, TPar *k, bool bilinear)
{
  TPar s, c, t;                     // sin(wc), cos(wc), tan((wc-PI)/4)
  //rsSinCos(wc, &s, &c);
  s  = rsSin(wc);
  c  = rsCos(wc);
  t  = (TPar) rsTan(0.25*(wc-PI));
  if(bilinear)
    *a = (c*t+s+t) / (s-(c+1)*t);  
     // todo: check, if feedback formula works for this, too. ...i don't think so. we need to 
     // evaluate the magnitude response of a bilinear 1st order filter at wc
  else
    *a = t / (s-c*t);
  *k = computeFeedbackFactor(fb, c, *a);
  // If the cutoff frequency goes to zero wc -> 0, then s -> 0, c -> 1, t -> -1, b -> 0, a -> -1.
  // The coefficient computation for the lowpass stages approaches this limit correctly, but the 
  // formula for the feedback factor k runs into numerical problems when wc -> 0. However, we know
  // from the analog filter that the correct limit for the feedback factor is k = 4*fb, since the
  // analog limit corresponds to infinite samplerate. Also, the formula for the compensation gain
  // becomes invalid. Maybe we should look at the magnitude response at zero frequency of the 
  // analog filter to get a gain formula from the feedback factor for this. Then, we could use a 
  // lower threshold for wc, below which these limiting case formulas are used. For the feedback 
  // computation, a lower limit of wc = 1.e-6 seems appropriate (for gain compensation, i did not 
  // yet figure it out)

  // ToDo: for optimization, we could use 1-dimensional tables for a and the k-multiplier 
  // (k = fb*multiplier) as a function of wc (in the range 0..2*pi) - that will solve also the
  // problem with formulas becoming numerically ill behaved for cutoff frequencies near zero
}

// internal functions:

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::updateCoefficients()
{
  TPar wc = 2 * (TPar)PI * cutoff / sampleRate;
  computeCoeffs(wc, resonance, s, &a, &k, &g, bilinear);
}


/*
ToDo:
-maybe make a version based on a cascade of 1st order highpasses
 -the mixing coeffs will then be different for the modes
 -we don't need gain-compensation for lowpass anyore (but then we need it for highpass)
-maybe make also versions based on a LP->HP->LP->HP chain
 -this may need a feedback factor with non-inverted sign
 -perhaps we should use BLT 1st order sections for this
 -maybe make a version based on 1st order allpass filters
-more generally, the stages could all have a different cutoff frequency and/or mode
 -maybe more flexible response types can be created by this (maybe shelf, peak?)
-Try to get arbitrary frequency responses out of it by forming the feedback signal also via a 
 linear combination, like (maybe) so: 
   y[0]  = in - (k[0]*y[0] + k[1]*y[1] + k[2]*y[2] + k[3]*y[3] + k[4]*y[4]);
 instead of:
   y[0]  = in - k*y[4];
 the coeffs for the k-values are found by equating:
           c0 + c1*s + c2*s^2 + c3*s^3 + c4*s^4
   H(s) = --------------------------------------  (not sure, if that's the right formula)
           k0 + k1*s + k2*s^2 + c3*s^3 + c4*s^4
 to some desired s-domain transfer function. But this would have to be a different design 
 procedure without the resonance parameter (i think). Or maybe it should be done in th z-domain.
 ...not yet sure

 -provide different morph modes:
   MORPH_LP_FLAT_HP:  LP_24, LP_18, LP_12, LP_6, FLAT, HP_6, HP_12, HP_18, HP_24
   MORPH_LP_BP_HP_24: LP_24, BP_6_18, BP_12_12, BP_18_6, HP_24
   MORPH_LP_BP_HP_18: LP_18, BP_6_12, BP_12_6, HP_18
   MORPH_LP_BP_HP_12: LP_12, BP_6_6, HP_12
   MORPH_LP_BP_HP:    LP_24, LP_18, LP_12, BP_6_12, BP_12_12, BP_12_6, HP_12, HP_18, HP_24

  for a 15 dB/oct lowpass, the prototype would be 1 / (1 + s^2.5), maybe we can approximate it
  resaonably by 1 / (1 + (s^2 + s^3)/2), see: https://www.desmos.com/calculator/gsbvrxiceb
  or: always use a 4th order Taylor series of the s^x term

*/
