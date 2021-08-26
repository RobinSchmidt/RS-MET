template<class TSig, class TPar>
rsLadderFilter<TSig, TPar>::rsLadderFilter()
{
  reset();
  setMode(Mode::LP_24);
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::setMode(int newMode)
{
  // Shorthands for convenience:
  using T = TPar;
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
rsComplex<TPar> rsLadderFilter<TSig, TPar>::getTransferFunctionAt(
  const rsComplex<TPar>& z, bool withGain)
{
  using T = TPar;
  TPar B0 = TPar(1) - B1; 
  rsComplex<T> G1, G2, G3, G4; // transfer functions of n-th stage output, n = 1..4
  rsComplex<T> H;              // transfer function with resonance   
  rsComplex<T> one(1, 0);
  G1 =  (1+a)*(B0 + B1/z) / (one + a/z);
  G2 = G1*G1;
  G3 = G2*G1;
  G4 = G3*G1;
  H  = (c[0] + c[1]*G1 + c[2]*G2 + c[3]*G3 + c[4]*G4) / (one + k * G4/z);
  if(withGain) return g*H;
  else         return   H;
}

template<class TSig, class TPar>
TPar rsLadderFilter<TSig, TPar>::getMagnitudeResponseAt(CRPar frequency, bool withG)
{
  using Cmp = rsComplex<TPar>;
  TPar w = 2 * TPar(PI) * frequency/sampleRate;
  Cmp j(0, 1);                             // imaginary unit
  Cmp z = rsExp(j*w);                      // location in the z-plane
  Cmp H = getTransferFunctionAt(z, withG); // H(z) at our z
  H *= rsConj(H);                          // magnitude-squared
  return rsSqrt(H.real());                 // imaginary part should be zero anyway

  // Why not just return rsAbs(H)?

  // I wonder, if a simpler formula is possible which avoids going through the complex transfer 
  // function -> computer algebra
}

// audio processing:

template<class TSig, class TPar>
inline TSig rsLadderFilter<TSig, TPar>::getSampleNoGain(CRSig in)
{
  // Apply feedback and saturation:
  //y[4] /= 1 + y[4]*y[4];   // (ad hoc) nonlinearity applied to the feedback signal
  // try: y[4] /= 1 / (y[4]*y[4])^n, see https://www.desmos.com/calculator/rjlnhzqecs
  // (x+x^3)/(1+x^4)  https://www.desmos.com/calculator/wb39hmukry
  // (x+0.25*x^3)/(1+x^4)  ...very straight in the middle, goes through (0.5,0.5)
  // (x+0.25*x^3)/(1+x^2)

  TSig tmp1 = y[0];           // needed for bilinear version only
  TSig tmp2;

  y[0]  = in - k*y[4];        // linear
  //y[0]  = rsClip(y[0], TSig(-1), TSig(+1));  // cheapest
  y[0] /= TSig(1) + y[0]*y[0]; // nonlinearity applied to input plus feedback signal (division could be interesting with complex signals)
  //y[0]  = rsNormalizedSigmoids<TSig>::softClipHexic(y[0]);  // most expensive
  // ToDo: let the user select the saturationMode: hardClip, tanh, softClip, etc. maybe with an 
  // optional DC (add before the saturation, subtract after) ..do that in subclasses - here, we just do the 
  // cheap hardclipping

  // Apply the 1st order stages:
  if(B1 != TPar(0)) 
  {
    // general case:
    TPar B0 = TPar(1) - B1;
    TPar b0 = (1+a)*B0;
    TPar b1 = (1+a)*B1;
    tmp2 = y[1]; y[1] = b0*y[0] + b1*tmp1 - a*y[1];  // tmp1 = old y[0]
    tmp1 = y[2]; y[2] = b0*y[1] + b1*tmp2 - a*y[2];  // tmp2 = old y[1]
    tmp2 = y[3]; y[3] = b0*y[2] + b1*tmp1 - a*y[3];  // tmp1 = old y[2]
                 y[4] = b0*y[3] + b1*tmp2 - a*y[4];  // tmp2 = old y[3]
  }
  else 
  {
    // allpole case:
    y[1] = y[0] + a * (y[0] - y[1]);
    y[2] = y[1] + a * (y[1] - y[2]);
    y[3] = y[2] + a * (y[2] - y[3]);
    y[4] = y[3] + a * (y[3] - y[4]);
  }

  // Form linear combination of the stage outputs:
  return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];

  // We should experiment with placing saturation at different points - and, of course, try more 
  // realistic functions. The y = x / (1 + x^2) that is currently used is not even a proper
  // sigmoid shape (it goes down to zero again) ...but maybe that's not necessarily a bad thing.
  // Maybe try y = b + (x-b) / (1 + a*x^2) with adjustable parameters a and b. "a" adjusts the 
  // amount of signal squasheing whereas "b" introduces asymmetry. This is different from just
  // subtracting b before the nonlinearity and adding it back after it, because the x^2 in 
  // denominator will be either affected or not

  // ToDo: test reso-modulation. In this case, it's particluarly interesting, if applying the 
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
TPar rsLadderFilter<TSig, TPar>::computeFeedbackFactor(
  CRPar fb, CRPar cosWc, CRPar a, TPar b1)
{
  TPar g2;
  if(b1 == TPar(0)) {  // simple allpole case
    TPar b = TPar(1) + a;
    g2 = b*b / (1 + a*a + 2*a*cosWc); 
  }
  else {               // general case
    TPar b0 = TPar(1) - b1;
    b0 *= (1+a);
    b1 *= (1+a);
    g2 = (b0*b0 + b1*b1 + 2*b0*b1*cosWc) / (1 + a*a + 2*a*cosWc);
    // optimize: compute reciprocal of g2....
  }
  return fb / (g2*g2); //...then this can be turned into a multiplication
}

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
  TPar *k, TPar *g, CRPar b1)
{
  computeCoeffs(wc, fb, a, k, b1);
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
void rsLadderFilter<TSig, TPar>::computeCoeffs(CRPar wc, CRPar fb, TPar *a, TPar *k, CRPar b1)
{
  TPar s, c, t;                     // sin(wc), cos(wc), tan((wc-PI)/4)
  //rsSinCos(wc, &s, &c);
  s  = rsSin(wc);
  c  = rsCos(wc);
  t  = (TPar) rsTan(0.25*(wc-PI));

  if(b1 == TPar(0))    // allpole case
    *a = t / (s-c*t);
  else {               // general case
    TPar b0 = TPar(1)-b1;
    *a = (b0*t + b1*c*t + b1*s) / (-b0*c*t + b0*s - b1*t); 
  }
  *k = computeFeedbackFactor(fb, c, *a, b1);
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
// In the general case, that may be in between regular (b0=1,b1=0) and bilinear (b0=b1=1), we may 
// choose b0,b1 arbitrarily and would compute a1 as:
//   a1 = (b0 t + b1 c t + b1 s) / (-b0 c t + b0 s - b1 t)
// maybe implement that, too. We could get rid of the conditionals

// internal functions:

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::updateCoefficients()
{
  TPar wc = 2 * (TPar)PI * cutoff / sampleRate;
  computeCoeffs(wc, resonance, s, &a, &k, &g, B1);
}


/*
ToDo:
-maybe make a version based on a cascade of 1st order highpasses and/or allpasses
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

-maybe rename to rsLadderFilterUDF and make a similar class for a ZDF ladder, maybe factor out
 a common baseclass

-Make the resonance adjustable in dB because this is what the .sfz format expects:
 https://sfzformat.com/tutorials/basic_sfz_file
 For this, we need a formula for the feedback-factor in terms of the resonance in dB. This should 
 take into account also the makeup gain, so there may be quite some algebra to churn through. At 
 the end, a formula should result that can be implemented as static member function just like
 resonanceDecayToFeedbackGain, maybe resonanceLevelToFeedbackGain

*/
