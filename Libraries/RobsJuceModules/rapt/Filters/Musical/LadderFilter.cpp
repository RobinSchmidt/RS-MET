template<class TSig, class TPar>
rsLadderFilter<TSig, TPar>::rsLadderFilter()
{
  sampleRate = 44100.0;
  cutoff     = 1000.0;
  resonance  = 0.0;
  setMode(LP_24);
  updateCoefficients();
  reset();
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
  using T = TPar;
  if( newMode >= 0 && newMode < NUM_MODES )
  {
    mode = newMode;
    switch(mode)
    {
    case FLAT:     { setMixingCoefficients(1,  0,   0,   0,  0); s = 0.5;    } break;
      // check, if the s value makes sense

    // lowpasses:
    case LP_6:     { setMixingCoefficients(0,  1,   0,   0,  0); s = T(1);    } break;
    case LP_12:    { setMixingCoefficients(0,  0,   1,   0,  0); s = T(1);    } break;
    case LP_18:    { setMixingCoefficients(0,  0,   0,   1,  0); s = T(1);    } break;
    case LP_24:    { setMixingCoefficients(0,  0,   0,   0,  1); s = T(1);    } break;

    // highpasses:
    case HP_6:     { setMixingCoefficients(1, -1,   0,   0,  0); s = T(0);    } break;
    case HP_12:    { setMixingCoefficients(1, -2,   1,   0,  0); s = T(0);    } break;
    case HP_18:    { setMixingCoefficients(1, -3,   3,  -1,  0); s = T(0);    } break;
    case HP_24:    { setMixingCoefficients(1, -4,   6,  -4,  1); s = T(0);    } break;

    // bandpasses:
    case BP_6_18:  { setMixingCoefficients(0,  0,   0,   4, -4); s = T(0.75); } break; // 0.75
    case BP_12_12: { setMixingCoefficients(0,  0,   4,  -8,  4); s = T(0.5);  } break; // 0.5
    case BP_18_6:  { setMixingCoefficients(0,  4, -12,  12, -4); s = T(0.25); } break; // 0.25
    case BP_6_12:  { setMixingCoefficients(0,  0,   3,  -3,  0); s = T(2./3); } break; // 0.66
    case BP_12_6:  { setMixingCoefficients(0,  3,  -6,   3,  0); s = T(1./3); } break; // 0.33
    case BP_6_6:   { setMixingCoefficients(0,  2,  -2,   0,  0); s = T(0.5);  } break; // 0.5
    // generally: s = slope2 / (slope1 + slope2) where slope2 is the lowpass slope




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
  std::complex<TPar> G1, G2, G3, G4; // transfer functions of n-th stage output, n = 1..4
  std::complex<TPar> H;              // transfer function with resonance   
  std::complex<TPar> one(1, 0);

  G1 = b / (one + a/z);
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

  // Create rational function object and return it (note again the binomial coeffs 1,4,6,4,1):
  using RF = RAPT::rsRationalFunction<T>;
  RF H(N, {a2*a2, T(4)*a*a2, T(6)*a2, T(4)*a + b4*k, T(1) });
  return H;
}
// Maybe make a function getCoeffsBA(TPar* b, TPar* a) which just fills the arrays (of length 5).
// Beware that the order must be reversed for a polynomial in z^-1 (maybe make it optional to be
// in this form)

template<class TSig, class TPar>
rsRationalFunction<TPar> rsLadderFilter<TSig, TPar>::getTransferFunctionOld()
{
  using RF = RAPT::rsRationalFunction<TPar>;
  TPar tol = 1024 * RS_EPS(TPar);  // ad hoc
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
  //y[4] /= 1 + y[4]*y[4];   // (ad hoc) nonlinearity applied to the feedback signal
  y[0]  = in - k*y[4];        // linear
  y[0]  = rsClip(y[0], TSig(-1), TSig(+1));
  //y[0] /= TSig(1) + y[0]*y[0]; // nonlinearity applied to input plus feedback signal (division could be interesting with complex signals)
  //y[0]  = rsNormalizedSigmoids<TSig>::softClipHexic(y[0]);
  y[1]  = b*y[0]  - a*y[1];
  y[2]  = b*y[1]  - a*y[2];
  y[3]  = b*y[2]  - a*y[3];
  y[4]  = b*y[3]  - a*y[4];
  return c[0]*y[0] + c[1]*y[1] + c[2]*y[2] + c[3]*y[3] + c[4]*y[4];

  // We should experiment with placing saturation at different points - and, of course, try more 
  // realistic functions. The y = x / (1 + x^2) that is currently used is not even a proper
  // sigmoid shape (it goes down to zero again) ...but maybe that's not necessarily a bad thing.
  // Maybe try y = b + (x-b) / (1 + a*x^2) with adjustable parameters a and b. "a" adjusts the 
  // amount of signal squasheing whereas "b" introduces asymmetry.

  // maybe try a different update euqation of the form:
  // y[i] = y[i-1] + coeff * (y[i] - y[i-1]);
  // as is used in the rsSmoothingFilter. Elan says, this responds better to modulation. i think,
  // coeff = a or coeff = -a


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
TPar rsLadderFilter<TSig, TPar>::computeFeedbackFactor(CRPar fb, CRPar cosWc, CRPar a, CRPar b)
{
  TPar g2 = b*b / (1 + a*a + 2*a*cosWc);
  return fb / (g2*g2);
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
void rsLadderFilter<TSig, TPar>::computeCoeffs(CRPar wc, CRPar fb, CRPar s, TPar *a, TPar *b, 
  TPar *k, TPar *g)
{
  computeCoeffs(wc, fb, a, b, k);
  //*g = 1 + *k; // this overall gain factor ensures unit gain at DC regardless of resonance
               // damn! this gain works only for the lowpass case...was it always like that?
               // what about the old formula?
  //*g = TPar(1);

  *g = 1 + s * *k;

  // 1+k is good for all lowpasses, 1 is good for all highpasses, the allpass needs 
  // perhaps 1 + k/2

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
}

template<class TSig, class TPar>
void rsLadderFilter<TSig, TPar>::computeCoeffs(CRPar wc, CRPar fb, TPar *a, TPar *b, TPar *k)
{
  TPar s, c, t;                     // sin(wc), cos(wc), tan((wc-PI)/4)
  //rsSinCos(wc, &s, &c);
  s  = rsSin(wc);
  c  = rsCos(wc);
  t  = (TPar) rsTan(0.25*(wc-PI));
  *a = t / (s-c*t);
  *b = 1 + *a;  
  *k = computeFeedbackFactor(fb, c, *a, *b);
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
  computeCoeffs(wc, resonance, s, &a, &b, &k, &g);
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



*/
