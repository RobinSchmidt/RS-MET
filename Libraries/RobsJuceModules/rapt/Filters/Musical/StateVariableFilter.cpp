// Construction/Destruction:

template<class TSig, class TPar>
rsStateVariableFilter<TSig, TPar>::rsStateVariableFilter()
{
  fs   = 44100.0;
  fc   = 1000.0;
  mode = LOWPASS; 
  G    = (TPar)SQRT2_INV;
  B    = 2.0;
  m    = 0.0;
  calcCoeffs();
  reset(); 
}

// Setup:

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  fs = newSampleRate;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setMode(int newMode)
{
  mode = newMode;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setFrequency(TPar newFrequency)
{
  fc = newFrequency;
  calcCoeffs();
}
 
template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setGain(TPar newGain)
{
  G = newGain;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setBandwidth(TPar newBandwidth)
{
  B = newBandwidth;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setMorph(TPar newMorph)
{
  m = newMorph;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::setupFromBiquad(
  CRPar b0, CRPar b1, CRPar b2, CRPar a1, CRPar a2)
{
  // Compute intermediate values. The square roots could be imaginary but when we form their 
  // quotient and product, it will become real:
  using Complex = std::complex<TPar>;  
  TPar    u1 = -TPar(1) - a1 - a2;     // could be negative
  TPar    u2 = -TPar(1) + a1 - a2;     // ...dito
  // rsAssert(u1*u2 >= 0);             // triggers in one of the unit tests
  Complex s1 = sqrt(Complex(u1));      // could be imaginary
  Complex s2 = sqrt(Complex(u2));      // ...dito  
  TPar    p  = real(s1 * s2);          // but their product should be real
  TPar    s  = TPar(1) / p;            // we actually need the product's reciprocal

  // Compute coeffs:
  g  = real(s1 / s2);                         // 16a, the quotient should also be real
  R2 = s * TPar(2) * (a2 - TPar(1));          // 16b
  cH = (b0 - b1 + b2) / (TPar(1) - a1 + a2);  // 16c, == -(b0-b1+b2) / s1 before taking the sqrt?
  cB = s * TPar(2) * (b2 - b0);               // 16d, but with a factor of -1 (why?)
  cL = (b0 + b1 + b2) / (TPar(1) + a1 + a2);  // 16e
  h  = 1 / (1 + R2*g + g*g);                  // factor for feedback precomputation

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

// Misc:

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::calcCoeffs()
{
  g = tan(TPar(PI) * fc/fs);  // embedded integrator gain (Fig 3.11), == tan(wc/2) I think

  switch(mode)
  {
  case BYPASS:
  {
    R2 = 1 / G;  // can we use an arbitrary value here, for example R2 = 1?
    cL = 1;
    cB = getBandpassScaler();
    cH = 1;
  }
  break;
  case LOWPASS:
  {
    R2 = 1 / G;
    cL = 1; cB = 0; cH = 0;
  }
  break;
  case HIGHPASS:
  {
    R2 = 1 / G;
    cL = 0; cB = 0; cH = 1;
  }
  break;
  case BANDPASS_SKIRT:
  {
    R2 = 1 / G;
    cL = 0; cB = 1; cH = 0;
    // why is the bandwidth-parameter not used here?
  }
  break;
  case BANDPASS_PEAK:
  {
    R2 = 2*bandwidthToR(B);
    cL = 0; cB = R2; cH = 0;
  }
  break;
  case BANDREJECT:
  {
    R2 = 2*bandwidthToR(B);
    cL = 1; cB = 0; cH = 1;
  }
  break;
  case BELL:
  {
    TPar fl = TPar(fc*pow(2, -B/2)); // lower bandedge frequency (in Hz)
    TPar wl = TPar(tan(PI*fl/fs));   // warped radian lower bandedge frequency /(2*fs)
    TPar r  = g/wl; r *= r;          // warped frequency ratio wu/wl == (wc/wl)^2 where wu is the 
                                     // warped upper bandedge, wc the center
    R2 = TPar(2*sqrt(((r*r+1)/r-2)/(4*G)));
    cL = 1; cB = R2*G; cH = 1;
  }
  break;
  case LOWSHELF:
  {
    TPar A = sqrt(G);
    g /= sqrt(A);                    // scale SVF-cutoff frequency for shelvers
    R2 = TPar(2*sinh(B*log(2.0)/2));
    cL = G; cB = R2*A; cH = 1;
  }
  break;
  case HIGHSHELF:
  {
    TPar A = sqrt(G);
    g *= sqrt(A);                    // scale SVF-cutoff frequency for shelvers
    R2 = TPar(2*sinh(B*log(2.0)/2));
    cL = 1; cB = R2*A; cH = G;
  }
  break;
  case ALLPASS:
  {
    R2 = 2*bandwidthToR(B);
    cL = 1; cB = -R2; cH = 1;
  }
  break;


  // new/under construction:
  case LowpassMVS:
  {
    //rsError("Not yet working");
    //R2 = 2*bandwidthToR(B);       // (R2 == 2*R == 1/Q) ....wrong?
    R2 = 1/G;
    TPar w0 = TPar(2*PI) * fc/fs;
    TPar b0, b1, b2, a1, a2;
    rsFilterDesignFormulas::mvLowpassSimple(w0, 1/R2, &b0, &b1, &b2, &a1, &a2);
    setupFromBiquad(b0, b1, b2, a1, a2);
  }
  break;
  // Needs tests - seems to lead to unstable filters -  




  // experimental - maybe we must find better curves for cL, cB, cH:
  case MORPH_LP_BP_HP:
  {
    R2 = 1 / G;
    TPar x  = 2*m-1;

    //double x2 = x*x;
    //cL = 0.5*(x2-x); cB = 1-x2; cH = 0.5*(x2+x); // nah - not good

    // better:
    cL = rsMax(-x, TPar(0)); cH = rsMax(x, TPar(0)); cB = 1-(cL+cH);
    cB = pow(cB, TPar(0.25));
      // freq-responses look good (on a linear scale), but we really have to check how it "feels" 
      // it would also be nice to get rid of the expensive pow-function and to replace it by 
      // something cheaper - the function should map the range 0...1 monotonically to itself

    // another (cheap) possibility:
    //cL = rsMax(-x, 0.0); /*cL *= cL;*/
    //cH = rsMax( x, 0.0); /*cH *= cH;*/
    //cB = 1-x*x;

      // bottom line: we need to test different versions for how they feel when tweaking the 
      // morph parameter

    // this scaling ensures constant magnitude at the cutoff point (we divide the coefficients by 
    // the magnitude response value at the cutoff frequency and scale back by the gain):
    TPar s = G * sqrt((R2*R2) / (cL*cL + cB*cB + cH*cH - 2*cL*cH));
    cL *= s; cB *= s; cH *= s;
  }
  break;

  }

  h = 1 / (1 + R2*g + g*g);  // factor for feedback precomputation
}

template<class TSig, class TPar>
TPar rsStateVariableFilter<TSig, TPar>::bandwidthToR(TPar B)
{
  TPar fl = fc*pow(TPar(2), TPar(-B/2)); // lower bandedge frequency (in Hz)
  TPar gl = tan(TPar(PI)*fl/fs);         // warped radian lower bandedge frequency /(2*fs)
  TPar r  = gl/g;            // ratio between warped lower bandedge- and center-frequencies
                             // unwarped: r = pow(2, -B/2) -> approximation for low
                             // center-frequencies
  return sqrt((1-r*r)*(1-r*r)/(4*r*r));
}

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::reset()
{
  s1 = s2 = 0.0;
}



/*

ToDo:

-Implement a getTransferFunctionAt(Complex z) function with the same API as corresponding functions
 in the Ladder and Biquad filters

-Implement formulas from at this paper (done)
 http://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf
 -> it has simpler formulas and even formulas that work for biquads with arbitrary coefficients
 -> make a function setupFromBiquad in the same way as in rsStateVectorFilter

-When done, maybe use these formulas for the biquad design:
 https://www.vicanek.de/articles/BiquadFits.pdf
 OK - the formulas have been implemented in biquadDesignVicanek in FilterExperiments.cpp. But maybe 
 actually using those to first compute biquad coeffs and then convert from biquad to svf coeffs is 
 unnecessarily expensive. Try to replicate Martin's approach using an expression for the SVF 
 magnitude response directly. I think, the R2 and g variables can be computed from Q and w0 and it 
 should be quite easy to derive the expressions. It then remains to derive expressions for cH, cB, 
 cL which may be a bit more involved - we'll see....

-Maybe have a look at this, too:
 https://zrna.org/akso/object/contrib/tiar/filter/ZDF-SVF-1
 The coefficient computation seems to start from a UDF design (Chamberlin) and then (Newton?) 
 iterate to refine for the ZDF case? Figure out!

-This paper approaches the same filter from a circuit-modeling perspective:
 https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf




*/