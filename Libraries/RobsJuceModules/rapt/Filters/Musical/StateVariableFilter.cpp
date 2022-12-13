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
  // Compute intermediate values. The square rooths could be imaginary but when we form their 
  // quotient and product, it will become real:
  using Complex = std::complex<TPar>;  
  TPar    u1 = -TPar(1) - a1 - a2;   // could be negative
  TPar    u2 = -TPar(1) + a1 - a2;   // ...dito
  Complex s1 = sqrt(Complex(u1));    // could be imaginary
  Complex s2 = sqrt(Complex(u2));    // ...dito  
  TPar    p  = real(s1 * s2);        // but their product should be real
  TPar    s  = TPar(1) / p;          // we actually need the product's reciprocal
  
  // Compute coeffs:
  g  = real(s1 / s2);                         // 16a, the quotient should also be real
  R2 = s * TPar(2) * (a2 - TPar(1));          // 16b
  cH = (b0 - b1 + b2) / (TPar(1) - a1 + a2);  // 16c, == -(b0-b1+b2) / s1   before taking the sqrt?
  cB = s * TPar(2) * (b0 - b2);               // 16d
  cL = (b0 + b1 + b2) / (TPar(1) + a1 + a2);  // 16e


  cB *= -1.0;  // Test - fudging - yes - this indeed fixes a problem!

  h = 1 / (1 + R2*g + g*g);  // factor for feedback precomputation


  
  // formulas from (Eq 16 a-e):
  // http://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf
  
  // ToDo:
  // -Figure out why we need the factor -1 for the cB coeff with respect to the formula in the 
  //  paper
}

// Misc:

template<class TSig, class TPar>
void rsStateVariableFilter<TSig, TPar>::calcCoeffs()
{



  g = tan( TPar(PI) * fc/fs);  // embedded integrator gain (Fig 3.11)

  switch( mode )
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

-Implement formulas from at this paper 
 http://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf
 -> it has simpler formulas and even formulas that work for biquads with arbitrary coefficients
 -> make a function setupFromBiquad in the same way as in rsStateVectorFilter

-When done, maybe use these formulas for the biquad design:
 https://www.vicanek.de/articles/BiquadFits.pdf

-Maybe have a look at this, too:
 https://zrna.org/akso/object/contrib/tiar/filter/ZDF-SVF-1
 The coefficient computation seems to start from a UDF design (Chamberlin) and then (Newton?) 
 iterate to refine for the ZDF case? Figure out!

-See also:
 https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf



*/