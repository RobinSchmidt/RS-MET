// Construction/Destruction:

StateVariableFilter::StateVariableFilter()
{
  fs   = 44100.0;
  fc   = 1000.0;
  mode = LOWPASS; 
  G    = SQRT2_INV;
  B    = 2.0;
  m    = 0.0;
  calcCoeffs();
  reset(); 
}

// Setup:

void StateVariableFilter::setSampleRate(double newSampleRate)
{
  fs = newSampleRate;
  fc = rsClipToRange(fc, 0.0, 0.5*fs);
  calcCoeffs();
}

void StateVariableFilter::setMode(int newMode)
{
  mode = newMode;
  calcCoeffs();
}

void StateVariableFilter::setFrequency(double newFrequency)
{
  fc = rsClipToRange(newFrequency, 0.0, 0.5*fs);
  calcCoeffs();
}
 
void StateVariableFilter::setGain(double newGain)
{
  G = rsClipToRange(newGain, 0.0, 100.0);
  calcCoeffs();
}

void StateVariableFilter::setBandwidth(double newBandwidth)
{
  B = rsClipToRange(newBandwidth, 0.0, 100.0);
  calcCoeffs();
}

void StateVariableFilter::setMorph(double newMorph)
{
  m = rsClipToRange(newMorph, 0.0, 1.0);
  calcCoeffs();
}

// Misc:

void StateVariableFilter::calcCoeffs()
{
  // \todo look at this paper - it has simpler formulas and even formulas that work for biquads 
  // with arbitrary coefiicients:
  // http://www.dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf


  g = tan(PI*fc/fs);  // embedded integrator gain (Fig 3.11)

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
      double fl = fc*pow(2, -B/2); // lower bandedge frequency (in Hz)
      double wl = tan(PI*fl/fs);   // warped radian lower bandedge frequency /(2*fs)
      double r  = g/wl; r *= r;    // warped frequency ratio wu/wl == (wc/wl)^2 where wu is the 
                                   // warped upper bandedge, wc the center
      R2 = 2*sqrt(((r*r+1)/r-2)/(4*G));
      cL = 1; cB = R2*G; cH = 1;
    }
    break;
  case LOWSHELF:
    {
      double A = sqrt(G);
      g /= sqrt(A);               // scale SVF-cutoff frequency for shelvers
      R2 = 2*sinh(B*log(2.0)/2);
      cL = G; cB = R2*A; cH = 1;
    }
    break;
  case HIGHSHELF:
    {
      double A = sqrt(G);
      g *= sqrt(A);               // scale SVF-cutoff frequency for shelvers
      R2 = 2*sinh(B*log(2.0)/2);
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
      double x  = 2*m-1;

      //double x2 = x*x;
      //cL = 0.5*(x2-x); cB = 1-x2; cH = 0.5*(x2+x); // nah - not good

      // better:
      cL = rsMax(-x, 0.0); cH = rsMax(x, 0.0); cB = 1-(cL+cH);
      cB = pow(cB, 0.25);
        // freq-responses look good (on a linear scale), but we really have to check how it "feels" - it
        // would also be nice to get rid of the expensive pow-function and to replace it by 
        // something cheaper - the function should map the range 0...1 monotonically to itself

      // another (cheap) possibility:
      //cL = rsMax(-x, 0.0); /*cL *= cL;*/
      //cH = rsMax( x, 0.0); /*cH *= cH;*/
      //cB = 1-x*x;
      
        // bottom line: we need to test different versions for how they feel when tweaking the 
        // morph parameter

      // this scaling ensures constant magnitude at the cutoff point (we divide the coefficients by 
      // the magnitude response value at the cutoff frequency and scale back by the gain):
      double s = G * sqrt((R2*R2) / (cL*cL + cB*cB + cH*cH - 2*cL*cH));
      cL *= s; cB *= s; cH *= s;
    }
    break;

  }

  h = 1 / (1 + R2*g + g*g);  // factor for feedback precomputation
}

double StateVariableFilter::bandwidthToR(double B)
{
  double fl = fc*pow(2, -B/2); // lower bandedge frequency (in Hz)
  double gl = tan(PI*fl/fs);   // warped radian lower bandedge frequency /(2*fs)
  double r  = gl/g;            // ratio between warped lower bandedge- and center-frequencies
                               // unwarped: r = pow(2, -B/2) -> approximation for low
                               // center-frequencies
  return sqrt((1-r*r)*(1-r*r)/(4*r*r));
}

void StateVariableFilter::reset()
{
  s1 = s2 = 0.0;
}
