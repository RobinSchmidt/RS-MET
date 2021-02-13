//=================================================================================================
// helper functions for damped-sine filter design:

template<class TPar, class TCof>
void rsDampedSineFilterCoeffs(
  TPar w, TPar A, TPar d, TPar p, TCof* b0, TCof* b1, TCof* a1, TCof* a2)
{
  TPar cw, sw, cp, sp, P;
  rsSinCos(w, &sw, &cw);
  rsSinCos(p, &sp, &cp);
  P   = exp(-1.0/d);              // = exp(-alpha), pole radius
  *a1 = TCof(-2*P*cw);            // = -2*P*cos(w)
  *a2 = TCof(P*P);                // = P^2
  *b0 = TCof(A*sp);               // = A*sin(p)
  *b1 = TCof(A*P*(sw*cp-cw*sp));  // = A*P*sin(w-p) via addition theorem
}
// i tried to bake the minus sign into the a-coeffs such that the difference equation can be 
// implemented with all plusses - but that didn't give any performance advantage, so i changed
// it back for consistency with DSP literature

template<class T>
T findDecayScalerLess1(T c)
{
  if(c <= 0.0 || c >= 1.0)
  {
    rsError("Function assumes 0 < c < 1");
    return 1.0;
  }

  // precomputations:
  T kp  = 1/c;                               // location of the the peak of g(k)
  T k   = 1 + 2*(kp-1);                      // initial guess for the zero of g(k)
  T eps = std::numeric_limits<T>::epsilon(); // relative tolerance
  int    i   = 0;                            // iteration counter

  // Newton iteration:
  T g, gp;      // g(k), g'(k)
  T kOld = 2*k; // ensure to enter the loop
  while(fabs(k-kOld) > k*eps && i < 1000)  // whoa! 1000? that seems way too high for production code!
  {
    kOld = k;
    g    = log(k) + c*(1-k); // g(k)
    gp   = 1/k - c;          // g'(k)
    k    = k - g/gp;         // Newton step
    i++;                     // count iteration
  }

  return k;

  // \todo: check this function in the range 0 <= c < 1, if all works well and the iteration count 
  // is always low, get rid of the iteration counter - it serves a purpose only during development

  // \todo find a refined formula for the initial guess by plotting the output against the input in
  // the range 0...1 and fit a polynomial (or other suitable function) to the data

  // ...i think, the goal should be that the Newton iteration converges in 2 or 3 steps - maybe we 
  // can then switch to a fixed number of steps, maybe we could also try Halley iteration instead,
  // where g''(k) = -1/k^2 - compute: kr = 1/k, gp = kr - c, gpp = -kr*kr
  // See:
  // https://en.wikipedia.org/wiki/Halley%27s_method
  // http://numbers.computation.free.fr/Constants/Algorithms/newton.html
}

template<class T>
void expDiffScalerAndTau2(T tau1, T tp, T* tau2, T* scaler)
{
  if(tp >= tau1)
  {
    rsError("assumes tp < tau1");
    *tau2   = tau1;
    *scaler = 1.0;
    return;
  }
  T a1 = 1/tau1;
  T c  = a1 * tp;
  T k  = findDecayScalerLess1(c);
  T a2 = k*a1;
  T hp = exp(-a1*tp) - exp(-a2*tp); // peak height

  *tau2   = 1/a2;
  *scaler = 1/hp;
}

//=================================================================================================
// class rsTwoPoleFilter:

template<class TSig, class TPar>
rsTwoPoleFilter<TSig, TPar>::rsTwoPoleFilter()
{
  a1 = a2 = 0;
  g  = 1;
  reset();
}

template<class TSig, class TPar>
void rsTwoPoleFilter<TSig, TPar>::setFrequencyAndDecay(TPar w, TPar d)
{
  TPar P = exp(-1/d);  // pole radius
  a1 = -2*P*cos(w);
  a2 = P*P;
}

template<class TSig, class TPar>
void rsTwoPoleFilter<TSig, TPar>::setOutputGain(TPar newGain)
{
  g = newGain;
}

template<class TSig, class TPar>
TPar rsTwoPoleFilter<TSig, TPar>::getMagnitudeAt(TPar w)
{
  return biquadMagnitudeAt(TPar(1), TPar(0), TPar(0), a1, a2, w);
  // optimize: a simpler (less general) formula may be used
  // ...i think, we should include the gain g - but before doing so, figure out, if the function
  // is used somewhere in order to compute and set a compensation gain - something like:
  //   g = filter.getMagnitudeAt(cutoff);
  //   filter.setOutputGain(1/g);
  // in such a case, this outlying code would break if we include the gain here.
  // maybe we should have two getMagnitudeAt functions - including and excluding the gain
}

template<class TSig, class TPar>
void rsTwoPoleFilter<TSig, TPar>::reset()
{
  y1 = y2 = 0;
}

//=================================================================================================
// class rsModalFilter:

template<class TSig, class TPar>
rsModalFilter<TSig, TPar>::rsModalFilter()
{
  b0 = 1.0;
  b1 = a1 = a2 = 0.0;
  reset();
}

template<class TSig, class TPar>
void rsModalFilter<TSig, TPar>::setModalParameters(TPar frequency, TPar amplitude, TPar decayTime, 
  TPar startPhase, TPar sampleRate)
{
  rsDampedSineFilterCoeffs(TPar(2*PI)*frequency/sampleRate, amplitude, decayTime*sampleRate,
    RAPT::rsDegreeToRadiant(startPhase), &b0, &b1, &a1, &a2);  
}

template<class TSig, class TPar>
void rsModalFilter<TSig, TPar>::copyCoefficientsFrom(const rsModalFilter &other)
{
  b0 = other.b0;
  b1 = other.b1;
  a1 = other.a1;
  a2 = other.a2;
}

template<class TSig, class TPar>
TPar rsModalFilter<TSig, TPar>::getDecayTime(TPar sampleRate)
{
  return -1.0 / (log(rsSqrt(a2))*sampleRate);
}

template<class TSig, class TPar>
std::complex<TPar> rsModalFilter<TSig, TPar>::getTransferFunctionAt(std::complex<TPar> z)
{
  return biquadTransferFunctionAt(b0, b1, TPar(0), a1, a2, z);
}

template<class TSig, class TPar>
TPar rsModalFilter<TSig, TPar>::getMagnitudeAt(TPar w)
{
  return biquadMagnitudeAt(b0, b1, TPar(0), a1, a2, w);
}

template<class TSig, class TPar>
void rsModalFilter<TSig, TPar>::reset()
{
  x1 = y1 = y2 = 0.0;
}

template<class TSig, class TPar>
void rsModalFilter<TSig, TPar>::processBlock(TSig in[], TSig out[], int blockSize)
{
  for(int n = 0; n < blockSize; n++)
    out[n] = getSample(in[n]);
}

/*
template<class TSig, class TPar>
void rsModalFilter<TSig, TPar>::processBlock(TSig in[], TSig out[], int blockSize)
{
  // under construction - still it is less efficient to call the block-based version compared to 
  // the sample-based version, also it is not yet entirely correct (the scale factor should be 
  // rendered into the state-variables - we need to swicth from g,b1 to b0,b1 coeffs

  out[0] = in[0] + b1*x1 - a1*y1 - a2*y2;
  x1 = in[0];  
  y2 = y1;
  y1 = out[0];

  if( blockSize < 2 )
    return;

  out[1] = in[1] + b1*x1 - a1*y1 - a2*y2;
  x1 = in[1];  
  y2 = y1;
  y1 = out[1];

  if( blockSize < 3 )
    return;

  for(int n = 2; n < blockSize; n++)
  {
    out[n] = in[n] + b1*in[n-1] - a1*out[n-1] - a2*out[n-2];
    //out[n] = g*y;  // we need to use b0, b1 instead of g, b1
  }

  x1 = in[blockSize-1];
  y2 = out[blockSize-2];
  y1 = out[blockSize-1];

  //scale(out, blockSize, g);


  //for(int n = 0; n < blockSize; n++)
  //{
  //  //out[n] = getSample(in[n]);
  //  double y = in[n] + b1*x1 - a1*y1 - a2*y2;
  //  x1 = in[n];
  //  y2 = y1;
  //  y1 = y;
  //  out[n] = g*y;
  //}

}
  */


//=================================================================================================
// class rsNonlinearModalFilter:

template<class TSig, class TPar>
rsNonlinearModalFilter<TSig, TPar>::rsNonlinearModalFilter()
{  
  a  = 0.0; 
  cr = 1.0;
  ci = 0.0;

  startPhase    = 0.0;
  amplitude     = 1.0;
  phaseModByIm  = 0.0;
  phaseModByAbs = 0.0;

  reset();
}

template<class TSig, class TPar>
void rsNonlinearModalFilter<TSig, TPar>::setModalParameters(TPar frequency, TPar amplitude, 
  TPar decayTime, TPar startPhase, TPar sampleRate)
{
  // compute recursion coeff:
  TPar w     = 2*PI*frequency/sampleRate;
  TPar alpha = 1.0 / (decayTime*sampleRate); 
  TPar r     = exp(-alpha);
  //a.setRadiusAndAngle(r, w); // write a function setPolar...or check if there already is one in std::complex
  a = std::polar(r, w);

  this->amplitude  = amplitude;
  this->startPhase = startPhase;

  // compute output weights:
  rsSinCos((PI/180.0)*startPhase, &cr, &ci);
  //cr *= amplitude;
  //ci *= amplitude;
}

template<class TSig, class TPar>
void rsNonlinearModalFilter<TSig, TPar>::setAmplitude(TPar newAmplitude)
{
  amplitude = newAmplitude;
}

template<class TSig, class TPar>
void rsNonlinearModalFilter<TSig, TPar>::setPhaseModulation(TPar newPhaseModulation)
{
  // todo: rename into phasemodualtionByIm, try to find a mor intuitive parametrization in terms
  // of the overtone spectrum (maybe a formula may be derived that relates a spectral slope to the 
  // modulation index or something) ...or maybe a slope (in dB/oct) can be related to the overall
  // level in dB ...so we may dial in a value in (dB/oct)/dB = dB^2/oct
  phaseModByIm = newPhaseModulation;
}

template<class TSig, class TPar>
void rsNonlinearModalFilter<TSig, TPar>::copyCoefficientsFrom(const rsNonlinearModalFilter &other)
{
  a  = other.a;
  cr = other.cr;
  ci = other.ci;
}

template<class TSig, class TPar>
void rsNonlinearModalFilter<TSig, TPar>::reset()
{
  z = 0.0;
}

//=================================================================================================
// class ModalFilterWithAttack:

template<class TSig, class TPar>
void rsModalFilterWithAttack<TSig, TPar>::setModalParameters(TPar frequency, TPar amplitude, 
  TPar attackTime, TPar decayTime, TPar startPhase, TPar sampleRate, TPar detuneFactor)
{
  rsAssert(attackTime < decayTime);  // attackTime >= decayTime will not work (because of math)

  TPar tau1, tau2, scaler;
  tau1 = decayTime;
  expDiffScalerAndTau2(tau1, attackTime, &tau2, &scaler);
  amplitude *= scaler;
  modalFilter1.setModalParameters(frequency,              amplitude, tau1, startPhase, sampleRate);
  modalFilter2.setModalParameters(frequency*detuneFactor, amplitude, tau2, startPhase, sampleRate);

  // \todo "detuneFactor" does not really work well because when attack and decay are very similar, 
  // the amplitude explodes (due to a high value of "scaler") - either remove this parameter or 
  // find a way to alleviate this (maybe the amplitude excess can be computed - math has to be 
  // worked out) the detune is supposed to introduce some roughness into the transient by detuning 
  // the quickly decaying sinusoid
}

template<class TSig, class TPar>
std::complex<TPar> rsModalFilterWithAttack<TSig, TPar>::getTransferFunctionAt(std::complex<TPar> z)
{
  return modalFilter1.getTransferFunctionAt(z) + modalFilter2.getTransferFunctionAt(z);
}

template<class TSig, class TPar>
TPar rsModalFilterWithAttack<TSig, TPar>::getLength(TPar decayLevel, TPar sampleRate)
{
  TPar td = modalFilter1.getDecayTime(sampleRate);       // decay time
  TPar a1 = 1.0 / td;
  TPar a2 = 1.0 / modalFilter2.getDecayTime(sampleRate);
  TPar tp = (log(a1)-log(a2))/(a1-a2);                   // peak time
  return tp + rsTauToDecayTime(td, decayLevel);
}

template<class TSig, class TPar>
void rsModalFilterWithAttack<TSig, TPar>::reset()
{
  modalFilter1.reset();
  modalFilter2.reset();
}

//=================================================================================================
// class ModalFilterWithAttack2:

template<class TSig, class TPar>
rsModalFilterWithAttack2<TSig, TPar>::rsModalFilterWithAttack2()
{
  //b0 = 1.0;
  b1 = b2 = b3 = a1 = a2 = a3 = a4 = 0.0;
  reset();
}

template<class TSig, class TPar>
void rsModalFilterWithAttack2<TSig, TPar>::setModalParameters(TPar frequency, TPar amplitude, 
  TPar attackTime, TPar decayTime, TPar startPhase, TPar sampleRate)
{
  rsAssert(attackTime < decayTime);  // attackTime >= decayTime will not work (because of math)

  TPar tau1, tau2, scaler;
  TPar w = 2*PI*frequency/sampleRate;
  TPar p = RAPT::rsDegreeToRadiant(startPhase);
  tau1 = decayTime;
  expDiffScalerAndTau2(tau1, attackTime, &tau2, &scaler);

  // compute coefficients for a parallel connection of filters:
  TPar a = amplitude * scaler;
  TPar b10, b11, a11, a12;
  TPar b20, b21, a21, a22;
  rsDampedSineFilterCoeffs(w, a, tau1*sampleRate, p, &b10, &b11, &a11, &a12);
  rsDampedSineFilterCoeffs(w, a, tau2*sampleRate, p, &b20, &b21, &a21, &a22);

  // convert parallel connection into a single 4th order filter (maybe factor out):
  //b0 = b10-b20; // b0 comes out as zero
  b1 = b11-b21+a21*b10-a11*b20;
  b2 = a21*b11-a11*b21+a22*b10-a12*b20;
  b3 = a22*b11-a12*b21;
  a1 = a21+a11;
  a2 = a22+a11*a21+a12;
  a3 = a11*a22+a12*a21;
  a4 = a12*a22;
}

template<class TSig, class TPar>
void rsModalFilterWithAttack2<TSig, TPar>::reset()
{
  w1 = w2 = w3 = w4 = 0.0;
}

//=================================================================================================
// class ModalFilterBank:

// construction/destruction:

template<class TSig, class TPar>
rsModalFilterBank<TSig, TPar>::rsModalFilterBank()
{
  sampleRate         = 44100.0;
  referenceFrequency = 440.0;
  referenceAttack    = 0.1;
  referenceDecay     = 1.0;
  numModes           = maxNumModes;
  modalFilters.reserve(maxNumModes);
  for(int m = 0; m < maxNumModes; m++)
    modalFilters.push_back(rsModalFilterWithAttack<TSig, TPar>());
}

template<class TSig, class TPar>
rsModalFilterBank<TSig, TPar>::~rsModalFilterBank()
{

}

// parameter settings:

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  calculateModalFilterCoefficients();
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setReferenceFrequency(TPar newFrequency)
{
  referenceFrequency = newFrequency;
  calculateModalFilterCoefficients();
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setReferenceAttack(TPar newAttack)
{
  referenceAttack = newAttack;
  calculateModalFilterCoefficients();
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setReferenceDecay(TPar newDecay)
{
  referenceDecay = newDecay;
  calculateModalFilterCoefficients();
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setModalParameters(std::vector<TPar> newFrequencies, 
  std::vector<TPar> newAmplitudes, std::vector<TPar> newAttackTimes, 
  std::vector<TPar> newDecayTimes, std::vector<TPar> newStartPhases)
{
  frequencies = newFrequencies;
  amplitudes  = newAmplitudes;
  attackTimes = newAttackTimes;
  decayTimes  = newDecayTimes;
  startPhases = newStartPhases;
  calculateModalFilterCoefficients();
}

// inquiry:

template<class TSig, class TPar>
std::complex<TPar> rsModalFilterBank<TSig, TPar>::getTransferFunctionAt(std::complex<TPar> z)
{
  std::complex<TPar> H = std::complex<TPar>(0, 0); // accumulator for H(z)
  for(int m = 0; m < getNumModes(); m++)
    H += modalFilters[m].getTransferFunctionAt(z);
  return H;
  // maybe later rename H to G and use H for transfer function with feedback, G for without
}

template<class TSig, class TPar>
TPar rsModalFilterBank<TSig, TPar>::getLength(TPar decayLevel)
{
  TPar max = 0.0;
  TPar tmp;
  for(size_t m = 0; m < decayTimes.size(); m++)
  {
    tmp = modalFilters[m].getLength(decayLevel, sampleRate);
    if( tmp > max )
      max = tmp;
  }
  return max;
  // actually, we could just look it up in our attckTime/decayTime members instead of letting the
  // filters compute their lengths
}

// audio processing:

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::processBlock(TSig in[], TSig out[], int blockSize)
{
  for(int n = 0; n < blockSize; n++)
    out[n] = getSample(in[n]);

  // \todo optimize this
}

// others:

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::reset()
{
  for(size_t m = 0; m < modalFilters.size(); m++)
    modalFilters[m].reset();
  out = TSig(0);
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::calculateModalFilterCoefficients()
{
  size_t nm = rsMin((size_t)numModes, frequencies.size(), amplitudes.size(), decayTimes.size());
  nm = rsMin(nm, startPhases.size());
  for(size_t m = 0; m < nm; m++)
  {
    modalFilters[m].setModalParameters(
      referenceFrequency * frequencies[m], 
      amplitudes[m],
      referenceAttack * attackTimes[m],
      referenceDecay * decayTimes[m], 
      startPhases[m], 
      sampleRate); 
  }
}

// static member functions:

template<class TSig, class TPar>
TPar rsModalFilterBank<TSig, TPar>::modeDecayTime(TPar f, TPar fc, TPar p)
{
  TPar k = pow(fc, -p);
  return (1-k) / ((1-k*2) + pow(f/fc, p));
}

template<class TSig, class TPar>
std::vector<TPar> rsModalFilterBank<TSig, TPar>::randomModePhases(
  const std::vector<TPar> &a, TPar randomness, int seed)
{
  // The formulas used here were obtained by requiring a1*sin(p1) + a2*sin(p2) = k and solving 
  // for p2. For the last value in case of an odd number of modes, it is: a*sin(p) = k.
  // In order to let one mode cancel the other at n = 0, we use k = 0. We need to ensure that the 
  // mode for which we adjust the phase has at least the same amplitude as the mode for which we
  // assume a fixed phase in order to be able to cancel in all cases for the phase - this is why 
  // there is this if(a[n] < a[n+1]) thing.

  rsRandomUniform(0.0, 1.0, seed);
  std::vector<TPar> p(a.size());
  TPar k = 0.0;    // target value for the sample at n = 0
  int    N = (int) p.size();  // upper limit for loop below
  if( rsIsOdd(p.size()) )
  {
    if( a[p.size()-1] == 0.0 )
      p[p.size()-1] = 0.0;
    else
      p[p.size()-1] = asin(k/a[p.size()-1]);
    N--;
  }
  for(int n = 0; n < N; n += 2)
  {
    if( fabs(a[n]) == 0.0 && fabs(a[n+1]) == 0.0 )
      continue;
    if( fabs(a[n]) < fabs(a[n+1]) )
    {
      p[n]   = rsRandomUniform(0.0, 2*PI*randomness);
      p[n+1] = asin((k-a[n]*sin(p[n]))/a[n+1]);
    }
    else // roles reversed
    {
      p[n+1] = rsRandomUniform(0.0, 2*PI*randomness);
      p[n]   = asin((k-a[n+1]*sin(p[n+1]))/a[n]);      
    }
  }
  rsArrayTools::scale(&p[0], (int)p.size(), 360.0/(2*PI));
  return p;
}

template<class TSig, class TPar>
std::vector<TPar> rsModalFilterBank<TSig, TPar>::modeDecayTimes(std::vector<TPar> f, TPar fc, TPar p)
{
  std::vector<TPar> d(f.size());
  for(size_t n = 0; n < d.size(); n++)
    d[n] = modeDecayTime(f[n], fc, p);
  return d;
}

template<class TSig, class TPar>
std::vector<TPar> rsModalFilterBank<TSig, TPar>::scaleAtIntervals(std::vector<TPar> v,
  int startIndex, int interval, TPar scaler)
{
  std::vector<TPar> r = v;
  for(size_t n = startIndex; n < r.size(); n += interval)
    r[n] *= scaler;
  return r;
}

/*
Idea 

use an adjustable mix of different imput signals: 
-unit impulse: pluck/strike, maybe use other kinds of impulse-like signals
-white noise: blowing, maybe use different kinds of noises - different probability densities,
 maybe coloring, maybe a bimodal density could be interesting as well

modeling a scraping input signal:
-use impulse train with adjustable random jitter
-amplitude for each pulse may also be randomized and/or be a function of the time passed since the
 previous pulse in order to normalize the energy over time (denser pulses should be more quiet)
-maybe use a box-filter to give each pulse a width - maybe that width can also be randomized (then
 the amplitude normalizer should also take into account the width)
-maybe use pairs or triples of pulses
-maybe add two or more of such scrape models


modeling transients:
-transients are modeled as superposition attack/decay envelope filters (i.e. zero frequency) with
 delay
-maybe an interative matching procedure can be used to find the parameters of the filters

-or maybe use modeling in terms of a general pole/zero model

*/