//=================================================================================================
// helper functions 

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

    // ToDo: Make this a static member function of rsModalFilterWithAttack 

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
  if(tp >= tau1)    // ToDo: Figure out if we need a tolerance! Like:  tp >= tau1 * (1-tol)
  {
    rsError("Assumes tp < tau1"); 
    *tau2   = tau1;
    *scaler = 1.0;
    return;
  }
  if(tp == 0.0)     // ToDo: Figure out if we need a tolerance! Like:  tp <= tol * tau1
  {
    rsError("Zero attack not yet implemented.");
    // ToDo: Handle the special case for when tp (time of the peak) is zero, i.e. the user requests
    // a zero attack time.

  }

  T a1 = 1/tau1;
  T c  = a1 * tp;
  T k  = findDecayScalerLess1(c);
  T a2 = k*a1;
  T hp = exp(-a1*tp) - exp(-a2*tp); // peak height

  *tau2   = 1/a2;
  *scaler = 1/hp;

  // ToDo: Make this a static member function of rsModalFilterWithAttack and document it. I think, 
  // it computes tau2 (the decay time constant for the second filter) from tau1 (the decay time 
  // constant of the first filter) and tp (the desired time instant of the peak). It also computes
  // the overall scaler for the whole signal so that the peak amplitude becomes 1.
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

  if(attackTime == TPar(0))
  {
    modalFilter1.setModalParameters(frequency, amplitude, decayTime, startPhase, sampleRate);
    modalFilter2.setModalParameters(frequency, TPar(0),   TPar(0),   startPhase, sampleRate);
  }
  else
  {
    TPar tau1, tau2, scaler;
    tau1 = decayTime;
    expDiffScalerAndTau2(tau1, attackTime, &tau2, &scaler);
    amplitude *= scaler;
    modalFilter1.setModalParameters(frequency, amplitude, tau1, startPhase, sampleRate);
    modalFilter2.setModalParameters(frequency*detuneFactor, amplitude, tau2, startPhase, sampleRate);
  }

  // ToDo: 
  // 
  // - "detuneFactor" does not really work well because when attack and decay are very similar, 
  //   the amplitude explodes (due to a high value of "scaler") - either remove this parameter or 
  //   find a way to alleviate this (maybe the amplitude excess can be computed - math has to be 
  //   worked out) the detune is supposed to introduce some roughness into the transient by 
  //   detuning the quickly decaying sinusoid
  //
  // - Try to get rid of the special case branch for attackTime == 0. Make sure that the general
  //   formulas that are used in the lower branch do the right thing in this case. We may need to 
  //   take care of divisions by zero. Make sure that they behave corrently in the limiting case.
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
  /*
  sampleRate         = 44100.0;
  referenceFrequency = 440.0;
  referenceAttack    = 0.1;
  referenceDecay     = 1.0;
  numModes           = maxNumModes;
  */
  // ToDo: Init these in the class header! ...done


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
void rsModalFilterBank<TSig, TPar>::setMaxNumModes(int newMax)
{
  maxNumModes = newMax;
  numModes = rsMin(numModes, maxNumModes);   // Truncate current numModes, if needed

  frequencies.resize(maxNumModes);
  amplitudes.resize(maxNumModes);
  attackTimes.resize(maxNumModes);
  decayTimes.resize(maxNumModes);
  startPhases.resize(maxNumModes);

  modalFilters.resize(maxNumModes);

  // ToDo: Maybe init the arrays to all zeros.
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setNumModes(int newNum)
{
  rsAssert(newNum <= maxNumModes, "Requested number of modes exceeds maximum");
  numModes = rsMin(newNum, maxNumModes);
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setReferenceFrequency(TPar newFrequency)
{
  referenceFrequency = newFrequency;
  calculateModalFilterCoefficients();
}

template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::setReferenceAmplitude(TPar newAmplitude)
{
  referenceAmplitude = newAmplitude;
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
void rsModalFilterBank<TSig, TPar>::setModeParams(int m, 
  TPar freq, TPar amp, TPar attack, TPar decay, TPar phase)
{
  rsAssert(m >= 0 && m < numModes);  // Maybe factor out into isValidModeIndex()

  frequencies[m] = freq;
  amplitudes[m]  = amp;
  attackTimes[m] = attack;
  decayTimes[m]  = decay;
  startPhases[m] = phase;

  updateFilterCoeffs(m);             // Maybe make that update optional
}


// DEPRECATE THIS:
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
void rsModalFilterBank<TSig, TPar>::updateFilterCoeffs(int m)
{
  modalFilters[m].setModalParameters(  
    referenceFrequency * frequencies[m], 
    referenceAmplitude * amplitudes[m],
    referenceAttack    * attackTimes[m],
    referenceDecay     * decayTimes[m], 
    startPhases[m], 
    sampleRate); 
}

// Maybe rename to updateModalFilterCoeffs()
template<class TSig, class TPar>
void rsModalFilterBank<TSig, TPar>::calculateModalFilterCoefficients()
{
  int nm = rsMin(numModes, (int)frequencies.size(), (int)amplitudes.size(), (int)decayTimes.size());
  nm = rsMin(nm, (int)startPhases.size());
  for(int m = 0; m < nm; m++)
    updateFilterCoeffs(m);

  // Old:
  //size_t nm = rsMin((size_t)numModes, frequencies.size(), amplitudes.size(), decayTimes.size());
  //nm = rsMin(nm, startPhases.size());
  //for(size_t m = 0; m < nm; m++)
  //  updateFilterCoeffs(m);
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

template<class TSig, class TPar>
bool rsModalFilterBank<TSig, TPar>::checkClassInvariants() const
{
  bool ok = true;

  //ok &= numModes    == getNumModes();  
  // Will later produce a stack overflow when we call checkClassInvariants there. The code is not 
  // active yet, though.

  ok &= numModes    <= maxNumModes;

  ok &= maxNumModes == frequencies.size();
  ok &= maxNumModes == amplitudes.size();
  ok &= maxNumModes == attackTimes.size();
  ok &= maxNumModes == decayTimes.size();
  ok &= maxNumModes == startPhases.size();

  ok &= maxNumModes == modalFilters.size();

  return ok;

  // ToDo: 
  //
  // - Maybe we should really switch to a arrays-of-struct design rather than the current 
  //   struct-of-arrays design. Then we could remove the check that all the arrays have the same 
  //   size
  //
  // - Verify that the vector sizes really correspond to the maxNumModes. I'm not sure anymore if
  //   I have intended it this way, i.e. using the vector sizes for maxNumModes and numModes to 
  //   only use a part of these vectors or if I wanted to use the sizes for the current number and
  //   the capacities for the max. The former makes more sense thoug, so that should be the way to
  //   do it.
  //
  // - Check that the attack times are less than the decay times with some safety margin.
}



//=================================================================================================
/*

Ideas

-Don't demand modes to start at zero amplitude - allow envelopes whose start-amplitude is nonzero.
 startAmp must be <= peakAmp though ...or oes it? Maybe if it's not, we can do a 2-stage decay? We 
 would need to define what the peakTime is then supposed to mean - it should still be the instant
 at which peakAmp is reached, although the term "peak" is then wrong. Being able to start from 
 nonzero seems to refelct some data seen in the real-world better and should not impact any 
 per-sample calculations

-Make a class rsKeyVelMap<T> that stores datapoints of type T for an arbitrary number of midi key 
 and velocities an lets the user retrieve the data for an arbitrary key/vel pair using (bilinear) 
 interpolation using the data, i.e. find the 4 closest available key/vel datapoints and form a
 weighted average. The type T is intended to be, for example, a specification of ModalFilterBank 
 data - but it should work for arbitrary types. The intention is to simulate multisamples but using
 synthesis data instead of samples. Maybe make it more general like rsBilinearMap or 
 rsBivariateData<TIdx, TData> with different possibilities for the index (e.g. float or int). 
 Should work similar to rsNodeBasedFunction but with node-locations in 2D.

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

- Additional ideas for per-mode parameters:
  - Tremolo freq and amount. Can be implemented via beating between two slightly detuned modal 
    filters or directly as amplitude modulation via sine-producing filter in the LFO-range (which
    then somehow needs to amplitude normalize the output to 1 - maybe a complex phasor based 
    implementation can do that)
  - Panning - or maybe better: LFO-modulated panning, i.e. stereo-tremolo. Maybe a cheaper version
    would be to (dynamically) pan even and odd modes such that not every mode needs it dedicated 
    Pan-LFO. But of course, having a dedicated pans LFO for each mode with its own frequency will
    create a much more complex modulation.
  - Delay. Could be implemented by a simple delay line

-Maybe it could even make sense to have a brickwall filter feature in rsModalFilterBank similar to 
 what we have in the wavetable oscillator class. So we could have functions like setMinModeIndex(), 
 setMaxModeIndex() or maybe setHighpass(), setLowpass(). Maybe the function names should be made 
 consistent with those in the wavetable osc.

Modeling transients:
-Transients are modeled as superposition attack/decay envelope filters (i.e. zero frequency) with
 delay
-Maybe an interative matching procedure can be used to find the parameters of the filters
-Or maybe use modeling in terms of a general pole/zero model
-Maybe model broadband transients by allpass filters

A general instrument based on modal synthesis could look like:
- output = modes + transient  where 
    modes = modal-bank -> dispersion allpass -> equalizer
    transient = transient allpass (bank?) -> equalizer
  so we would feed an exciter input into a parallel connection of a modal bank (with post 
  processing) and a mdoule responsibel for creating transients. Maybe for sustained excitation 
  signals like noise, they should not go into the transient module - or maybe with an envelope 
  applied. We could use a differential envelope detector to extract a transient from an input noise
  generator.

- Create a variation that provides a a two-stage decay envelope. Instead of a single decay 
  pameteter, we would have two: EarlyDecay and LateDecay and a DecayMix parameter. The way this 
  would be implemented is by creating a weighted sum of two decaying modes with two different 
  decay times. The computation of the decay-time for the subtracted attack mode would have to use
  a more complex formula. To achieve this, we would now need 3 modal filters for a single mode:
  Two with the two different decay times and one for the (subtracted) attack portion.

- Provide a delay parameter for each mode by prepending or appending a delay line.

- Allow mode beating by using two slightly detuned versiond of the modal filter. The parameters 
  should be: BeatFreq, BeatAmount. See the experiments about sine beating for how to realize it.
  We want it to look like amplitude modulation. Maybe we could also have some stereo-beating
  effect, i.e. different phases in the envelope in left and right signal. 

- Maybe the carrier wave itself could also have a different phase in left and right signal.

- Each mode should also have a Pan parameter. 

- Maybe it could also have a waveshape parameter. We could pass the waveshape through a 
  waveshaper, like tanh. It may also be interesting to pass groups of modes through waveshapers 
  to let them interact nonlinearly. Maybe the octaves of the fundamental together with the octaves
  of mode 3, i.e. sines at f0 and 3*f0, would produce some nice powerchord like sound?

- Provide a different implementation that splits the computation into two multiplicative parts
  used for the sinusoid and the attack/decay envelope. The envelope may operate on the absolute
  value of the inputs signal. Such a splitting may lead to more efficient algorithms compared to
  producing differently enveloped sinuosoids dierectly and adding them.

- Full list of possible per mode parameters: Freq, Amp, Pan, PhaseLeft, PhaseRight, EarlyDecay, 
  LateDecay, DecayMix, Attack, Delay, BeatFreq, BeatAmount, StereoBeat, FreqByAmp,
  The FreqByAmp feature could be implemented in a way similar to rsModalFilterNonLinear but 
  instead of computing the instantaneous amplitude using sqrt, we may run a filter with the same 
  envelope settings and zero freq in parallel. Essentially, we bump the phase a bit further by 
  additional multiplication of the state with another rotor which depends on the instantaneous 
  amplitude..I think. Need to think about it more - it's tricky to do efficiently but it would be
  a really nice fetaure. It would even nicer if we could generally modulate the frequency somehow.
  With an envelope and an LFO.

*/