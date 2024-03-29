
//-------------------------------------------------------------------------------------------------
// Analysis:

template<class T>
void rsSingleSineModeler<T>::analyzeAmpAndPhase(const T* x, int N, T* a, T* p) const
{
  if(algo == Algorithm::ampViaPeaks)  {
    sigToAmpsViaPeaks(x, N, a, ampEnvPrecision);
    sigAndAmpToPhase(x, a, N, p);
    return;  }

  if(algo == Algorithm::freqViaFormula || algo == Algorithm::freqViaZeros) {
    sigToFreq(x, N, p);
    sigAndFreqToPhaseAndAmp(x, p, N, p, a);
    return;  }

  rsError("Unknown algorithm");
}

template<class T>
void rsSingleSineModeler<T>::analyzeAmpAndFreq(const T* x, int N, T* a, T* w) const
{
  if(algo == Algorithm::ampViaPeaks) {
    analyzeAmpAndPhase(x, N, a, w);         // w temporarily used for phase
    phaseToFreq(w, N, w);    
    return; }

  if(algo == Algorithm::freqViaFormula || algo == Algorithm::freqViaZeros) {
    sigToFreq(x, N, w);
    sigAndFreqToPhaseAndAmp(x, w, N, w, a);
    phaseToFreq(w, N, w);
    // sigAndFreqToAmp(x, w, N, a); instead of sigAndFreqToPhaseAndAmp -> difference does not 
    // work. Hmm - it seems, we can not use the omegas from the freq-estimation pass. We need 
    // indeed compute the phases from the originally etsimated omegas and difference them - why?
    return;  }

  rsError("Unknown algorithm");
}

template<class T>
void rsSingleSineModeler<T>::analyzeAmpFreqAndPhaseMod(const T* x, int N, T* a, T* w, T* pm) const
{
  if(algo == Algorithm::ampViaPeaks) {
    analyzeAmpAndPhase(x, N, a, pm);        // pm (phase-mod) temporarily used for phase itself
    phaseToFreq(pm, N, w);
    smoothFreqs(w, N, freqMedianOrder, freqAverageOrder);
    phaseAndFreqToPhaseMod(pm, w, N, pm);   // convert phase to phase-mod
    // Maybe we should make a decision here whether or not to call this function - if it's not 
    // called (because smoothing is off), we can just fill the pm-array with zeros.
    // Note that calling analyzeAmpAndFreq here instead of the 1st two lines won't work because 
    // phaseAndFreqToPhaseMod needs phase and freq as input - the 2 lines look very similar to the 
    // ones in analyzeAmpAndFreq, but the arguments to the called functions are different.
    return;  }

  if(algo == Algorithm::freqViaFormula || algo == Algorithm::freqViaZeros){
    sigToFreq(x, N, w);
    if(freqMedianOrder > 0 && freqAverageOrder > 0) {
      sigAndFreqToPhaseAndAmp(x, w, N, pm, a);
      smoothFreqs(w, N, freqMedianOrder, freqAverageOrder);
      phaseAndFreqToPhaseMod(pm, w, N, pm);  }
    else {
      sigAndFreqToPhaseAndAmp(x, w, N, w, a);
      phaseToFreq(w, N, w);
      rsArrayTools::fillWithZeros(pm, N);    }
    return; }

  rsError("Unknown algorithm");
}

template<class T>
void rsSingleSineModeler<T>::sigToFreq(const T* x, int N, T* w) const
{
  if(algo == Algorithm::freqViaFormula)
    sigToFreqViaFormula(x, N, w);
  else
    sigToFreqViaZeros(x, N, w);
}


//-------------------------------------------------------------------------------------------------
// Synthesis:

template<class T>
void rsSingleSineModeler<T>::synthesizeFromAmpAndPhase(const T* a, const T* p, int N, T* y)
{
  for(int n = 0; n < N; n++)
    y[n] = a[n] * sin(p[n]);
}

template<class T>
void rsSingleSineModeler<T>::synthesizeFromAmpAndFreq(const T* a, const T* w, int N, T* y)
{
  T wi = T(0); // integrated w
  for(int n = 0; n < N; n++) {
    //wi += w[n];    // maybe wrap-around at 2*pi, maybe use wi = fmod(wi + w[n], 2*PI)
    wi = fmod(wi + w[n], 2*PI);
    y[n] = a[n] * sin(wi); }
}

template<class T>
void rsSingleSineModeler<T>::synthesizeFromAmpFreqAndPhaseMod(
  const T* a, const T* w, const T* pm, int N, T* y)
{
  T wi = T(0); // integrated w
  for(int n = 0; n < N; n++) {
    wi += w[n];
    wi = fmod(wi + w[n], 2*PI);  // this makes the unit test fail
    //y[n] = a[n] * sin(wi + pm[n]);   // why is this commented out?
  }
}
// compare numerical error with and without fmod - use longer arrays and maybe a high-freq sinewave
// (like w = 0.99*PI) that drives the unwrapped phase up quickly

template<class T>
void rsSingleSineModeler<T>::synthesizeFromAmpAndFreqNew(const T* a, const T* w, int N, T* y, T p0)
{
  T wi = T(0);     // integrated w
  y[0] = a[0] * sin(p0);
  for(int n = 1; n < N; n++) {
    wi = fmod(wi + T(0.5)*(w[n-1]+w[n]), T(2*PI));
    y[n] = a[n] * sin(wi + p0); }
}


// maybe we should base everything on cosine for consistency with the rsSinusoidalModel - but maybe
// we should use the sine there




//-------------------------------------------------------------------------------------------------
// internal sub-algorithms:


template<class T>
bool rsSingleSineModeler<T>::handlePhaseAmpEdgeCase(T y0, T w, T* a, T* p)
{
  const T margin = 1.e-12;  
  if( rsDistanceToMultipleOf(w, PI) <= margin ) {
    if(y0 > 0)        *p = +PI/2;
    else if(y0 < 0)   *p = -PI/2;
    else              *p =  0;
    *a = rsAbs(y0);
    return true; }
  else
    return false;
  // margin chosen ad hoc - do tests, what is best - should be different for float and double maybe
  // some multiple of the epsilon? or maybe a power? maybe PI * pow(eps, 1.5) or something?
}

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaForward(T y0, T yR, T w, T* a, T* p)
{
  if( handlePhaseAmpEdgeCase(y0, w, a, p) )
    return;
  T s, c, sR;
  rsSinCos(w, &s, &c);
  *p = atan2(y0*s, yR-y0*c);
  s  = sin(*p);
  sR = sin(*p + w);
  rsAssert(s != 0 || sR != 0); // should be guaranteed by our edge case handling above
  if( rsAbs(s) > rsAbs(sR) )
    *a = y0 / s;
  else
    *a = yR / sR;
}
// what if s and sR are both close to zero?, what if w == 0, but y0 != yR
// copy documentation from rsSineAmplitudeAndPhase

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaBackward(T y0, T yL, T w, T* a, T* p)
{
  if( handlePhaseAmpEdgeCase(y0, w, a, p) )
    return;
  T s, c, sL;
  rsSinCos(w, &s, &c);
  *p = atan2(y0*s, y0*c-yL);
  s  = sin(*p);
  sL = sin(*p - w);
  rsAssert(s != 0 || sL != 0);
  if( rsAbs(s) > rsAbs(sL) )
    *a = y0 / s;
  else
    *a = yL / sL;
}

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaCentral1(T yL, T y0, T yR, T w, T* a, T* p)
{
  if( handlePhaseAmpEdgeCase(y0, w, a, p) )
    return;

  T s, c, sL, sR;
  rsSinCos(w, &s, &c);
  T pB = atan2(y0*s, y0*c-yL);  // phase from backward estimation
  T pF = atan2(y0*s, yR-y0*c);  // phase from forward estimation

  *p   = rsInterpolateWrapped(pB, pF, 0.5, T(-PI), T(+PI));
  // todo: is the some way to estimate the accuracy of pB and pF and give more weight to the more 
  // accurate one? Currently, we just use an unweighted average, i.e. 0.5*(pB+pF) - maybe try to 
  // predict yL, yR from y0, w and a and give more weight to the one for which the prediction 
  // matches the actual value better...or maybe we can jointly estimate p,a by minimizing a squared
  // error:
  //   err = (yL - a*sin(p-w))^2 + (yR - a*sin(p+w))^2 = min

  // compute sines at sample indices n-1, n, n+1:
  sL = sin(*p - w);
  s  = sin(*p    );
  sR = sin(*p + w);

  // compute amplitudes based on the 3 samples:
  T aL = yL / sL;
  T a0 = y0 / s;
  T aR = yR / sR;

  // compute weights for the 3 amplitudes for how much they should contribute to the final 
  // amplitude - they idea is that those values, for which we a had a division by (close to) zero
  // in the previous step should contribute less:
  T wL = T(1) / rsAbs(aL);
  T w0 = T(1) / rsAbs(a0);
  T wR = T(1) / rsAbs(aR);
  T wS = wL + w0 + wR;    // sum of weights

  // compute final amplitude as weighted average over the 3 estimates:
  *a = (wL*aL + w0*a0 + wR*aR) / wS;  // what if wS == 0?

  // this algo is quite ad-hoc - i'm not sure, if it's good - tests are needed. Maybe try the 
  // joint estimation of p and a by minimization of the error defined above and compare results. 
  // They should give same results in the case of an fixed-freq sine but with a sine-sweep, they
  // may give different results, so we may use sweep to assess the quality of both approaches
  // if the (linear) sweep does not reveal any quality difference, maybe try a quadratic sweep 
  // and/or introduce and amplitude fade, too and measure which of the formulas estimates the
  // actual instantaneous phase and amp better - maybe try it with using the correct instantaneous
  // freq and with its estimated value via the freq-formula
}

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaCentral2(T yL, T y0, T yR, T w, T* a, T* p)
{
  if(handlePhaseAmpEdgeCase(y0, w, a, p))
    return;
  T sw, cw; rsSinCos(w, &sw, &cw);  // sw = sin(w), cw = cos(w)
  *p   = atan2(cw*(yL+yR) - 2*y0*(cw*cw-sw*sw), sw*(yR-yL));
  T sp = sin(*p);
  *a   = y0 / sp;
  // todo: handle sp == 0 and if *a < 0, make it positive and invert the phase
  // maybe if p is too close to zero, use the other formula from above?


  // This formula was derived (see SineParameters.txt) from:
  //   yL = a * sin(p-w)
  //   y0 = a * sin(p)
  //   yR = a * sin(p+w)
  // solving the middle equation for a:
  //   a = y0 / sin(p)
  // and combining the yL,yR equations into the minimization problem:
  //   E = (yL - a*sin(p-w))^2 + (yR - a*sin(p+w))^2 = min
  // taking the partial derivative of the error function E with respect to p and setting it to 
  // zero. This gives the atan2-based formula for p - but we had to swap/rotate the argument 
  // compared to the formula derived there -> figure out, why
}

template<class T>
void rsSingleSineModeler<T>::sigToFreqViaFormula(const T* x, int N, T* w)
{
  // The algorithm uses rsSineFrequency as its core to estimate the frequency at each sample. 
  // However, it was observed, that this function gives unreliable results, whenever there's a 
  // zero-crossing, so we first compute the (expected) reliabilities of the measurements at each 
  // sample and then actually use a weigthed sum of the estimates at the sample and at its two 
  // neighbours, where there weights are determined by the reliabilities. This way, the unreliable
  // etsimates at the zero-crossings will be more or less replaced by a weighted average of the
  // estimates at neighbour samples. The reliability itself is computed as ratio of a sample's 
  // absolute value to the average of the absolute values of its neighbours. We use the w array to
  // temporarily store the reliabilities and the overwrite it with the actual frequency estimates 
  // in a second pass.

  rsAssert(x != w);
  rsAssert(N >= 3);
  T small1 = 1.e-8;  // ad-hoc - do tests what value is best
  T small2 = 1.e-8;  // dito - best choice could be the same as small1 but maybe not
  int n;

  // compute reliabilities:
  T* r = w;                         // re-use w-array for temporary storage of reliabilities
  for(n = 1; n < N-1; n++) {
    T num = rsAbs(x[n]);
    T den = T(0.5) * (rsAbs(x[n-1]) + rsAbs(x[n+1]));
    if(den >= small1*num)
      r[n] = num / den;
    else
      r[n] = T(0); }
  r[0] = 0; r[N-1] = 0;
  // factor out into omegaFormulaReliabilities(*x, N, *r)

  // compute radian frequencies:
  T rL = r[0], rC = r[1], rR, rS;
  T wL = T(0), wC = freqFormula(x[0], x[1], x[2]), wR;
  for(n = 1; n < N-2; n++) {
    rR = r[n+1];                              // reliability of right neighbour sample
    if(rR > T(0))
      wR = freqFormula(x[n], x[n+1], x[n+2]); // frequency of right neighbour sample
    else
      wR = T(0);
    rS = rL + rC + rR;                     // reliability sum of all 3 - used as normalizer
    if( rS > small2 )                      // sum is large enough as denominator
      w[n] = (rL*wL + rC*wC + rR*wR) / rS; //   -> use weighted sum of neighbour estimates
    else                                   // all 3 reliabilities are too small
      w[n] = w[n-1];                       //   -> repeat previous value
    rL = rC; rC = rR; wL = wC; wC = wR; }  // updates for next iteration

  rS = rL + rC;                     // handle n = N-2 ouside the loop
  if( rS > small2 ) 
    w[n] = (rL*wL + rC*wC) / rS;
  else
    w[n] = w[n-1];
  w[0] = w[1]; w[N-1] = w[N-2];     // handle edges at n = 0 and n = N-1

  // factor out into omegaFormulaOmegas(*x, *r, N, *w)
  // r == w is allowed, is also x == w allowed? not so important but would be nice
}

template<class T>
void rsSingleSineModeler<T>::sigToFreqViaZeros(const T* x, int N, T* w)
{
  for(int n = 0; n < N; n++)
    w[n] = rsSineFrequencyAt(x, N, n, false);
}

template<class T>
void rsSingleSineModeler<T>::sigToAmpsViaPeaks(const T* x, int N, T* a, int precision)
{
  // todo: take a shadowing-time parameter and use a peak-shadower

  // Algo:
  // -obtain shadowed version of abs(x)
  // -find peaks
  //  -maybe refine their values by using the maximum through a parabola
  // -connect them by linear interpolation


  // old version, doesn't allocate:
  if(precision <= 1)
  {
    T* y = a;
    for(int n = 0; n < N; n++)
      y[n] = rsAbs(x[n]);         // todo: apply shadower here (shadows are casted only rightward)
    connectPeaks(y, N, a, precision == 1);  
    return;
  }

  // new version, allowing for higher precision:
  using Vec = std::vector<T>;
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = rsAbs(x[n]);        // use shadower
  connectPeaks(&y[0], x, N, a, precision); 
  // todo: pass y as xTest and x as xInterpolate

  //rsError("High preicison not yet implemented");
}
// get rid of the duplication

template<class T>
void rsSingleSineModeler<T>::sigAndAmpToPhase(const T* x, const T* a, int N, T* p) const
{
  for(int n = 0; n < N; n++) {
    rsAssert(a[n] >= T(0), "amplitudes must be nonegative");  // maybe try to lift that restriction
    rsAssert(a[n] >= rsAbs(x[n]), "amplitudes can't be less than signal values");
    if(a[n] == 0)
      p[n] = 0;
    else
      p[n] = asin(x[n] / a[n]); }

  // dispatch between different unreflection algorithms:
  if(phaseAlgo == PhaseUnreflectAlgorithm::fromSignalSlope)
    unreflectPhaseFromSig(x, p, N);
  else
    unreflectPhaseFromSigAndAmp(x, a, p, N); // this dispatches further
}

/*
// old implementation - todo: compare numerical errors with new
template<class T>
void rsSingleSineModeler<T>::phaseToFreq(const T* p, int N, T* w)
{
  for(int n = 0; n < N; n++)
    w[n] = rsWrapToInterval(p[n], 0.0, 2*PI); // without it, the unit test fails
  rsArrayTools::unwrap(w, N, 2*PI); // look at code comment there - optimize!
  rsArrayTools::difference(w, N);

  //rsPlotArrays(N, p, w);
}
*/

template<class T>
void rsSingleSineModeler<T>::phaseToFreq(const T* p, int N, T* w)
{
  // This code basically does something like first unwrapping the phase array and then taking the
  // difference of the unwrapped array but with the unwrapping done on the fly. I think, this 
  // should avoid some numerical errors since we avoid that the unwrapped phase grows too large.
  T pOld = 0;
  for(int n = 0; n < N; n++) {
    T pNew = p[n];
    int k  = rsUnwrapFactor(pNew, pOld, 2*PI, 0); // maybe we should pass the k from the previous iteration as last argument?
    pNew  += 2*k*PI;
    w[n]   = pNew - pOld;
    pOld   = pNew; }

  //rsPlotArrays(N, p, w);
}

// old:
template<class T>
void rsSingleSineModeler<T>::sigAndFreqToPhaseAndAmp(const T* x, const T* w, int N, T* p, T* a)
{
  for(int n = 0; n < N-1; n++)
    phaseAndAmpFormulaForward(x[n], x[n+1], w[n], &a[n], &p[n]);
  phaseAndAmpFormulaBackward(x[N-1], x[N-2], w[N-1], &a[N-1], &p[N-1]);
}
// todo: use forward formula only for 1st sample (index 0) and for(n = 1; n < N-1; n++) use the 
// central formula

/*
// new - makes indentity-resynthesis tests fail - why?:
template<class T>
void rsSingleSineModeler<T>::sigAndFreqToPhaseAndAmp(const T* x, const T* w, int N, T* p, T* a)
{
  phaseAndAmpFormulaForward(x[0], x[1], w[0], &a[0], &p[0]);
  for(int n = 1; n < N-1; n++)
    phaseAndAmpFormulaCentral1(x[n-1], x[n], x[n+1], w[n], &a[n], &p[n]);
  phaseAndAmpFormulaBackward(x[N-1], x[N-2], w[N-1], &a[N-1], &p[N-1]);
}
// in the context of this function, it seems to be a bad idea that w[0] represents the initial 
// phase - this will also mess up a[0] and p[0]
// the failure of the unit test in not just about numerical precision - the output is totally wrong
// ..maybe if we use this function for analysis, we have to use a different synthesis formula, too?
// ..i think, it would be better anyway, to use trapezoidal integration of the frequency to obtain
// the phase and represent the start-phase as extra-value when doing synthesis from freqs only 
// (without additional phase-modulation). That makes transfomations of the data more convenient and
// also does not require the phase/amp data to accomodate for the special meaning of the first 
// frequency datapoint
*/

template<class T>
void rsSingleSineModeler<T>::sigAndFreqToAmp(const T* x, const T* w, int N, T* a)
{
  T dummy; // for the phase output argument of the called function which we are not interested in
  for(int n = 0; n < N-1; n++)
    phaseAndAmpFormulaForward(x[n], x[n+1], w[n], &a[n], &dummy);
  phaseAndAmpFormulaBackward(x[N-1], x[N-2], w[N-1], &a[N-1], &dummy);

  // what about the last value? use the backward formula - maybe use forward formula only for 
  // sample 0 and a central formula for 1...N-2
}

/*
template<class T>
void rsSingleSineModeler<T>::freqToPhase(const T* w, int N, T* p, bool wrap)
{
  rsArrayTools::cumulativeSum(w, N, p);  // is this correct?
}
*/
// maybe make a loop that wraps on the fly during accumulation - might be numerically better than
// wrapping after cumulating everything because we avoid the intermediate results to grow large

template<class T>
void rsSingleSineModeler<T>::phaseAndFreqToPhaseMod(const T* p, const T* w, int N, T* pm)
{
  T wi = T(0);
  for(int n = 0; n < N; n++) {
    wi += w[n];
    pm[n] = p[n]-wi; }

  for(int n = 0; n < N; n++)
    pm[n] = rsWrapToInterval(pm[n], -PI, PI); 
  // todo: make wrapping optional and/or maybe allow for returning unwrapped phase-mod data (the 
  // raw result from above is neither wrapped nor unwrapped, i think) - unwrapped will be more 
  // useful anyway, as we may want to filter the pm-values before resynthesis
}


template<class T>
void rsSingleSineModeler<T>::smoothFreqs(T* w, int N, int medianOrder, int averageOrder)
{
  for(int i = 0; i < medianOrder;  i++) 
    rsArrayTools::movingMedian3pt( w, N, w);
  for(int i = 0; i < averageOrder; i++) 
    rsArrayTools::movingAverage3pt(w, N, w);
}

template<class T>
inline void lerpPeaks(const T* y, int nL, int nR, T tL, T tR, T yL, T yR, T* a)
{
  for(int i = nL; i < nR; i++)
  {
    T ai = rsLinToLin(T(i), tL, tR, yL, yR); // ToDo: optimize!
    a[i] = rsMax(ai, rsAbs(y[i]));                  // we want no stick-outs!
  }
}
// lerp peaks from nL to nR (not including nR)
// maybe make member

template<class T>
void rsSingleSineModeler<T>::connectPeaks(const T* y, int N, T* a, bool useParabola)
{
  // Make these function parameters - these determine whether we use the height of a parabolic
  // interpolant for the peak instead of the array value itself. We may also use the actual 
  // position/time of the peak in the linear interpolation loop, but that might not always be
  // a good idea, even when the interpolant is used for the height:
  bool parabolicHeight = useParabola; 
  bool parabolicTime   = useParabola;  // makes sense only, if parabolicHeight is also true
  // todo: introduce a peak-precision parameter: 0 - take peak sample directly, 1: parabolic
  // 2: quartic, 3: sixtic, etc. - see code for zero-crossing finder for reference
  // rename parabolicTime to exactTime - but then it may not be possible to use it in place anymore
  // maybe make a second function for higher order peak finding

  T smalll = 1.e-8;              // ad hoc - tests needed, what's best
  int nL = 0,     nR;            // index of current left and right peak
  T   tL = T(nL), tR;            // position or time of current left and right peak
  T   yL = y[0],  yR;            // amplitude or height of current left and right peak
  for(int n = 1; n < N-1; n++) {
    if(y[n] >= y[n-1] && y[n] >= y[n+1]) {             // there's a peak or plateau at y[n]...
      nR = n; tR = T(nR); yR = y[n];
      if(parabolicHeight) {
        using Poly = rsPolynomial<T>;
        T c[3]; Poly::fitQuadratic_m1_0_1(c, &y[n-1]);  // c = polynomial coeffs of parabola
        //if(c[2] != 0) {                                 // TODO: use a tolerance
        if(rsAbs(c[2]) >= rsAbs(smalll*c[1])) {         // avoid div-by-close-to-zero
          T dt = Poly::quadraticExtremumPosition(c);    // time offset of peak between -1..+1
          yR   = Poly::evaluate(dt, c, 2);              // height of peak
          if(parabolicTime)                             // we may or may not use the time offset..
            tR += dt;       }}                          // ..in the linear interpolation loop below
      lerpPeaks(y, nL, nR, tL, tR, yL, yR, a);
      nL = nR; tL = tR; yL = yR; }}                     // update for next iteration
  nR = N-1; tR = T(nR); yR = y[nR];          // handle last peak to last sample
  lerpPeaks(y, nL, nR, tL, tR, yL, yR, a);
}
// ToDo: factor out the parabolic refinement and allow it to use higher order polynomials - it's
// really important to estimate the amplitude accurately - otherwise we'll get jaggies in the
// instantaneous frequency measurements
// quadraticExtremumPosition computes c[1]/c[2], so the tolerance should be based on the 
// ratio |c[1]| and |c[2]| - if abs(c[2]) < (small * c1), skip the step
// actually, this should also distinguish between xt, xi (x used for the test and x used for the
// interpolation)
// -to handle the left and right edge, it might be better to use linear extrapolation from existing 
//  peaks rather than just using the first and last sample - that would work better for artificial, 
//  pure sine tones - but maybe for more natural sounds, the current strategy is actually the 
//  better way? using linear extrapolation would complicate the code

template<class T>
void rsSingleSineModeler<T>::exactPeakPositionAndHeight(
  const T* x, int N, int n0, int precision, T* pos, T* height)
{
  rsAssert(n0 >= 0 && n0 < N, "n0 is out of range");

  static const int maxPrecision = 4;
  rsAssert(precision <= maxPrecision); // for higher precisions, we need to allocate a larger a-array below

  int p = rsMin(precision, n0, (N-1)-n0);       // todo: verify this formula
  if(p == 0 || n0 == 0 || n0 == N-1) {
    *pos    = T(n0);
    *height = x[n0];
    return; }

  // First, fit a parabola and find its maximum:
  using Poly = rsPolynomial<T>;
  T a[ 2*maxPrecision+1];                       // polynomial coeffs of fitted polynomial
  T ad[2*maxPrecision];                         // ...and its derivative
  Poly::fitQuadratic_m1_0_1(a, &x[n0-1]);
  T dt = Poly::quadraticExtremumPosition(a);    // delta t - time offset of peak from n0
  if(p == 1 || n0 == 1 || n0 == N-2) {
    *pos    = T(n0) + dt;
    *height = Poly::evaluate(dt, a, 2);
    return; }

  // Next, fit a polynomial of degree 2*precision to some number of samples near the peak and find 
  // its maximum using Newton iteration, using the location of the peak of the parabola as initial
  // guess:
  int degree = 2*p;
  Poly::interpolant(a, T(-p), T(1), &x[n0-p], degree+1); 
  // degree+1 because it is number of datapoints, function allocates - maybe avoid this by using ad 
  // as workspace for the function (then we need to allocate 2*maxPrecision+1 for the ad array also 
  // (i think) - make a unit test that tests all possible precisions)


  Poly::derivative(a, ad, degree);
  dt = Poly::rootNear(dt, ad, degree-1, T(-1), T(1));
  *pos    = T(n0) + dt;
  *height = Poly::evaluate(dt, a, degree);
}
// maybe move to somewhere in the Analysis section - maybe together with the algo to find 
// zero-crossings into a class rsFeatureFinder - maybe it could also find locations of other 
// features such as center-of gravity, etc. - or rsTimeDomainFeatureFinder - the peaks could also
// be used to increase the time-resolution of pitch-estimation to quarter cycles - we would use:
// locations of upward-zero -> peak -> downward-zero -> trough for each cycle (makes sense only
// for bandpass signals, i.e. signals that have sinusoidal shape)

template<class T>
void rsSingleSineModeler<T>::connectPeaks(
  const T* xt, const T* xi, int N, T* env, int precision)
{
  rsAssert(xi != env && xt != env);  // this does not work in place

  int nL = 0,            nR;    // index of current left and right peak
  T   tL = T(nL),        tR;    // position or time of current left and right peak
  T   xL = rsAbs(xi[0]), xR; 
  for(int n = 1; n < N-1; n++)
  {
    if(xt[n] >= xt[n-1] && xt[n] >= xt[n+1])
    {
      nR = n; 
      exactPeakPositionAndHeight(xi, N, nR, precision, &tR, &xR); // allocates
      xR = rsAbs(xR);
      lerpPeaks(xi, nL, nR, tL, tR, xL, xR, env);
      nL = nR; tL = tR; xL = xR;
    }
  }

  nR = N-1; 
  tR = T(nR); 
  xR = rsAbs(xi[nR]);       // handle last peak to last sample
  lerpPeaks(xi, nL, nR, tL, tR, xL, xR, env);
}
// for using higher order interpolation, it's actually a bad idea to use an array of absolute 
// values - when the intepolant-width exceeds the width of the sinusoids lobe, we should actually
// use the original signal itself - but then we may have to look for maxima *and* minima

template<class T>
T refinePhase(T p, T pL, T pR, int n) // n is only passed for debugging
{
  // Idea: compare phase p to the average of its two neighbours pL, pR and also compare and 
  // appropriately reflected phase to this average. Return either original or reflected value,
  // depending on their distances to this average - we want the smaller distance.

  if( pL > pR )  
    pL -= 2*PI;                        // avoids undesirable adjustments around wrap-arounds

  T pi2 = 0.5*PI;
  T pa  = T(0.5)*(pL+pR);
  T pm  = p;

  //if(p < pi2 && pR >= pi2)             // transition from zone 1 to zone 2
  if(pL < pi2 && pR >= pi2)             // transition from zone 1 to zone 2
    pm =  PI - p;
  //else if(p < -pi2 && pR >= -pi2)      // transition from zone 3 to zone 4
  else if(pL < -pi2 && pR >= -pi2)      // transition from zone 3 to zone 4
    pm = -PI - p;

  if(rsAbs(pa-pm) < rsAbs(pa-p))
    return pm;
  else
    return p;
}
// make member
// do we need more branches? what about the reverse transitions from zone 2 to 1 and from 4 to 3?

template<class T>
void refinePhase(T* p, int N)
{
  T pi2 = 0.5*PI;
  for(int n = 1; n < N-1; n++)
    p[n] = refinePhase(p[n], p[n-1], p[n+1], n);
}
// make member
// could it make sense to run this function multiple times until the phase array converges? maybe
// try this with white-noise inputs - the bandpass noise is too well-behaved for that

template<class T>
void rsSingleSineModeler<T>::unreflectPhaseFromSig(const T* x, T* p, int N)
{
  for(int n = 1; n < N; n++) {
    if(x[n] >= 0) {
      if(x[n] < x[n-1])  p[n] =  PI - p[n];  }   // x is positive and going down -> zone 2
    else {
      if(x[n] < x[n-1])  p[n] = -PI - p[n]; }}   // x is negative and going down -> zone 3

  // Post-process - compenaste too late transitions:
  refinePhase(p, N);     // maybe rename 
}
// zone 1: 0...pi/2, zone 2: pi/2...pi, zone 3: -pi/...-pi/2, zone 4: -pi/2...0
// can too early transitions also happen? i've not yet seen one

//-------------------------------------------------------------------------------------------------
// under construction:

template<class T>
void rsSingleSineModeler<T>::unreflectPhaseFromSigAndAmp(const T* x, const T* a, T* p, int N) const
{
  // estimate frequencies from preliminary phases:
  std::vector<double> w(N);
  phaseToFreq(p, N, &w[0]);  
  for(int n = 0; n < N; n++)
    w[n] = rsAbs(w[n]);        // enforce positive freqs

  // The so obtained w-array has errors whereever p has maxima or minima (possibly with a 1-sample 
  // delay), so we mark all w-values at the extrema of p as ivalid (we use the impossible value of 
  // -1 as mark):
  for(int n = 1; n < N-1; n++) {
    if(rsArrayTools::isPeakOrValley(p, n)) {
      T wa = (w[n-1] + w[n] + w[n+1]) / 3;     // 3-point average
      T dL = rsAbs(wa - w[n-1]);
      T dC = rsAbs(wa - w[n]);
      T dR = rsAbs(wa - w[n+1]);
      if(dC > dL && dC > dR)  w[n]   = -1;
      if(dL > dC && dL > dR)  w[n-1] = -1;     // this seems to never happen -> check this!
      if(dR > dC && dR > dL)  w[n+1] = -1; }}

  // Replace the marked invalid w-values by the average of their neighbour values:
  for(int n = 1; n < N-1; n++) {
    if(w[n] == -1)
      w[n] = 0.5*(w[n-1] + w[n+1]); }

  // phaseToFreq results in (unwrapped version of) w[n] = p[n] - p[n-1], i.e. a backward estimate, 
  // but we want a central estimate, so we do a 2-point forward moving average to fix that:
  for(int n = 0; n < N-1; n++)
    w[n] = 0.5*(w[n] + w[n+1]);
  // ...todo: verify, if that's the correct way to do it

  // maybe factor out everything up to here into a function reflectedPhasesToFreqs


  // (secondary) dispatch (primary was done by caller):
  if(phaseAlgo == PhaseUnreflectAlgorithm::fromFreq)
    unreflectPhaseFromFreq(&w[0], p, N);
  else
    unreflectPhaseFromSigAmpAndFreq(x, a, &w[0], p, N, true);
}
// todo: maybe use a workspace for w - avoid allocation here - maybe we can completely avoid memory
// allocation (but it's not tertibly important to do so - this is offline code anyway - it would 
// just be nice)

// input: a phase p assumed to be in the range -pi/2...+pi/2 as returned, for example, by asin
// output: a possible other phase q for which sin(q) == sin(p) due to reflection symmetry of the 
//         sine (this has nothing to do with periodic shifting)
template<class T>
T getReflectedSinePhase(T p)
{
  if(p >= 0) return  PI - p;
  else       return -PI - p;
}
// make member, maybe rename to getAlternativeSinePhase...maybe without the get

template<class T>
void rsSingleSineModeler<T>::unreflectPhaseFromFreq(const T* w, T* p, int N)
{
  for(int n = 1; n < N; n++)
  {
    T wa = 0.5*(w[n-1] + w[n]);     // average freq between sample n-1 and sample n
    T pp = p[n-1] + wa;             // predicted phase at sample n from phase at sample n-1
    if( pp >= PI ) pp -= 2*PI;      // wrap around pp, if necessarry

    // compute the possible alternative phase at sample n resulting from reflection:
    T pr = getReflectedSinePhase(p[n]);

    // if it gives a closer match between predicted and measured phase, reflect p[n]:
    if( rsPhaseDistance(pp, pr) < rsPhaseDistance(pp, p[n]) ) 
      p[n] = pr;
    // else leave p[n] as is
  }
}

template<class T>
void rsSingleSineModeler<T>::unreflectPhaseFromSigAmpAndFreq(
  const T* x, const T* a, const T* w, T* p, int N, bool central)
{
  // Idea: for each sample compare predicted signal from original and alternative version of phase
  // to actual signal value and choose the one, for which the difference is smaller.


  int n = 1;     // start at p[1], p[0] is the anchor and never changed - should stay close to 0
  if(central) {
    // use both (forward/backward) prediction errors for all samples except the last, if desired:
    while(n < N-1) {
      T wL = 0.5*(w[n-1] + w[n]);   // average freq between sample n-1 and sample n
      T wR = 0.5*(w[n] + w[n+1]);   // average freq between sample n and sample n+1

      // left and right amplitudes used in the (backward and forward) prediction - i'm not sure, if 
      // it's a good idea to use an average here - maybe using aL = a[n-1], aR = a[n+1] could be 
      // better - we'll see:
      T aL = 0.5*(a[n-1] + a[n]);
      T aR = 0.5*(a[n] + a[n+1]);

      // compute prediction error eo using original p[n]:
      T eL = x[n-1] - aL * sin(p[n] - wL);  // left (backward) prediction error
      T eR = x[n+1] - aR * sin(p[n] + wR);  // right (forward) prediction error
      T eo = rsAbs(eL) + rsAbs(eR);         // total prediction error using original p[n]

      // compute prediction error ea using alternative to p[n]:
      T pa = getReflectedSinePhase(p[n]);
      eL = x[n-1] - aL * sin(pa - wL);      // left predication error
      eR = x[n+1] - aR * sin(pa + wR);      // right prediction error
      T ea = rsAbs(eL) + rsAbs(eR);         // total prediction error using alternative phase

      // replace p[n] with its alternative, if this reduces the prediction error:
      if(ea < eo)
        p[n] = pa;
      n++;
    }
  }

  // use backward prediction error for all remaining samples (which are either all samples except 
  // the first (index 0) or only the last sample, i.e. n == 1 or n == N-1):
  while(n < N)
  {
    T wL = 0.5*(w[n-1] + w[n]);
    T aL = 0.5*(a[n-1] + a[n]);
    T pa = getReflectedSinePhase(p[n]);
    T eo = x[n-1] - aL * sin(p[n] - wL);
    T ea = x[n-1] - aL * sin(pa   - wL);
    if(rsAbs(ea) < rsAbs(eo))
      p[n] = pa;
    n++;
  }
}



// pNew: preliminary phase at sample n, pOld: phase at sample n-1,w: frequency between sample 
// n-1 and n, rteurn: adjusted phase at sample n
template<class T>
T adjustPhase(T pNew, T pOld,  T w)
{
  T pTgt = pOld + w;  // target phase at sample n assuming freq w between n-1 and n
  if(pTgt > PI)
    pTgt -= 2*PI;   // wrap around

  if(w < PI/2)
  {
    // w < PI/2 - we may have transitions: 1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1

    if(pNew > 0 && (pTgt > PI/2 || pOld > PI/2 ))
      pNew =  PI - pNew;                                 // transition from zone 1 to zone 2

    if(pNew < 0 && (pTgt < -PI/2 || pOld > PI/2 ))
      pNew = -PI - pNew;                                 // transition from zone 3 to zone 4


  }
  else 
  {
    // w >= PI/2 - we may have transitions: 1 -> 3, 2 -> 4, 3 -> 1, 4 -> 2

    rsError("Not yet implemented");

  }

  if(pTgt > PI)    // hwy test pTgt and not pNew?
    pNew -= 2*PI;  
    
  return pNew;
}


template<class T>
void rsSingleSineModeler<T>::unreflectPhase2(const T* w, T* p, int N)
{
  // not yet tested
  // under construction - uses linear extrapolation from current phase and omega approach as 
  // oppsoed to using the input signal values to determine the zone

  for(int n = 1; n < N; n++)
    p[n] = adjustPhase(p[n], p[n-1], 0.5*(w[n-1]+w[n]));
  return;


  double pi = PI;

  for(int n = 1; n < N; n++)
  {
    T po = p[n-1];                   // old phase at sample n-1
    T pa = p[n];                     // actual, measured phase at sample n

    //T pt = p[n-1] + w[n];
    T pt = p[n-1] + 0.5*(w[n-1]+w[n]);
    // predicted, target phase at sample n ..should we use w[n-1] or w[n]? - maybe best would be 
    // 0.5*(w[n-1]+w[n]) because that represents the *average* measured frequency between previous 
    // and current sample - the integration algo that we actually use for synthesis is not really
    // relevant here because the input w-array is the measured freq and not the differenced phase 
    // that we would use in synthesis (btw. maybe we could use trapzoidal integration in synthesis, 
    // too - we just have to take care, how to represent the initial phase in this case)

 

    if(pt > PI)
      pt -= 2*PI;  // test


    //if(p[n] > 0 && pt > PI/2)
    if(p[n] > 0 && (pt > PI/2 || p[n-1] > PI/2 ))
      p[n] =  PI - p[n];            // transition from zone 1 to zone 2
    //else if(p[n] < 0 && pt < -PI/2)
    else if(p[n] < 0 && (pt < -PI/2 || p[n-1] > PI/2 ))
      p[n] = -PI - p[n];            // transition from zone 3 to zone 4
    // or should we use p[n] <= 0 

    if(pt > PI)                     // wrap seems to happen 1 sample too early - will be repaired below
    //if(p[n] > PI)           
      p[n] -= 2*PI;                 // transition for zone 2 to zone 3 (wrap-around)
  }

  //for(int n = 1; n < N; n++) {      // repair the too early wrap-around
  //  if(p[n] < -PI)
  //    p[n] += PI; }


  //for(int n = 1; n < N; n++) {      // repair the too early wrap-around
  //  if(p[n] < -PI)
  //    p[n] += 2*PI; }

  //for(int n = 1; n < N; n++) {      // repair the too early wrap-around
  //  if(p[n] < -2*PI)
  //    p[n] += 2*PI; }

}
// maybe try the same approach but without using the w-array by estimating w from p[n] and p[n-1]
// as w = p[n] - p[n-1]...with some care about wrapping...or maybe we need p[n-1] - p[n-2] because
// p[n] is actually the value, we want to compute...

// this seems to work only, if w < pi/2



template<class T>
void rsSingleSineModeler<T>::unreflectPhase3(const T* w, T* p, int N)
{
  // initially (i.e. at n = 0), the zone is either 1 or 4:
  int zone;
  if(p[0] >= 0) zone = 1;
  else          zone = 4;

  // pre-process to wrap phases < 0 around:
  for(int n = 0; n < N; n++)
    if(p[n] < 0)
      p[n] += 2*PI;

  // actual unreflection, using phases in the range 0 <= p < 2*pi:
  for(int n = 1; n < N; n++)
  {
    zone = newPhaseZone(p[n], p[n-1], 0.5*(w[n-1]+w[n]), zone);
    if(zone == 2) p[n] =  PI - p[n];
    if(zone == 3) p[n] = -PI - p[n] + 2*PI;
  }

  // post-process to wrap phases >= pi around:
  for(int n = 0; n < N; n++)
    if(p[n] >= PI)
      p[n] -= 2*PI;
}
// zones for the phase p:
// 1: 0      <= p < pi/2
// 2: pi/2   <= p < pi
// 3: pi     <= p < 3*pi/2    or  -pi   <= p < -pi/2
// 4: 3*pi/2 <= p < 2*pi      or  -pi/2 <= p < 0

// numerically, it would be better to avoid the pre/post-process and the +2*PI in the zone==3 
// branch...but for development, it's more convenient that way - maybe try to get rid of it later

// maybe remove oldZone parameter, maybe add xNew, xOld, amp instead - they may be useful to figure 
// out if phase has moved forward or backward between samples n-1 and n
template<class T>
int rsSingleSineModeler<T>::newPhaseZone(T pNew, T pOld, T w, int oldZone)
{
  rsAssert(w >= 0 && w <= PI); 
  // the analysis should produce values in that range, right? ...nope: only the zero-crossing based
  // freq-estimation algo is guaranteed to produce nonnegative values the (acos-based) formula may
  // totally produce negative frequencies

  T pTgt = pOld + w;  // target phase
  if(pTgt >= 2*PI)
    pTgt -= 2*PI;

  if(pNew >= 0) {                     // new zone = 1 or 2
    if(pTgt < PI/2)   return 1;
    else              return 2; }
  else {                              // new zone = 3 or 4
    if(pTgt < 3*PI/2) return 3;
    else              return 4; }
  // oldZone is not even used here - maybe remove the parameter - we may not need it - in fact, the
  // pOld variable (and therefore also pTgt together with w) fully contains that information, so 
  // it's redundant

  // but to accomodate negative phase-modulation, maybe we should decide first, if the phase was 
  // more likely to have made a forward or backward step by comparing, how close xNew is to 
  // a*sin(pOld+w) and a*sin(pOld-w) - so we need also to pass xNew and amplitude a[n] to this 
  // function - if this works can be tested with phase-modulated signals


  rsError("Case not handled");
  return 0; // preliminary
}
// input: 
// pNew: preliminary, measured phase at sample n (range: -pi/2 <= pNew < pi/2) 
//       ..or -pi/2 < pNew <= pi/2?
// pOld: already adjusted phase at sample n-1 (range: 0 <= pOld < 2*pi)
// w: radian frequency between sample n-1 and sample n
// oldZone: zone of pOld (1,2,3,4)
//
// output: the new zone, i.e. the zone of pNew



/*

ToDo:
-check for edge cases in freqFormula (see unit test - i think, there may be some nan-results for
 certain inputs), phaseAndAmpFormulaBackward - implement and use phaseAndAmpFormulaCentral
-Figure out and document to which sample instant the instantaneous freq w[n] actually applies. 
 Maybe it applies to sample n-1 or n or n+1 or n-0.5 or n+0.5? Or maybe it represents the *average* 
 instantaneous freq for the time interval between two sample instants? I think, it's probably the
 average inst freq between n and n+1 but i'm not sure...may also be the interval between n-1 and n.


Other ideas for phase unreflection:
-minimize the sum of the distances to left and right neighbour (i think, this may be equivalent to
 minimzing the distance to their midpoint, as we do now)
-minimize the distance to the phase predicted by linearly extrapolating from two left neighbours
 this may actually not need the x-data - it can take a raw phase-array containing reflected
 phases and turn it into an unreflected one
-maybe an algorithm for this could also take the amplitude array as input - i don't know, if tha 
 information could be useful for unreflection

How about estimating amplitude and frequency first and then computing the phase from two successive
samples:

  x[n]   = x0 = a * sin(p)
  x[n+1] = xR = a * sin(p+w)

can this be solved for p?

x0 / sin(p) = xR / sin(p+w)

*/