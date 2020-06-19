
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

    //rsArrayTools::difference(w, N);
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

      //rsArrayTools::difference(w, N); // does not take wrapping into account
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
    wi += w[n];
    y[n] = a[n] * sin(wi); }
}

template<class T>
void rsSingleSineModeler<T>::synthesizeFromAmpFreqAndPhaseMod(
  const T* a, const T* w, const T* pm, int N, T* y)
{
  T wi = T(0); // integrated w
  for(int n = 0; n < N; n++) {
    wi += w[n];
    y[n] = a[n] * sin(wi + pm[n]); }
}
// maybe we should base everything on cosine for consistency with the rsSinusoidalModel - but maybe
// we should use the sine there


//-------------------------------------------------------------------------------------------------
// internal sub-algorithms:

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaForward(T y0, T yR, T w, T* a, T* p)
{
  const T margin = 1.e-12;  
  // ad hoc - do tests, what is best - should be different for float and double maybe some 
  // multiple of the epsilon? or maybe a power? maybe PI * pow(eps, 1.5) or something?

  if( rsDistanceToMultipleOf(w, PI) <= margin )
  {
    *a = rsAbs(y0);
    if(y0 > 0)        *p = +PI/2;
    else if(y0 < 0)   *p = -PI/2;
    else              *p =  0;
    return;
  }




  T s, c, sR;
  rsSinCos(w, &s, &c);

  T x = yR-y0*c;
  T y = y0*s;

  *p = atan2(y, x);  
  // computes atan(y/x) - so we should make sure that x != 0 - can this occur? ..yes - when
  // yR == y0*cos(w) - but this is no problem - atan2 handles this case itself

  //*p = atan2(y0*s, yR-y0*c);


  s  = sin(*p);
  sR = sin(*p + w);
  if( rsAbs(s) > rsAbs(sR) )
    *a = y0 / s;
  else {
    if(sR == 0.0)
    {  
      *a = 0.0;  
      // may not be a good idea - we don't want to enforce zero samples
      // needs a tolerance - i think, we already ruled that out with the checks above - at least
      // one of the sines must be > 0 because w is guranteed to be not too close to a half-period
      // of the sine
    }
    else
      *a = yR / sR; }
  // what if s and sR are both close to zero?, what if w == 0, but y0 != yR
}
// copy documentation from rsSineAmplitudeAndPhase

template<class T>
void rsSingleSineModeler<T>::phaseAndAmpFormulaBackward(T y0, T yL, T w, T* a, T* p)
{
  T sw = sin(w);
  T cw = cos(w);

  *p = atan2(-y0*sw, yL-y0*cw) + PI;  // avoid the addition of PI by rotating arg to atan2
  *a = y0 / sin(*p);

  int dummy = 0;
  // ...not yet finished...needs a switch to avoid div-by-zero
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

  //T power = 1.0; // experimental - doesn't seem to help


  // that's the old, imprecise version
  if(precision <= 1)
  {
    T* y = a;
    for(int n = 0; n < N; n++)
      y[n] = rsAbs(x[n]);         // todo: apply shadower here (shadows are casted only rightward)
    connectPeaks(y, N, a);  
    // todo: pass false, i precision = 0 - this should indicate to not use parabolic interpolation
    return;
  }

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
void rsSingleSineModeler<T>::sigAndAmpToPhase(const T* x, const T* a, int N, T* p)
{
  for(int n = 0; n < N; n++)
  {
    //p[n] = asin(x[n] / a[n]);  // this may produce NaN!

    rsAssert(a[n] >= T(0), "amplitudes must be nonegative");  // maybe try to lift that restriction
    rsAssert(a[n] >= rsAbs(x[n]), "amplitudes can't be less than signal values");

    if(a[n] == 0)
      p[n] = 0;
    else
      p[n] = asin(x[n] / a[n]); 
  }
  unreflectPhase(x, p, N);
}

template<class T>
void rsSingleSineModeler<T>::phaseToFreq(const T* p, int N, T* w)
{
  for(int n = 0; n < N; n++)
    w[n] = rsWrapToInterval(p[n], 0.0, 2*PI); // is this needed? try without!
  rsArrayTools::unwrap(w, N, 2*PI); // look at code comment there - optimize!
  rsArrayTools::difference(w, N);

  // todo: do the unwrapping on the fly
}

template<class T>
void rsSingleSineModeler<T>::sigAndFreqToPhaseAndAmp(const T* x, const T* w, int N, T* p, T* a)
{
  for(int n = 0; n < N-1; n++)
    phaseAndAmpFormulaForward(x[n], x[n+1], w[n], &a[n], &p[n]);
  phaseAndAmpFormulaBackward(x[N-1], x[N-2], w[N-1], &a[N-1], &p[N-1]);
}

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
void rsSingleSineModeler<T>::connectPeaks(const T* y, int N, T* a)
{
  // Make these function parameters - these determine whether we use the height of a parabolic
  // interpolant for the peak instead of the array value itself. We may also use the actual 
  // position/time of the peak in the linear interpolation loop, but that might not always be
  // a good idea, even when the interpolant is used for the height:
  bool parabolicHeight = true; 
  bool parabolicTime   = true;  // makes sense only, if parabolicHeight is also true
  // todo: introduce a peak-precision parameter: 0 - take peak sample directly, 1: parabolic
  // 2: quartic, 3: sixtic, etc. - see code for zero-crossing finder for reference
  // rename parabolicTime to exactTime - but then it may not be possible to use it in place anymore
  // maybe make a second function for higher order peak finding

  int nL = 0,     nR;            // index of current left and right peak
  T   tL = T(nL), tR;            // position or time of current left and right peak
  T   yL = y[0],  yR;            // amplitude or height of current left and right peak
  for(int n = 1; n < N-1; n++) {
    if(y[n] >= y[n-1] && y[n] >= y[n+1]) {             // there's a peak or plateau at y[n]...
      nR = n; tR = T(nR); yR = y[n];
      if(parabolicHeight) {
        using Poly = rsPolynomial<T>;
        T c[3]; Poly::fitQuadratic_m1_0_1(c, &y[n-1]);  // c = polynomial coeffs of parabola
        if(c[2] != 0) {                                 // TODO: use a tolerance
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
// actually, this shoudl also dsitinguish between xt, xi (x used for the test and x used for the
// interpolation)

template<class T>
void rsSingleSineModeler<T>::exactPeakPositionAndHeight(
  const T* x, int N, int n0, int precision, T* pos, T* height)
{
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
  Poly::interpolant(a, T(-p), T(1), &x[n0-p], degree+1); // +1 bcs it takes number of datapoints, allocates
  Poly::derivative(a, ad, degree);
  dt = rsPolynomial<T>::rootNear(dt, ad, degree-1, T(-1), T(1));
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
void rsSingleSineModeler<T>::unreflectPhase(const T* x, T* p, int N)
{
  for(int n = 1; n < N; n++) {
    if(x[n] >= 0) {
      if(x[n] < x[n-1])  p[n] =  PI - p[n];  }   // x is positive and going down -> zone 2
    else {
      if(x[n] < x[n-1])  p[n] = -PI - p[n]; }}   // x is negative and going down -> zone 3

  // Post-process - compenaste too late transitions:
  refinePhase(p, N);     // maybe rename 
  //refinePhase(p, N);  // doesn't make a difference
}
// zone 1: 0...pi/2, zone 2: pi/2...pi, zone 3: -pi/...-pi/2, zone 4: -pi/2...0
// can too early transitions also happen? i've not yet seen one


/*

ToDo:

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