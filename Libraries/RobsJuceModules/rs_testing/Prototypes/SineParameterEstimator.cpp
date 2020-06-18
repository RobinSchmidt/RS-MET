
template<class T>
void rsSineParameterEstimator<T>::sigToOmegasViaFormula(const T* x, int N, T* w)
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
  T wL = T(0), wC = omegaFormula(x[0], x[1], x[2]), wR;
  for(n = 1; n < N-2; n++) {
    rR = r[n+1];                               // reliability of right neighbour sample
    if(rR > T(0))
      wR = omegaFormula(x[n], x[n+1], x[n+2]); // frequency of right neighbour sample
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
void rsSineParameterEstimator<T>::sigToAmpsViaPeaks(const T* x, int N, T* a, int precision)
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

template<class T>
void rsSineParameterEstimator<T>::sigAndAmpToPhase(const T* x, const T* a, int N, T* p)
{
  for(int n = 0; n < N; n++)
    p[n] = asin(x[n] / a[n]);
  unreflectPhase(x, p, N);
}

template<class T>
void rsSineParameterEstimator<T>::phaseToFreq(const T* p, int N, T* w)
{
  rsAssert(p != w);                 // hmm - i guess, it would actually work in place - try it!
  for(int n = 0; n < N; n++)
    w[n] = rsWrapToInterval(p[n], 0.0, 2*PI); // is this needed? try without!
  rsArrayTools::unwrap(w, N, 2*PI); // look at code comment there - optimize!
  rsArrayTools::difference(w, N);
}

template<class T>
void rsSineParameterEstimator<T>::phaseAndFreqToPhaseMod(const T* p, const T* w, int N, T* pm)
{
  T wi = w[0];
  for(int n = 1; n < N; n++)
  {
    wi += w[n];
    pm[n] = rsWrapToInterval(p[n]-wi, -PI, PI); 
    // make optional and/or maybe allow for returning unwrapped phase-mod - this will be more
    // useful anyway, as we may want to filter the pm-values before resynthesis
  }
}

template<class T>
void rsSineParameterEstimator<T>::synthesizeFromAmpFreqPhaseMod(
  const T* a, const T* w, const T* pm, int N, T* y)
{
  T wi = w[0]; // integrated w
  for(int n = 1; n < N; n++)
  {
    wi += w[n];
    y[n] = a[n] * sin(wi + pm[n]);
  }
}


//-------------------------------------------------------------------------------------------------
// internal sub-algorithms:

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
void rsSineParameterEstimator<T>::connectPeaks(const T* y, int N, T* a)
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
void rsSineParameterEstimator<T>::exactPeakPositionAndHeight(
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
void rsSineParameterEstimator<T>::connectPeaks(
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
void rsSineParameterEstimator<T>::unreflectPhase(const T* x, T* p, int N)
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
-in connectPeaks, we need a higher order polynomial fit to estimate the amplitude more accurately
 to get rid of the frequency jaggies

Other ideas for phase unreflection:
-minimize the sum of the distances to left and right neighbour (i think, this may be equivalent to
 minimzing the distance to their midpoint, as we do now)
-minimize the distance to the phase predicted by linearly extrapolating from two left neighbours
-maybe an algorithm for this could also take the amplitude array as input - i don't know, if tha 
 information could be useful for unreflection

Maybe this class should also have the synthesis functions - that may mean, we need another name -
maybe rsSineRecreator or rsSineRepresenter something - it represents *any* signal as a 
time-varying sinewave of the form:

  x[n] = a[n] * sin( p[n] )

with intantaneous amplitude a[n] and instantaneous phase p[n], where p[n] can be either given 
directly or it can be further split into instantaneous frequency w[n] (omega) an and instantaneous 
phase-modulation pm[n] term like:

  p[n] = sum_{k=0}^{n-1} w[k]  +  pm[n]   ...check, if the sum should run to n-1 or to n

..the goal should be that the analysis functions compute either a[n],p[n] or a[n],w[n],pm[n] such
that the corresponding synthesis functions exactly reconstruct any signal - regardless whether it
is actually a sinusoid or soemthing completely different. 

*/