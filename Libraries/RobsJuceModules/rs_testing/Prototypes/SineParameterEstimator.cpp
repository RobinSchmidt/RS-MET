
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
void rsSineParameterEstimator<T>::sigToAmpsViaPeaks(const T* x, int N, T* a)
{
  // todo: take a shadowing-time parameter and use a peak-shadower

  // Algo:
  // -obtain shadowed version of abs(x)
  // -find peaks
  //  -maybe refine their values by using the maximum through a parabola
  // -connect them by linear interpolation

  T *y = a;                     // re-use a for temporary storage
  for(int n = 0; n < N; n++)
    y[n] = rsAbs(x[n]);         // todo: apply shadower here (shadows are casted only rightward)
  connectPeaks(y, N, a);
}

template<class T>
void rsSineParameterEstimator<T>::sigAndAmpToPhase(const T* x, const T* a, int N, T* p)
{
  for(int n = 0; n < N; n++)
    p[n] = asin(x[n] / a[n]);
  unreflectPhase(x, p, N);
}

//-------------------------------------------------------------------------------------------------
// internal sub-algorithms:

template<class T>
inline void lerpPeaks(const T* y, int nL, int nR, T tL, T tR, T yL, T yR, T* a)
{
  for(int i = nL; i < nR; i++)                 // lerp peaks from nL to nR
  {
    T ai = rsLinToLin(T(i), tL, tR, yL, yR);
    a[i] = rsMax(ai, y[i]);  // we want no stick-outs!
    //a[i] = ai;
  }
}
// make member or maybe move elsewhere

template<class T>
void rsSineParameterEstimator<T>::connectPeaks(const T* y, int N, T* a)
{
  // Make these function parameters - these determine whether we use the height of a parabolic
  // interpolant for the peak instead of the array value itself. We may also use the actual 
  // position/time of the peak in the linear interpolation loop, but that might not always be
  // a good idea, even when the interpolant is used for the height:
  bool parabolicHeight = true; 
  bool parabolicTime   = true;  // makes sense only, if parabolicHeight is also true

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

      //for(int i = nL; i < nR; i++)                 // lerp peaks from nL to nR
      //  a[i] = rsLinToLin(T(i), tL, tR, yL, yR);   // ToDo: optimize!


      nL = nR; tL = tR; yL = yR; }}                // update for next iteration

  nR = N-1; tR = T(nR); yR = y[nR];          // lerp from last peak to last sample
  lerpPeaks(y, nL, nR, tL, tR, yL, yR, a);
  //for(int i = nL; i < nR; i++)
  //  a[i] = rsLinToLin(T(i), tL, tR, yL, yR); // optimize!
}
// quadraticExtremumPosition computes c[1]/c[2], so the tolerance should be based on the 
// ratio |c[1]| and |c[2]| - if abs(c[2]) < (small * c1), skip the step
// maybe the lerp-loop can be factored out - maybe it should be added to 
// rsArrayTools::lerp(T *y, int iL, int iR, T min, T max) ..actually, it may be useful, if it
// writes the lerped value into a[i] only if it is >= y[i], otherwise, write y[i] - this will
// ensure, that nothing can stick out - when we do this, we actually can use parabolicTime without
// risk

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

  if(p < pi2 && pR >= pi2)             // transition from zone 1 to zone 2
    pm =  PI - p;
  else if(p < -pi2 && pR >= -pi2)      // transition from zone 3 to zone 4
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
}
// zone 1: 0...pi/2, zone 2: pi/2...pi, zone 3: -pi/...-pi/2, zone 4: -pi/2...0
// can too early transitions also happen? i've not yet seen one


/*

Other ideas for phase unreflection:
-minimize the sum of the distances to left and right neighbour (i think, this may be equivalent to
 minimzing the distance to their midpoint, as we do now)
-minimize the distance to the phase predicted by linearly extrapolating from two left neighbours

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