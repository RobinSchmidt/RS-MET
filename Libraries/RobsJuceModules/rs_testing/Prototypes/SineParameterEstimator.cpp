
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

  bool parabolicHeight = true; 
  bool parabolicTime   = false; // makes sense only, if parabolicHeight is also true

  T *y = a;                     // re-use a for temporary storage
  for(int n = 0; n < N; n++)
    y[n] = rsAbs(x[n]);         // todo: apply shadower here (shadows are casted only rightward)


  // factor out into connectPeaks(y, N, a)
  int nL = 0,     nR;            // index of current left and right peak
  T   tL = T(nL), tR;            // position or time of current left and right peak
  T   yL = y[0],  yR;            // amplitude or height of current left and right peak
  for(int n = 1; n < N-1; n++)
  {
    if(y[n] >= y[n-1] && y[n] >= y[n+1])  // there's a peak or plateau at y[n]...
    {
      nR = n; tR = T(nR); yR = y[n];

      if(parabolicHeight) 
      {
        using Poly = rsPolynomial<T>;
        T c[3]; Poly::fitQuadratic_m1_0_1(c, &y[n-1]);  // c = polynomial coeffs of parabola
        if(c[2] != 0)  // TODO: use a tolerance
        {
          T dt = Poly::quadraticExtremumPosition(c); // time offset of peak between -1..+1
          yR   = Poly::evaluate(dt, c, 2);           // height of peak
          if(parabolicTime)                          // we may or may not use the time offset..
            tR += dt;                                // ..in the linear interpolation loop below
        }
        // quadraticExtremumPosition computes c[1]/c[2], so the tolerance should be based on the 
        // ratio |c[1]| and |c[2]| - if abs(c[2]) < (small * c1), skip the step
      }

      for(int i = nL; i < nR; i++)                 // lerp peaks from nL to nR
        a[i] = rsLinToLin(T(i), tL, tR, yL, yR);   // ToDo: optimize!
      nL = nR; tL = tR; yL = yR;                   // update for next iteration
    }
  }

  nR = N-1;
  tR = T(nR);
  yR = y[nR];
  for(int i = nL; i < nR; i++)
    a[i] = rsLinToLin(T(i), tL, tR, yL, yR); // optimize!

}