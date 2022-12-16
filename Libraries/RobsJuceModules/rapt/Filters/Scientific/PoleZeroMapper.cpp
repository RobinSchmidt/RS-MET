template<class T>
void rsPoleZeroMapper<T>::sLowpassToLowshelf(Complex *z, Complex *p, T *k, Complex *zNew, 
  Complex *pNew, T *kNew, int N, T G0, T G)
{
  if( G0 == 0.0 )
  {
    rsArrayTools::copy(z, zNew, N);
    rsArrayTools::copy(p, pNew, N);
    *kNew = *k * G;
    return;
  }

  T kTmp;
  Complex *zTmp = new Complex[2*N];
  Complex *pTmp = new Complex[2*N];
  T *a  = new T[N+1];
  T *b  = new T[N+1];
  T *bS = new T[2*N+1];  // shelving magnitude squared numerator coeffs

  // we design a boost filter in any case - a dip filter will later be obtained by swapping 
  // poles and zeros:
  bool dip = false;
  if( G < G0 )
  {
    dip = true;
    G   = T(1) / G;
    G0  = T(1) / G0;
  }

  // scale poles/zeros/gain of the lowpass such that the resulting low-shelving filter has a gain 
  // of sqrt(G*G0) at unit frequency:
  //T GB = sqrt(G0*G);                       // bandwidth gain
  //T GC = sqrt((GB*GB-G0*G0)/(G*G-G0*G0));  // desired lowpass gain at unit frequency
  //PrototypeDesigner::scaleToMatchGainAtUnity(z, p, k, zTmp, pTmp, &kTmp, N, GC);  // commented for test

  // for debug:
  rsArrayTools::copy(z, zTmp, N);
  rsArrayTools::copy(p, pTmp, N);
  kTmp = *k;

  // obtain magnitude-squared numerator polynomial for shelving filter:
  rsPolynomial<T>::complexRootsToRealCoeffs(zTmp, b, N);
  rsPolynomial<T>::complexRootsToRealCoeffs(pTmp, a, N);
  rsPrototypeDesigner<T>::shelvingMagSqrNumeratorFromLowpassTransfer(b, a, kTmp, N, G0, G, bS);

  // obtain zeros, poles and gain of the new shelving filter:
  //findPolynomialRoots(bS, 2*N, &zTmp[-1]);
  //int numLeftZeros = onlyLeftHalfPlane(zTmp, zTmp, 2*N);
  rsPrototypeDesigner<T>::getLeftHalfPlaneRoots(bS, zTmp, 2*N);  // test - replaces the 2 commented lines above - seems to work

  // if we make a dip-filter, poles and zeros exchange roles:
  if( dip == false )
  {
    rsArrayTools::copy(zTmp, zNew, N);
    rsArrayTools::copy(pTmp, pNew, N);
    *kNew = sqrt(fabs(bS[2*N]));
  }
  else
  {
    rsArrayTools::copy(zTmp, pNew, N);
    rsArrayTools::copy(pTmp, zNew, N);
    *kNew = T(1) / sqrt(fabs(bS[2*N]));
  }

  delete[] zTmp;
  delete[] pTmp;
  delete[] a;
  delete[] b;
  delete[] bS;
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToLowpass(Complex* z, Complex* p, T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N, T wc)
{
  for(int n = 0; n < N; n++) zNew[n] = wc * z[n];
  for(int n = 0; n < N; n++) pNew[n] = wc * p[n];
  int nz = rsGetNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew = *k * pow(wc, N-nz);
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToHighpass(Complex *r, Complex *rNew, int N, T wc)
{
  for(int n = 0; n < N; n++)
  {
    if( rsIsInfinite(r[n]) )
      rNew[n] = Complex(0.0, 0.0);
    else
      rNew[n] = wc / r[n];
  }
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToHighpass(Complex* z, Complex* p,T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N, T wc)
{
  sLowpassToHighpass(z, zNew, N, wc);
  sLowpassToHighpass(p, pNew, N, wc);
  Complex n = rsProductOfFiniteFactors(z, N);
  Complex d = rsProductOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).real();

  // sign-factor: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p);
  // ..or maybe not? ..consult octave/matlab implementation
  int nz = rsGetNumFiniteValues(z, N);
  if( rsIsOdd(N+nz) )
    *kNew = - *kNew;
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToBandpass(Complex* r, Complex* rNew, int N, T wc, T bw)
{
  Complex tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1      = T(0.5) * bw * r[n];
    tmp2      = -bw * r[n];

    // this is common with the bandreject case - factor out:
    tmp2      = T(0.25) * tmp2*tmp2;
    tmp2      = sqrt(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToBandpass(Complex* z, Complex* p, T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N, T wl, T wu)
{
  T wc = sqrt(wl*wu);  // center (radian) frequency
  T bw = wu-wl;        // bandwidth
  sLowpassToBandpass(z, zNew, N, wc, bw);
  sLowpassToBandpass(p, pNew, N, wc, bw);
  for(int n = 0; n < N; n++)
  {
    if( rsIsInfinite(z[n]) )
    {
      zNew[n]   = Complex(0.0,       0.0);
      zNew[N+n] = Complex(RS_INF(T), 0.0);
    }
  }
  int nz = rsGetNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew  = *k * pow(bw, N-nz);
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToBandreject(Complex* r, Complex* rNew, int N, T wc, T bw)
{
  Complex tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1      = T(0.5) * bw / r[n];
    tmp2      = -bw / r[n];

    // this is common with the bandpass case - factor out:
    tmp2      = T(0.25) * tmp2*tmp2;
    tmp2      = sqrt(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}

template<class T>
void rsPoleZeroMapper<T>::sLowpassToBandreject(Complex* z, Complex* p, T* k, Complex* zNew, 
  Complex* pNew, T* kNew, int N, T wl, T wu)
{
  T wc = sqrt(wl*wu);  // center (radian) frequency
  T bw = wu-wl;        // bandwidth
  sLowpassToBandreject(z, zNew, N, wc, bw);
  sLowpassToBandreject(p, pNew, N, wc, bw);
  Complex sz = Complex(0.0, wc);
  for(int n = 0; n < N; n++)
  {
    if( rsIsInfinite(z[n]) )
    {
      zNew[n]   =  sz;
      zNew[N+n] = -sz;
    }
  }
  Complex n = rsProductOfFiniteFactors(z, N);
  Complex d = rsProductOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).real();

  // sign-factor - see m-file: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p);
  int nz = rsGetNumFiniteValues(z, N);
  if( rsIsOdd(N+nz) )
    *kNew = - *kNew;
}

template<class T>
void rsPoleZeroMapper<T>::prototypeToAnalogLowpass(Complex *poles, int numPoles, Complex *zeros,
  int numZeros, T* gain, T targetCutoff)
{
  int i;
  for(i = 0; i < numPoles; i++)
    poles[i] = poles[i] * targetCutoff;
  for(i = 0; i < numZeros; i++)
    zeros[i] = zeros[i] * targetCutoff;
  *gain *= 1; // preliminary - \todo adjust gain: *gain /= targetCutoff; ?
}

template<class T>
void rsPoleZeroMapper<T>::sPlanePrototypeToLowpass(Complex* prototypePoles, Complex* prototypeZeros, 
  Complex* targetPoles, Complex* targetZeros, int prototypeOrder, T targetCutoff)
{
  for(int i = 0; i < prototypeOrder; i++)
  {
    targetPoles[i] = prototypePoles[i] * targetCutoff;
    targetZeros[i] = prototypeZeros[i] * targetCutoff;
  }
}

template<class T>
void rsPoleZeroMapper<T>::prototypeToAnalogHighpass(Complex* poles, int numPoles, Complex* zeros, 
  int numZeros, T* gain, T targetCutoff)
{
  int i;
  for(i = 0; i < numPoles; i++)
    poles[i] = targetCutoff / poles[i];
  for(i = 0; i < numZeros; i++)
    zeros[i] = targetCutoff / zeros[i];

  // there are numPoles-numZeros zeros at infinity in the lowpass-prototype and for for each such 
  // zero at infinity, we obtain a zero at s=0 in the highpass filter:
  for(i = numZeros; i < numPoles; i++)
    zeros[i] = Complex(0.0, 0.0);

  *gain *= 1; // preliminary \todo adjust gain
}

template<class T>
void rsPoleZeroMapper<T>::sPlanePrototypeToHighpass(Complex* prototypePoles, Complex* prototypeZeros, 
  Complex *targetPoles, Complex *targetZeros, int prototypeOrder, T targetCutoff)
{
  for(int i = 0; i < prototypeOrder; i++)
  {
    targetPoles[i] = targetCutoff / prototypePoles[i];
    if( rsIsInfinite(prototypeZeros[i]) )
      targetZeros[i] = Complex(0.0, 0.0);
    else
      targetZeros[i] = targetCutoff / prototypeZeros[i];
  }
}

template<class T>
void rsPoleZeroMapper<T>::prototypeToAnalogHighShelv(Complex* poles, int numPoles, Complex* zeros, 
  int numZeros, T* /*gain*/, T targetCutoff)
{
  int i;
  for(i = 0; i < numPoles; i++)
    poles[i] = zeros[i] * targetCutoff;
  for(i = 0; i < numZeros; i++)
    zeros[i] = poles[i] * targetCutoff;

  // \todo adjust gain - huh= this is not supposed to work when overwriting the poles in the first 
  // array
}

template<class T>
void rsPoleZeroMapper<T>::sPlanePrototypeToHighShelv(Complex* prototypePoles, Complex* prototypeZeros, 
  Complex* targetPoles, Complex* targetZeros, int prototypeOrder, T targetCutoff)
{
  for(int i = 0; i < prototypeOrder; i++)
  {
    targetPoles[i] = prototypeZeros[i] * targetCutoff;
    targetZeros[i] = prototypePoles[i] * targetCutoff;
  }
}

template<class T>
void rsPoleZeroMapper<T>::prototypeToAnalogBandpass(Complex* poles, int numPoles, Complex* zeros, 
  int numZeros, T* /*gain*/, T targetLowCutoff, T targetHighCutoff)
{
  T wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  T bw = targetHighCutoff - targetLowCutoff;     // bandwidth
  Complex* tmpPoles = new Complex[2*numPoles];
  Complex* tmpZeros = new Complex[2*numPoles];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  int     i;
  Complex tmp1, tmp2;
  for(i = 0; i < numPoles; i++)
  {
    tmp1                     = T(0.5) * bw * poles[i];
    tmp2                     = -bw * poles[i];
    tmp2                     = T(0.25) * tmp2*tmp2;
    tmp2                     = sqrt(tmp2 - wc*wc);
    //tmpPoles[2*i]   = tmp1 + tmp2;
    //tmpPoles[2*i+1] = tmp1 - tmp2;
    tmpPoles[i]              = tmp1 + tmp2;
    tmpPoles[2*numPoles-i-1] = tmp1 - tmp2;
  }
  for(i = 0; i < numZeros; i++)
  {
    tmp1                     = T(0.5) * bw * zeros[i];
    tmp2                     = -bw * zeros[i];
    tmp2                     = T(0.25) * tmp2*tmp2;
    tmp2                     = sqrt(tmp2 - wc*wc);
    //tmpZeros[2*i]   = tmp1 + tmp2;
    //tmpZeros[2*i+1] = tmp1 - tmp2;
    tmpZeros[i]              = tmp1 + tmp2;
    tmpZeros[2*numPoles-i-1] = tmp1 - tmp2;
  }

  // put additional zeros at 0.0 and infinity when there are less finite zeros than finite poles:
  //for(i=numZeros; i<(numPoles-numZeros); i++)
  for(i = numZeros; i < numPoles; i++)
  {
    //tmpZeros[2*i]   = 0.0;
    //tmpZeros[2*i+1] = INF;
    tmpZeros[i]              = 0.0;
    tmpZeros[2*numPoles-i-1] = RS_INF(T);
  }

  //int         order = 2 * max(numPoles, numZeros);

  // copy the content of the temporary arrays into the output arrays:
  for(i = 0; i < 2*numPoles; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // \todo: adjust gain

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;
}

template<class T>
void rsPoleZeroMapper<T>::sPlanePrototypeToBandpass(Complex* prototypePoles, Complex* prototypeZeros, 
  Complex* targetPoles, Complex* targetZeros, int prototypeOrder, T targetLowerCutoff,                      
  T targetUpperCutoff)
{
  T wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  T bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  Complex tmp1, tmp2;
  int k;
  for(k = 0; k < prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = T(0.5) * bw * prototypePoles[k];
    tmp2 = -bw * prototypePoles[k];
    tmp2 = T(0.25) * tmp2*tmp2;
    tmp2 = sqrt(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    if( rsIsInfinite(prototypeZeros[k]) )
    {
      targetZeros[2*k]   = RS_INF(T);
      targetZeros[2*k+1] = 0.0;
    }
    else
    {
      tmp1 = T(0.5) * bw * prototypeZeros[k];
      tmp2 = -bw * prototypeZeros[k];
      tmp2 = T(0.25) * tmp2*tmp2;
      tmp2 = sqrt(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  //// re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  //k = 1;
  //while( k+1 < 2*prototypeOrder )
  //{
  //  rsSwap(targetPoles[k], targetPoles[k+1]); // old
  //  rsSwap(targetZeros[k], targetZeros[k+1]);
  //  k += 4;
  //}
  //// doesn't work anymore since proting to rapt...

  // ...instead, this is needed now:
  sSortFilterRoots(targetPoles, 2*prototypeOrder);
  sSortFilterRoots(targetZeros, 2*prototypeOrder);
    // of course, it would be highly desirable to be able to get rid of a full-blown sorting 
    // procedure here - but it's going to be re-written anyway...
}

template<class T>
void rsPoleZeroMapper<T>::prototypeToAnalogBandstop(Complex* poles, int numPoles, Complex* zeros, 
  int numZeros, T* /*gain*/, T targetLowCutoff, T targetHighCutoff)
{
  T wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  T bw = targetHighCutoff - targetLowCutoff;     // bandwidth

  int         order = 2 * rsMax(numPoles, numZeros);
  Complex* tmpPoles = new Complex[order];
  Complex* tmpZeros = new Complex[order];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  int     i;
  Complex tmp1, tmp2;
  for(i = 0; i < numPoles; i++)
  {
    tmp1            = T(0.5) * bw / poles[i];
    tmp2            = -bw / poles[i];
    tmp2            = T(0.25) * tmp2*tmp2;
    tmp2            = sqrt(tmp2 - wc*wc);
    tmpPoles[2*i]   = tmp1 + tmp2;
    tmpPoles[2*i+1] = tmp1 - tmp2;
  }
  for(i = 0; i < numZeros; i++)
  {
    tmp1            = T(0.5) * bw / zeros[i];
    tmp2            = -bw / zeros[i];
    tmp2            = T(0.25) * tmp2*tmp2;
    tmp2            = sqrt(tmp2 - wc*wc);
    tmpZeros[2*i]   = tmp1 + tmp2;
    tmpZeros[2*i+1] = tmp1 - tmp2;
  }

  // put additional zeros at s=j*wc and s=-j*wc when there are less finite zeros than finite poles:
  Complex j(0.0, 1.0); // imaginary unit
  for(i = numZeros; i < (numPoles-numZeros); i++)
  {
    tmpZeros[2*i]   =  j*wc;
    tmpZeros[2*i+1] = -j*wc;
  }

  // copy the content of the temporary arrays into the output arrays:
  for(i = 0; i < order; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;

  // \todo: adjust gain
}

template<class T>
void rsPoleZeroMapper<T>::sPlanePrototypeToBandreject(Complex* prototypePoles, 
  Complex* prototypeZeros, Complex* targetPoles, Complex* targetZeros, int prototypeOrder, 
  T targetLowerCutoff, T targetUpperCutoff)
{
  T wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  T bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  Complex tmp1, tmp2;
  int k;
  for(k = 0; k < prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = T(0.5) * bw / prototypePoles[k];
    tmp2 = -bw / prototypePoles[k];
    tmp2 = T(0.25) * tmp2*tmp2;
    tmp2 = sqrt(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    Complex j(0.0, 1.0); // imaginary unit
    if( rsIsInfinite(prototypeZeros[k]) )
    {
      targetZeros[2*k]   =  j*wc;
      targetZeros[2*k+1] = -j*wc;
    }
    else
    {
      tmp1 = T(0.5) * bw / prototypeZeros[k];
      tmp2 = -bw / prototypeZeros[k];
      tmp2 = T(0.25) * tmp2*tmp2;
      tmp2 = sqrt(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  //// re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  //k = 1;
  //while( k+1 < 2*prototypeOrder )
  //{
  //  rsSwap(targetPoles[k], targetPoles[k+1]);
  //  int kp = (k-1)/4;
  //  if( !isInfinite(prototypeZeros[kp]) )
  //    rsSwap(targetZeros[k], targetZeros[k+1]); // don't want to swap the generated zeros at +-j*wc
  //  k += 4;
  //}
  //// doesn't work anymore since porting to rapt..

  // instead, this is needed now:
  sSortFilterRoots(targetPoles, 2*prototypeOrder);
  sSortFilterRoots(targetZeros, 2*prototypeOrder);
  // todo: rewrite these design routines, clean up code and adopt conventions that avoid the need
  // for re-ordering
}

template<class T>
void rsPoleZeroMapper<T>::bilinearAnalogToDigital(Complex* poles, int numPoles, Complex* zeros, 
  int numZeros, T sampleRate, T* /*gain*/)
{
  int     i;
  Complex z;
  T scaler = T(0.5)/sampleRate;
  Complex one(1, 0);

  for(i = 0; i < numPoles; i++)
  {
    // we don't check against infinity (as we do with zeros) because infinite poles should actually 
    // never occur in practice
    z        = scaler * poles[i];
    poles[i] = (one+z)/(one-z);
  }

  for(i = 0; i < numZeros; i++)
  {
    if( rsIsInfinite(zeros[i]) )
      zeros[i] = Complex(-1.0, 0.0);
    else
    {
      z        = scaler * zeros[i];
      zeros[i] = (one+z)/(one-z);
    }
  }
  // maybe factor out into bilinearAnalogToDigital(Complex* roots, int numRoots, T sampleRate)
}

template<class T>
void rsPoleZeroMapper<T>::zLowpassToLowpass(Complex* /*z*/, Complex* /*p*/, T* /*k*/,
  Complex* /*zNew*/, Complex* /*pNew*/, T* /*kNew*/, int /*N*/, T /*wc*/, T /*wt*/)
{
  // not yet implemented
  //int dummy = 0;
}

template<class T>
void rsPoleZeroMapper<T>::sSortFilterRoots(Complex* r, int N)
{
  rsHeapSort(r, N, &sFilterRootBefore);
  // todo: make sure that roots with positive imag parts are before those with negative 
  // ...but maybe not - later, we don't care about the sign anyway (i think)
}

template<class T>
bool rsPoleZeroMapper<T>::sFilterRootBefore(const Complex& r1, const Complex& r2)
{
  if(abs(r1.imag()) > abs(r2.imag()))
    return true;   // sort roots by (descending, absolute) frequency
  return false;
}

/*
ToDo:
-Implement also the discrete versions of the frequency transformations (Constantinides formulas)
-Implement this: https://vicanek.de/articles/BiquadFits.pdf and maybe MZTi, too
 -Maybe use as matching frequencies, 0, fs/2, fc, fu, fu/2 where fu is the lower bandedge freq.
  If there is no such thing, maybe use fu = fc/2. Avoid using the upper bandedge because that may 
  actually go above fs/2
-Maybe we can also derive formulas that match the magnitude only at 4 frequencies and match the
 phase at fc. ...or maybe match magnitude at 3 frequencies and phase at fc and fs/2
-Maybe rename bilinearAnalogToDigital to bilinear_s2z - or have a seperate class 
 TimeDiscretizaionMapper (or something) and the methods can be named bilinear (for s -> z) and
 invBilinear (for z -> s)
-Implement also Impulse-Invariant, (improved) Matched-Z, 3-point Magnitude Match, 
 5-point Magnitude Match, etc.
-Implement Phase-Invariant-Method (PIM) and Magnitude-Invariant-Method (MIM) to map poles and zeros
 from s-plane to z-plane (and maybe back, if possible), 
 see https://soar.wichita.edu/handle/10057/1564


-Maybe wite a similar class that operates on biquads rather than poles and zeros. It could have 
 functions like bilinearAnalogToDigital(T B0, T B1, T B2, T A0, T A1, T A2, T& b0, ...)


*/
