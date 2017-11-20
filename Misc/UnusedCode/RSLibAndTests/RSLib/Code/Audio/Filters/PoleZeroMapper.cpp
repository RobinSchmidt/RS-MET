using namespace RSLib;

void rsPoleZeroMapper::sLowpassToLowshelf(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double G0, double G)
{
  if( G0 == 0.0 )
  {
    rsCopyBuffer(z, zNew, N);
    rsCopyBuffer(p, pNew, N);
    *kNew = *k * G;
    return;
  } 

  double kTmp;
  rsComplexDbl *zTmp = new rsComplexDbl[2*N];
  rsComplexDbl *pTmp = new rsComplexDbl[2*N];
  double *a     = new double[N+1];
  double *b     = new double[N+1];
  double *bS    = new double[2*N+1];  // shelving mgnitude squared numerator coeffs

  // we design a boost filter in any case - a dip filter will later be obtained by swapping poles 
  // and zeros:
  bool dip = false;
  if( G < G0 )
  {
    dip = true;
    G   = 1.0 / G;
    G0  = 1.0 / G0;
  }

  // scale poles/zeros/gain of the lowpass such that the resulting low-shelving filter has a gain 
  // of sqrt(G*G0) at unit frequency:
  double GB = sqrt(G0*G);                       // bandwidth gain  
  double GC = sqrt((GB*GB-G0*G0)/(G*G-G0*G0));  // desired lowpass gain at unit frequency

  //rsPrototypeDesigner::scaleToMatchGainAtUnity(z, p, k, zTmp, pTmp, &kTmp, N, GC);  
    // commented for test

  // for debug:
  rsCopyBuffer(z, zTmp, N);
  rsCopyBuffer(p, pTmp, N);
  kTmp = *k;

  // obtain magnitude-squared numerator polynomial for shelving filter:
  rootsToCoeffs(zTmp, b, N);
  rootsToCoeffs(pTmp, a, N);
  rsPrototypeDesigner::shelvingMagSqrNumeratorFromLowpassTransfer(b, a, kTmp, N, G0, G, bS);

  // obtain zeros, poles and gain of the new shelving filter:
  //findPolynomialRoots(bS, 2*N, &zTmp[-1]);
  //int numLeftZeros = onlyLeftHalfPlane(zTmp, zTmp, 2*N);
  rsPrototypeDesigner::getLeftHalfPlaneRoots(bS, zTmp, 2*N);  
    // test - replaces the 2 commented lines above - seems to work

  // if we make a dip-filter, poles and zeros exchange roles:
  if( dip == false )
  {
    rsCopyBuffer(zTmp, zNew, N);
    rsCopyBuffer(pTmp, pNew, N);
    *kNew = sqrt(fabs(bS[2*N]));
  }
  else
  {
    rsCopyBuffer(zTmp, pNew, N);
    rsCopyBuffer(pTmp, zNew, N);
    *kNew = 1.0 / sqrt(fabs(bS[2*N]));
  }

  delete[] zTmp;
  delete[] pTmp;
  delete[] a;
  delete[] b;
  delete[] bS;
}

void rsPoleZeroMapper::sLowpassToLowpass(rsComplexDbl *z, rsComplexDbl *p, double *k,                                      
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wc)
{
  for(int n = 0; n < N; n++)
    zNew[n] = wc * z[n];
  for(int n = 0; n < N; n++)
    pNew[n] = wc * p[n];
  int nz = rsGetNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew = *k * pow(wc, N-nz);
}

void rsPoleZeroMapper::sLowpassToHighpass(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc)
{
  for(int n = 0; n < N; n++)
  {
    if( r[n].isInfinite() )
      rNew[n] = rsComplexDbl(0.0, 0.0);
    else
      rNew[n] = wc / r[n];
  }
}
void rsPoleZeroMapper::sLowpassToHighpass(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wc)
{
  sLowpassToHighpass(z, zNew, N, wc);
  sLowpassToHighpass(p, pNew, N, wc);
  rsComplexDbl n = rsProductOfFiniteFactors(z, N);  
  rsComplexDbl d = rsProductOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).re;  

  // sign-factor: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p); 
  // ..or maybe not? ..consult octave/matlab implementation
  int nz = rsGetNumFiniteValues(z, N);
  if( rsIsOdd(N+nz) )
    *kNew = - *kNew;
}

void rsPoleZeroMapper::sLowpassToBandpass(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc, 
  double bw)
{
  rsComplexDbl tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1 = 0.5 * bw * r[n];
    tmp2 = -bw * r[n];

    // this is common with the bandreject case - factor out:
    tmp2      = 0.25 * tmp2*tmp2;
    tmp2      = rsSqrtC(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}
void rsPoleZeroMapper::sLowpassToBandpass(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wl, double wu)
{
  double wc = sqrt(wl*wu);  // center (radian) frequency
  double bw = wu-wl;        // bandwidth
  sLowpassToBandpass(z, zNew, N, wc, bw);
  sLowpassToBandpass(p, pNew, N, wc, bw);
  for(int n = 0; n < N; n++)
  {
    if( z[n].isInfinite() )
    {
      zNew[n]   = rsComplexDbl(0.0, 0.0);
      zNew[N+n] = rsComplexDbl(rsInfDouble, 0.0);
    }
  }
  int nz = rsGetNumFiniteValues(z, N); // number of finite zeros in prototype
  *kNew  = *k * pow(bw, N-nz); 
}

void rsPoleZeroMapper::sLowpassToBandreject(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc, 
  double bw)
{
  rsComplexDbl tmp1, tmp2;
  for(int n = 0; n < N; n++)
  {
    tmp1 = 0.5 * bw / r[n];
    tmp2 = -bw / r[n];

    // this is common with the bandpass case - factor out:
    tmp2      = 0.25 * tmp2*tmp2;
    tmp2      = rsSqrtC(tmp2 - wc*wc);
    rNew[n]   = tmp1 + tmp2;
    rNew[N+n] = tmp1 - tmp2;
  }
}
void rsPoleZeroMapper::sLowpassToBandreject(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wl, double wu)
{
  double wc = sqrt(wl*wu);  // center (radian) frequency
  double bw = wu-wl;        // bandwidth
  sLowpassToBandreject(z, zNew, N, wc, bw);
  sLowpassToBandreject(p, pNew, N, wc, bw);
  rsComplexDbl sz = rsComplexDbl(0.0, wc);
  for(int n = 0; n < N; n++)
  {
    if( z[n].isInfinite() )
    {
      zNew[n]   =  sz;
      zNew[N+n] = -sz;
    }
  }
  rsComplexDbl n = rsProductOfFiniteFactors(z, N);  
  rsComplexDbl d = rsProductOfFiniteFactors(p, N);
  *kNew     = *k * (n/d).re;  

  // sign-factor - see m-file: kNew = k * (-1)^(length(p) + length(z)) * prod(z) / prod(p);
  int nz = rsGetNumFiniteValues(z, N);
  if( rsIsOdd(N+nz) )
    *kNew = - *kNew;
}

void rsPoleZeroMapper::prototypeToAnalogLowpass(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double *gain, double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = poles[i] * targetCutoff;
  for(i=0; i<numZeros; i++)
    zeros[i] = zeros[i] * targetCutoff;

  // \todo adjust gain
}

void rsPoleZeroMapper::sPlanePrototypeToLowpass(rsComplexDbl *prototypePoles, 
  rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
  int prototypeOrder, double targetCutoff)
{
  for(int i = 0; i < prototypeOrder; i++)
  {
    targetPoles[i] = prototypePoles[i] * targetCutoff;
    targetZeros[i] = prototypeZeros[i] * targetCutoff;
  }
}

void rsPoleZeroMapper::prototypeToAnalogHighpass(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double *gain, double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = targetCutoff / poles[i];
  for(i=0; i<numZeros; i++)
    zeros[i] = targetCutoff / zeros[i];

  // there are numPoles-numZeros zeros at infinity in the lowpass-prototype and for for each such 
  // zero at infinity, we obtain a zero at s=0 in the highpass filter:
  for(i=numZeros; i<numPoles; i++)
    zeros[i] = rsComplexDbl(0.0, 0.0);

  // \todo adjust gain
}

void rsPoleZeroMapper::sPlanePrototypeToHighpass(rsComplexDbl *prototypePoles, 
  rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
  int prototypeOrder, double targetCutoff)
{
  for(int i=0; i<prototypeOrder; i++)
  {
    targetPoles[i] = targetCutoff / prototypePoles[i];
    if( prototypeZeros[i].isInfinite() )
      targetZeros[i] = rsComplexDbl(0.0, 0.0);
    else
      targetZeros[i] = targetCutoff / prototypeZeros[i];
  }
}

void rsPoleZeroMapper::prototypeToAnalogHighShelv(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double *gain, double targetCutoff)
{
  int i;
  for(i=0; i<numPoles; i++)
    poles[i] = zeros[i] * targetCutoff;
  for(i=0; i<numZeros; i++)
    zeros[i] = poles[i] * targetCutoff;

  // \todo adjust gain - huh? this is not supposed to work when overwriting the poles in the first 
  // array
}

void rsPoleZeroMapper::sPlanePrototypeToHighShelv(rsComplexDbl *prototypePoles, 
  rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
  int prototypeOrder, double targetCutoff)
{
  for(int i=0; i<prototypeOrder; i++)
  {
    targetPoles[i] = prototypeZeros[i] * targetCutoff;
    targetZeros[i] = prototypePoles[i] * targetCutoff;
  }
}

void rsPoleZeroMapper::prototypeToAnalogBandpass(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double *gain, double targetLowCutoff, double targetHighCutoff)
{
  double wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  double bw = targetHighCutoff - targetLowCutoff;     // bandwidth
  rsComplexDbl* tmpPoles = new rsComplexDbl[2*numPoles];
  rsComplexDbl* tmpZeros = new rsComplexDbl[2*numPoles];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  int     i;
  rsComplexDbl tmp1, tmp2;
  for(i=0; i<numPoles; i++)
  {
    tmp1                     = 0.5 * bw * poles[i];
    tmp2                     = -bw * poles[i];
    tmp2                     = 0.25 * tmp2*tmp2;
    tmp2                     = rsSqrtC(tmp2 - wc*wc);
    //tmpPoles[2*i]   = tmp1 + tmp2;
    //tmpPoles[2*i+1] = tmp1 - tmp2;
    tmpPoles[i]              = tmp1 + tmp2;
    tmpPoles[2*numPoles-i-1] = tmp1 - tmp2;
  }
  for(i=0; i<numZeros; i++)
  {
    tmp1                     = 0.5 * bw * zeros[i];
    tmp2                     = -bw * zeros[i];
    tmp2                     = 0.25 * tmp2*tmp2;
    tmp2                     = rsSqrtC(tmp2 - wc*wc);
    //tmpZeros[2*i]   = tmp1 + tmp2;
    //tmpZeros[2*i+1] = tmp1 - tmp2;
    tmpZeros[i]              = tmp1 + tmp2;
    tmpZeros[2*numPoles-i-1] = tmp1 - tmp2;
  }

  // put additional zeros at 0.0 and infinity when there are less finite zeros than finite poles:
  //for(i=numZeros; i<(numPoles-numZeros); i++)
  for(i=numZeros; i<numPoles; i++)
  {
    //tmpZeros[2*i]   = 0.0;
    //tmpZeros[2*i+1] = rsInfDouble;
    tmpZeros[i]              = 0.0;
    tmpZeros[2*numPoles-i-1] = rsInfDouble;
  }

  //int         order = 2 * max(numPoles, numZeros);

  // copy the content of the temporary arrays into the output arrays:
  for(i=0; i<2*numPoles; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // \todo: adjust gain

  // for debug:
  rsComplexDbl pDbg[40];
  rsComplexDbl zDbg[40];
  int k;
  for(k=0; k<2*numPoles; k++)
  {
    pDbg[k] = poles[k];
    zDbg[k] = zeros[k];
  }

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;
}

void rsPoleZeroMapper::sPlanePrototypeToBandpass(rsComplexDbl *prototypePoles, 
  rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
  int prototypeOrder, double targetLowerCutoff, double targetUpperCutoff)
{
  double wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  double bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  rsComplexDbl tmp1, tmp2;
  int k;
  for(k=0; k<prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = 0.5 * bw * prototypePoles[k];
    tmp2 = -bw * prototypePoles[k];
    tmp2 = 0.25 * tmp2*tmp2;
    tmp2 = rsSqrtC(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    if( prototypeZeros[k].isInfinite() )
    {
      targetZeros[2*k]   = rsInfDouble;
      targetZeros[2*k+1] = 0.0;
    }
    else
    {
      tmp1 = 0.5 * bw * prototypeZeros[k];
      tmp2 = -bw * prototypeZeros[k];
      tmp2 = 0.25 * tmp2*tmp2;
      tmp2 = rsSqrtC(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  // re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  k = 1;
  while( k+1 < 2*prototypeOrder )
  {
    rsSwap(targetPoles[k], targetPoles[k+1]);
    rsSwap(targetZeros[k], targetZeros[k+1]);
    k += 4;
  }
}

void rsPoleZeroMapper::prototypeToAnalogBandstop(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double *gain, double targetLowCutoff, 
  double targetHighCutoff)
{
  double wc = sqrt(targetLowCutoff*targetHighCutoff); // center (radian) frequency
  double bw = targetHighCutoff - targetLowCutoff;     // bandwidth

  int         order = 2 * rsMax(numPoles, numZeros);
  rsComplexDbl* tmpPoles = new rsComplexDbl[order];
  rsComplexDbl* tmpZeros = new rsComplexDbl[order];

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  int     i;
  rsComplexDbl tmp1, tmp2;
  for(i=0; i<numPoles; i++)
  {
    tmp1            = 0.5*bw / poles[i];
    tmp2            = -bw / poles[i];
    tmp2            = 0.25 * tmp2*tmp2;
    tmp2            = rsSqrtC(tmp2 - wc*wc);
    tmpPoles[2*i]   = tmp1 + tmp2;
    tmpPoles[2*i+1] = tmp1 - tmp2;
  }
  for(i=0; i<numZeros; i++)
  {
    tmp1            = 0.5*bw / zeros[i];
    tmp2            = -bw / zeros[i];
    tmp2            = 0.25 * tmp2*tmp2;
    tmp2            = rsSqrtC(tmp2 - wc*wc);
    tmpZeros[2*i]   = tmp1 + tmp2;
    tmpZeros[2*i+1] = tmp1 - tmp2;
  }

  // put additional zeros at s=j*wc and s=-j*wc when there are less finite zeros than finite poles:
  rsComplexDbl j(0.0, 1.0); // imaginary unit
  for(i=numZeros; i<(numPoles-numZeros); i++)
  {
    tmpZeros[2*i]   =  j*wc;
    tmpZeros[2*i+1] = -j*wc;
  }

  // copy the content of the temporary arrays into the output arrays:
  for(i=0; i<order; i++)
  {
    zeros[i] = tmpZeros[i];
    poles[i] = tmpPoles[i];
  }

  // free dynamically allocate memory:
  delete[] tmpPoles;
  delete[] tmpZeros;

  // \todo: adjust gain
}

void rsPoleZeroMapper::sPlanePrototypeToBandreject(rsComplexDbl *prototypePoles, 
  rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
  int prototypeOrder, double targetLowerCutoff, double targetUpperCutoff)
{
  double wc = sqrt(targetLowerCutoff*targetUpperCutoff); // center (radian) frequency
  double bw = targetUpperCutoff - targetLowerCutoff;     // bandwidth

  // for each pole (or zero) in the prototype, we obtain a pair of poles (or zeros) in the 
  // bandpass-filter:
  rsComplexDbl tmp1, tmp2;  
  int k;
  for(k=0; k<prototypeOrder; k++)
  {
    // transform poles:
    tmp1 = 0.5 * bw / prototypePoles[k];
    tmp2 = -bw / prototypePoles[k];
    tmp2 = 0.25 * tmp2*tmp2;
    tmp2 = rsSqrtC(tmp2 - wc*wc);
    targetPoles[2*k]   = tmp1 + tmp2;
    targetPoles[2*k+1] = tmp1 - tmp2;

    // transform zeros:
    rsComplexDbl j(0.0, 1.0); // imaginary unit
    if( prototypeZeros[k].isInfinite() )
    {
      targetZeros[2*k]   =  j*wc;
      targetZeros[2*k+1] = -j*wc;
    }
    else
    {
      tmp1 = 0.5 * bw / prototypeZeros[k];
      tmp2 = -bw / prototypeZeros[k];
      tmp2 = 0.25 * tmp2*tmp2;
      tmp2 = rsSqrtC(tmp2 - wc*wc);
      targetZeros[2*k]   = tmp1 + tmp2;
      targetZeros[2*k+1] = tmp1 - tmp2;
    }
  }

  // re-arrange poles and zeros such that complex conjugate pairs (again) occupy successive slots:
  k = 1;
  while( k+1 < 2*prototypeOrder )
  {
    rsSwap(targetPoles[k], targetPoles[k+1]);
    int kp = (k-1)/4;
    if( !prototypeZeros[kp].isInfinite() )
      rsSwap(targetZeros[k], targetZeros[k+1]); // don't want to swap the generated zeros at +-j*wc
    k += 4;
  }
}

void rsPoleZeroMapper::bilinearAnalogToDigital(rsComplexDbl *poles, int numPoles, 
  rsComplexDbl *zeros, int numZeros, double sampleRate, double *gain)
{
  int i;
  rsComplexDbl z;
  double scaler = 0.5/sampleRate;

  for(i=0; i<numPoles; i++)
  {
    // we don't check against infinity (as we do with zeros) because infinite poles should 
    // actually never occur in practice
    z        = scaler * poles[i];
    poles[i] = (1.0+z)/(1.0-z);
  }

  for(i=0; i<numZeros; i++)
  {
    if( zeros[i].isInfinite() )
      zeros[i] = rsComplexDbl(-1.0, 0.0);
    else
    {
      z        = scaler * zeros[i];
      zeros[i] = (1.0+z)/(1.0-z);
    }
  }
}

void rsPoleZeroMapper::zLowpassToLowpass(rsComplexDbl *z, rsComplexDbl *p, double *k, 
  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wc, double wt)
{
  rsError("not yet implemented");
  // doo stufff....
}
