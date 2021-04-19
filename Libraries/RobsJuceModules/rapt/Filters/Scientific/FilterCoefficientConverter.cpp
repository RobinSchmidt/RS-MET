template<class T>
void rsFilterCoefficientConverter<T>::directFormToLatticeFir(T *directFormCoeffs, int order,
  T *reflectionCoeffs)
{
  int     N    = order;               // number of direct-form and reflection coefficients
  T* a    = new T[N+1];
  T* aOld = new T[N+1];
  T* k    = reflectionCoeffs;
  T  scaler;
  int     i, m;

  // the recursion assumes the leading FIR coefficient to be unity, so we normalize the FIR-vector 
  // by its leading coefficient and remember it as overall gain-factor:
  T  gain = directFormCoeffs[0];

  // copy the sign inverted, gain normalized direct-form coefficients into an internal array 
  // (Eq 6.50):
  scaler  = T(1) / gain;
  a[0]    = directFormCoeffs[0] * scaler;
  aOld[0] = a[0];
  for(i = 1; i <= N; i++)
  {
    a[i]    = -directFormCoeffs[i] * scaler;
    aOld[i] = a[i];
  }

  // do the recursion (Eq. 6.56c):
  for(i = N; i >= 1; i--)
  {
    k[i] = a[i];

    scaler = T(1) / (T(1) - k[i]*k[i]);
    for(m = 1; m <= (i-1); m++)
      a[m] = (aOld[m] + k[i]*aOld[i-m]) * scaler;

    // update the aOld-array:
    for(m = 0; m <= N; m++)
      aOld[m] = a[m];
  }

  // store the gain factor in the first element of the reflection coefficient vector:
  k[0] = gain;

  // invert the sign of the reflection coeffs to follow the MatLab convention (but leave the 
  // additional gain factor in k[0] as is):
  for(i = 1; i <= N; i++)
    k[i] = -k[i];

  delete[] a;
  delete[] aOld;
}

template<class T>
void rsFilterCoefficientConverter<T>::latticeToDirectFormFir(T *reflectionCoeffs, int order,
  T *directFormCoeffs)
{
  int     N    = order;            // number of direct-form and reflection coefficients
  T* a    = directFormCoeffs;
  T* aOld = new T[N+1];
  T* k    = reflectionCoeffs;
  int     i, m;

  // copy the direct-form coefficients into an internal array:
  for(i = 0; i <= N; i++)
    aOld[i] = a[i];

  // do the recursion (Eq. 6.54a 6.54b, tweaked to accomodate for MatLabs sign-convention):
  for(i = 1; i <= N; i++)
  {
    a[i] = -k[i];

    for(m = 1; m <= (i-1); m++)
      a[m] = aOld[m] + k[i]*aOld[i-m];

    // update the aOld-array:
    for(m = 0; m <= N; m++)
      aOld[m] = a[m];
  }

  // in the first element of the reflection coefficient vector is the overall gain-factor - we 
  // need to scale the FIR vector with this factor and also use it as leading coefficient, we 
  // also need to invert the sign:
  for(i = 1; i <= N; i++)
    a[i] = -a[i] * k[0];
  a[0] = k[0];

  delete[] aOld;
}
// todo:
// -to which book do these Eq. comments refer? i think, it was the oppenheim/schafer? -> figure 
//  out and add to the comments. 
// -i think, the recursions can be used also for the IIR part - if so, rename the functions 
//  accordingly. 
// -i think, these functions do the same thing as matlab's tf2latc/latc2tf functions - maybe add 
//  that info to the comments


template<class T>
void rsFilterCoefficientConverter<T>::polesAndZerosToBiquadCascade(Complex *poles, int numPoles,
  Complex *zeros, int numZeros, T *b0, T *b1, T *b2, T *a1, T *a2,
  bool lastStageIsFirstOrder)
{
  int order = 2*rsMax(numPoles, numZeros);
  if(lastStageIsFirstOrder)
    order -= 1;

  int numBiquads;
  if(rsIsEven(order))
    numBiquads = order/2;
  else
    numBiquads = (order+1)/2;

  int b;
  for(b = 0; b < numBiquads; b++)
  {
    b0[b] = T(1);
    b1[b] = T(-2) * zeros[b].real();
    b2[b] = zeros[b].real() * zeros[b].real() + zeros[b].imag() * zeros[b].imag();
    a1[b] = T(-2) * poles[b].real();
    a2[b] = poles[b].real() * poles[b].real() + poles[b].imag() * poles[b].imag();
  }

  // overwrite the coefficients of the last stage, when it must be a first order stage:
  if(lastStageIsFirstOrder)
  {
    b1[numBiquads-1] = -zeros[numBiquads-1].real();
    b2[numBiquads-1] = 0.0;
    a1[numBiquads-1] = -poles[numBiquads-1].real();
    a2[numBiquads-1] = 0.0;
  }
}

template<class T>
void rsFilterCoefficientConverter<T>::polesAndZerosToBiquadCascade(Complex *poles, Complex *zeros,
  int order, T *b0, T *b1, T *b2, T *a1, T *a2)
{
  int numBiquads;
  if(rsIsEven(order))
    numBiquads = order/2;
  else
    numBiquads = (order+1)/2;
  // use (order+1)/2 regardless - if order is even, truncation will take place and the result is 
  // the same

  int b;
  //for(b=0; b<numBiquads; b++)
  for(b = 0; b < order/2; b++)
  {
    b0[b] = 1.0;
    b1[b] = -(zeros[2*b]+zeros[2*b+1]).real();
    b2[b] =  (zeros[2*b]*zeros[2*b+1]).real();
    a1[b] = -(poles[2*b]+poles[2*b+1]).real();
    a2[b] =  (poles[2*b]*poles[2*b+1]).real();

    /*
    b1[b] = -2.0 * zeros[2*b].re;
    b2[b] = zeros[2*b].re * zeros[2*b].re + zeros[2*b].im * zeros[2*b].im;
    a1[b] = -2.0 * poles[2*b].re;
    a2[b] = poles[2*b].re * poles[2*b].re + poles[2*b].im * poles[2*b].im;
    */
  }

  // overwrite the coefficients of the last stage, when it must be a first order stage:
  if(rsIsOdd(order))
  {
    b0[b]            = 1.0;
    b1[numBiquads-1] = -zeros[order-1].real();
    b2[numBiquads-1] = 0.0;
    a1[numBiquads-1] = -poles[order-1].real();
    a2[numBiquads-1] = 0.0;
  }
}

template<class T>
void rsFilterCoefficientConverter<T>::biquadCascadeToDirectForm(int numBiquads, T *b0,
  T *b1, T *b2, T *a1, T *a2, T *b, T *a)
{
  // aquire and initialize memory for intermediate results and accumulation:
  int         N      = 2*numBiquads+1;      // number of direct form coefficients
  long double *tmp   = new long double[N];
  long double *aAccu = new long double[N];
  long double *bAccu = new long double[N];
  long double *aQuad = new long double[3];
  long double *bQuad = new long double[3];
  // i used long double to avoid roundoff error accumulation - but on 64 bit systems, "long double"
  // is the same as "double", so this doesn't help anything anymore - so get rid of these temporary
  // accumulation buffers and avoid memory allocation


  int i;
  for(i = 0; i < N; i++)
  {
    aAccu[i] = 0.0;
    bAccu[i] = 0.0;
  }
  aAccu[0] = 1.0;
  bAccu[0] = 1.0;

  // calculate the direct form coefficients by recursively applying convolution to the result of 
  // the previous convolution with the new quadratic factor:
  for(i = 0; i < numBiquads; i++)
  {
    // aquire the current quadratic factor:
    aQuad[0] = 1.0;
    aQuad[1] = a1[i];
    aQuad[2] = a2[i];
    bQuad[0] = b0[i];
    bQuad[1] = b1[i];
    bQuad[2] = b2[i];

    // convolve the current quadratic factor with the result of the previous convolution:
    rsArrayTools::copy(aAccu, tmp, N);
    rsArrayTools::convolve(tmp, N-2, aQuad, 3, aAccu);
    rsArrayTools::copy(bAccu, tmp, N);
    rsArrayTools::convolve(tmp, N-2, bQuad, 3, bAccu);
      // this can be optimized .... xLength does not need always be N-2 (can be shorter in early 
      // iterations)
  }

  // copy (and typecast) the accumulated coefficents into the output arrays:
  for(i = 0; i < N; i++)
  {
    a[i] = (T)aAccu[i];
    b[i] = (T)bAccu[i];
  }

  // release temporarily aquired memory:
  delete[] tmp;
  delete[] aAccu;
  delete[] bAccu;
  delete[] aQuad;
  delete[] bQuad;
}

template<class T>
T rsFilterCoefficientConverter<T>::getBiquadMagnitudeAt(T b0, T b1, T b2, T a1, T a2, T omega)
{
  T cosOmega  = cos(omega);
  T cos2Omega = cos(T(2)*omega);
  T num       = b0*b0 + b1*b1 + b2*b2 + T(2)*cosOmega*(b0*b1 + b1*b2) + T(2)*cos2Omega*b0*b2;
  T den       = T(1)  + a1*a1 + a2*a2 + T(2)*cosOmega*(a1    + a1*a2) + T(2)*cos2Omega*a2;
  return sqrt(num/den);
}

template<class T>
void rsFilterCoefficientConverter<T>::normalizeBiquadStages(T *b0, T *b1, T *b2, 
  T *a1, T *a2, T omega, int numStages, T gainFactor)
{
  T m; // magnitude of current stage
  T scaler;
  for(int s = 0; s < numStages; s++)
  {
    // divide all the feedforward coefficients by the current magnitude at omega:
    m       = getBiquadMagnitudeAt(b0[s], b1[s], b2[s], a1[s], a2[s], omega);
    scaler  = gainFactor/m;
    b0[s]  *= scaler;
    b1[s]  *= scaler;
    b2[s]  *= scaler;
  }
}

