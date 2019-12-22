void rsFilterCoefficientConverter::directFormToLatticeFir(double *directFormCoeffs, int order,
  double *reflectionCoeffs)
{
  int     N    = order;               // number of direct-form and reflection coefficients
  double* a    = new double[N+1];
  double* aOld = new double[N+1];
  double* k    = reflectionCoeffs;
  double  scaler;
  int     i, m;

  // the recursion assumes the leading FIR coefficient to be unity, so we normalize the FIR-vector 
  // by its leading coefficient and remember it as overall gain-factor:
  double  gain = directFormCoeffs[0];

  // copy the sign inverted, gain normalized direct-form coefficients into an internal array 
  // (Eq 6.50):
  scaler  = 1.0 / gain;
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

    scaler = 1.0 / (1.0 - k[i]*k[i]);
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

void rsFilterCoefficientConverter::latticeToDirectFormFir(double *reflectionCoeffs, int order,
  double *directFormCoeffs)
{
  int     N    = order;            // number of direct-form and reflection coefficients
  double* a    = directFormCoeffs;
  double* aOld = new double[N+1];
  double* k    = reflectionCoeffs;
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

void rsFilterCoefficientConverter::polesAndZerosToBiquadCascade(Complex *poles, int numPoles,
  Complex *zeros, int numZeros, double *b0, double *b1, double *b2, double *a1, double *a2,
  bool lastStageIsFirstOrder)
{
  int order = 2*RAPT::rsMax(numPoles, numZeros);
  if(lastStageIsFirstOrder)
    order -= 1;

  int numBiquads;
  if(RAPT::rsIsEven(order))
    numBiquads = order/2;
  else
    numBiquads = (order+1)/2;

  int b;
  for(b = 0; b < numBiquads; b++)
  {
    b0[b] = 1.0;
    b1[b] = -2.0 * zeros[b].re;
    b2[b] = zeros[b].re * zeros[b].re + zeros[b].im * zeros[b].im;
    a1[b] = -2.0 * poles[b].re;
    a2[b] = poles[b].re * poles[b].re + poles[b].im * poles[b].im;
  }

  // overwrite the coefficients of the last stage, when it must be a first order stage:
  if(lastStageIsFirstOrder)
  {
    b1[numBiquads-1] = -zeros[numBiquads-1].re;
    b2[numBiquads-1] = 0.0;
    a1[numBiquads-1] = -poles[numBiquads-1].re;
    a2[numBiquads-1] = 0.0;
  }
}

void rsFilterCoefficientConverter::polesAndZerosToBiquadCascade(Complex *poles, Complex *zeros,
  int order, double *b0, double *b1, double *b2, double *a1, double *a2)
{
  int numBiquads;
  if(RAPT::rsIsEven(order))
    numBiquads = order/2;
  else
    numBiquads = (order+1)/2;

  int b;
  //for(b=0; b<numBiquads; b++)
  for(b = 0; b < order/2; b++)
  {
    b0[b] = 1.0;
    b1[b] = -(zeros[2*b]+zeros[2*b+1]).re;
    b2[b] =  (zeros[2*b]*zeros[2*b+1]).re;
    a1[b] = -(poles[2*b]+poles[2*b+1]).re;
    a2[b] =  (poles[2*b]*poles[2*b+1]).re;

    /*
    b1[b] = -2.0 * zeros[2*b].re;
    b2[b] = zeros[2*b].re * zeros[2*b].re + zeros[2*b].im * zeros[2*b].im;
    a1[b] = -2.0 * poles[2*b].re;
    a2[b] = poles[2*b].re * poles[2*b].re + poles[2*b].im * poles[2*b].im;
    */
  }

  // overwrite the coefficients of the last stage, when it must be a first order stage:
  if(RAPT::rsIsOdd(order))
  {
    b0[b]            = 1.0;
    b1[numBiquads-1] = -zeros[order-1].re;
    b2[numBiquads-1] = 0.0;
    a1[numBiquads-1] = -poles[order-1].re;
    a2[numBiquads-1] = 0.0;
  }
}

void rsFilterCoefficientConverter::biquadCascadeToDirectForm(int numBiquads, double *b0,
  double *b1, double *b2, double *a1, double *a2, double *b, double *a)
{
  // aquire and initialize memory for intermediate results and accumulation:
  int         N      = 2*numBiquads+1;      // number of direct form coefficients
  long double *tmp   = new long double[N];
  long double *aAccu = new long double[N];
  long double *bAccu = new long double[N];
  long double *aQuad = new long double[3];
  long double *bQuad = new long double[3];
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
    RAPT::rsArrayTools::copy(aAccu, tmp, N);
    RAPT::rsArrayTools::convolve(tmp, N-2, aQuad, 3, aAccu);
    RAPT::rsArrayTools::copy(bAccu, tmp, N);
    RAPT::rsArrayTools::convolve(tmp, N-2, bQuad, 3, bAccu);
      // this can be optimized .... xLength does not need always be N-2 (can be shorter in early 
      // iterations)
  }

  // copy (and typecast) the accumulated coefficents into the output arrays:
  for(i = 0; i < N; i++)
  {
    a[i] = (double)aAccu[i];
    b[i] = (double)bAccu[i];
  }

  // release temporarily aquired memory:
  delete[] tmp;
  delete[] aAccu;
  delete[] bAccu;
  delete[] aQuad;
  delete[] bQuad;
}

double rsFilterCoefficientConverter::getBiquadMagnitudeAt(double b0, double b1, double b2, 
  double a1, double a2, double omega)
{
  double cosOmega  = cos(omega);
  double cos2Omega = cos(2.0*omega);
  double num       = b0*b0 + b1*b1 + b2*b2 + 2.0*cosOmega*(b0*b1 + b1*b2) + 2.0*cos2Omega*b0*b2;
  double den       = 1.0   + a1*a1 + a2*a2 + 2.0*cosOmega*(a1    + a1*a2) + 2.0*cos2Omega*a2;
  return sqrt(num/den);
}

void rsFilterCoefficientConverter::normalizeBiquadStages(double *b0, double *b1, double *b2, 
  double *a1, double *a2, double omega, int numStages, double gainFactor)
{
  double m; // magnitude of current stage
  double scaler;
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
