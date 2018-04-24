using namespace RSLib;

double RSLib::rsCentroid(double *x, int N)
{
  double num = 0.0;
  double den = 0.0;
  for(int n = 0; n < N; n++)
  {
    num += (n+1)*x[n];
    den += x[n];
  }
  return num/den - 1.0;
}

double RSLib::rsCentroidOfEnergy(double *x, int N)
{
  double num = 0.0;
  double den = 0.0;
  double x2;
  for(int n = 0; n < N; n++)
  {
    x2   = x[n]*x[n];
    num += (n+1)*x2;
    den += x2;
  }
  return num/den - 1.0;
}
// maybe make a generalization that uses pow(fabs(x[n]), a) instead of just x[n] or x^2[n]

/*
void RSLib::rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p)
{
  // The nonlinear system of 2 equations:
  // y0 = a * sin(p)
  // y1 = a * sin(p+w)
  // with known y0, y1, w and unknown a, p is solved for a and p:
  double s, c;
  rsSinCos(w, &s, &c);
  *p = atan2(y0*s, y1-y0*c);
  s  = sin(*p);
  if( fabs(s) > EPS )
    *a = y0 / s;
  else
    *a = y1 / sin(*p + w);
}
*/

void RSLib::rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p)
{
  // The nonlinear system of 2 equations:
  // y0 = a * sin(p)
  // y1 = a * sin(p+w)
  // with known y0, y1, w and unknown a, p is solved for a and p:
  double s, c, s2;
  rsSinCos(w, &s, &c);
  *p = atan2(y0*s, y1-y0*c);
  s  = sin(*p);
  s2 = sin(*p + w);
  if( fabs(s) > fabs(s2) )
    *a = y0 / s;
  else
  {
    if(s2 == 0.0)
      *a = 0.0;
    else
      *a = y1 / s2;
  }
  // less efficient than commented version above (always computes sin(*p + w) instead of just when
  // necessarry) - but might be more accurate numerically? ...tests are needed
}

double RSLib::rsSineFrequency(double y0, double y1, double y2, double small)
{
  rsAssert( fabs(y1) > small * (fabs(y0)+fabs(y2)), "y1 (numerically) zero is not allowed");

  // There's a recursion for the sine y[n] = a1*y[n-1] - y[n-2] where a1 = 2*cos(w) and the states
  // y[n-1], y[n-2] are initialized as y[n-1] = A * sin(p - w), y[n-2] = A * sin(p - 2*w) which in 
  // our notation here translates to y2 = a1*y1 - y0. This leads to a1 = (y0+y2)/y1 and 
  // w = acos(a1/2):
  return acos(0.5*(y0+y2)/y1); 
}

double RSLib::rsSineFrequencyAtCore(double *x, int N, int n0, double small)
{
  if( n0 <= 0 )
    n0 = 1;
  else if( n0 >= N-1 )
    n0 = N-2;
  if( fabs(x[n0]) <= small * (fabs(x[n0-1])+fabs(x[n0+1])) )
  {
    if( n0 == N-2 )
      n0--;
    else
      n0++;
  }
  return rsSineFrequency(x[n0-1], x[n0], x[n0+1]);
}


double refineFrequencyEstimate(double *x, int N, int n, double w)
{
  double c = 2.0 / 3.0;  // found empirically
  double w1;             // w measured with 1 sample offset
  double wr;             // refined w

  //bool forward = true; 
  bool forward = n+2 < N-1; // x[n+2] must be a valid index

  if( forward )
  {
    w1 = rsSineFrequency(x[n], x[n+1], x[n+2]);
    //wr =  c*w + (1-c)*w1;
    wr = w1 + c*(w-w1); // equivalent to c*w + (1-c)*w1
  }
  else
  {
    //w1 = rsSineFrequency(x[n-2], x[n-1], x[n]);
    //wr =  w - c*(w-w1); // no!
    //wr =  w1 - c*(w-w1); // no!
    //wr =  c*w + (1-c)*w1; // no
    //wr =  (1-c)*w + c*w1; // even worse
    wr =  w; // for backward measurement, i did not yet find a good formula, so we just return
              // the input value
  }

  // at high w, we may get outliers that are worse than the unrefined values, so we return the 
  // refined value only if the relative difference to the input value is below some threshold 
  // (which was found empirically):
  double thresh = 0.0001;
  double relativeDelta = (w-wr)/w;
  if( rsAbs(relativeDelta) < thresh )
    return wr;
  else
    return w;
}

double RSLib::rsSineFrequencyAt(double *x, int N, int n0, bool refine)
{
  // find indices of maxima/minima before and after the sample-index in question, taking into 
  // account cases where our index n0 is not surrounded by peaks/valley in which case we take the 
  // two to the left or right to n0 and use these for extrapolation
  int nL = rsFindPeakOrValleyLeft( x, N, n0);
  int nR = rsFindPeakOrValleyRight(x, N, n0);
  if( nL == -1 )
  {
    nL = nR;
    nR = rsFindPeakOrValleyRight(x, N, nL+1);
  }
  if( nR == -1 )
  {
    nR = nL;
    nL = rsFindPeakOrValleyLeft(x, N, nR-1);
  }

  //if( nL == nR )
  //  return rsSineFrequencyAtCore(x, N, nL);

  // compute instantaneous frequencies at nL and nR:
  double wL = rsSineFrequencyAtCore(x, N, nL); 
  double wR = rsSineFrequencyAtCore(x, N, nR); 

  if( refine == true )
  {
    wL = refineFrequencyEstimate(x, N, nL, wL);
    wR = refineFrequencyEstimate(x, N, nR, wR);
  }

  // compute frequency at n0 by linear interpolation (or extrapolation):
  double w  = rsInterpolateLinear((double)nL, (double)nR, wL, wR, (double)n0);
  return w;



  // in case of upward sweeps, the frequency will be slightly underestimated and in case of 
  // downward sweeps it will be overestimated. maybe this bias can be predicted and compensated 
  // for...

  // experimental: try to remove measurement bias by predicting and compensating (to use one of these,
  // comment out the "return f" above and uncomment one of the formulas here:
  double dw = wR - wL;
  double rw = wR / wL;
  //return w + 0.5*dw*w;
  //return w + 0.5*dw*w*rw; // error around 10/10^5 for 1000-2000 sweep, biased to -
  //return w + 0.25*dw*w*rw; // error around 5/10^5 for 1000-2000 sweep, biased to +
  //return w + 0.3*dw*w*rw; // error around 5/10^5 for 1000-2000 sweep, unbiased
  //return w + 0.15*dw*w*rw; // seems unbiased for 100->200, but has ouliers at the maxima
  //return w + 0.5*rw*w - 0.5*w;
  //return w + 0.5*rw; // no
}

double RSLib::rsSinePhaseAt(double *x, int N, int n0, double w)
{
  rsAssert(N >= 2, "Signal too short");

  //double w = rsSineFrequencyAt(x, N, n0);
  double a, p;
  if( n0 <= N-2 )
    rsSineAmplitudeAndPhase(x[n0], x[n0+1], w, &a, &p);
  else
  {
    rsSineAmplitudeAndPhase(x[n0-1], x[n0], w, &a, &p);
    p = rsWrapToInterval(p+w, -PI, PI);
  }
  return p;
}

double RSLib::rsSinePhaseAt(double *x, int N, int n0)
{
  double w = rsSineFrequencyAt(x, N, n0);
  return rsSinePhaseAt(x, N, n0, w);
}



// move all these somewhere else:
template <class T>
inline bool isUpwardZeroCrossing(T x, int N, int n0)
{
  if( x[n0] <= 0 && x[n0+1] > 0 )
    return true;
  else
    return false;
}

template <class T>
inline bool isDownwardZeroCrossing(T x, int N, int n0)
{
  if( x[n0] >= 0 && x[n0+1] < 0 )
    return true;
  else
    return false;
}

int RSLib::rsFindZeroNear(double *x, int N, int n0, int searchDirection, bool upward, 
  bool downward)
{
  int n;
  bool up, down;
  if( searchDirection == +1 )
  {
    for(n = n0; n < N-1; n++)
    {
      up   = isUpwardZeroCrossing(  x, N, n) && upward;
      down = isDownwardZeroCrossing(x, N, n) && downward;
      if( up || down )
        return n;
    }
  }
  else if( searchDirection == -1 )
  {
    for(n = n0; n >= 0; n--)
    {
      up   = isUpwardZeroCrossing(  x, N, n) && upward;
      down = isDownwardZeroCrossing(x, N, n) && downward;
      if( up || down )
        return n;
    }
  }
  else
  {
    rsError("searchDirection parameter must be -1 or +1");
  }
  return -1;
}

// intergrate into RSLib:
double rsSubSamplePrecisionZeroCrossing(double *x, int N, int n0, int p)
{
  // maybe factor out a function rsFractionalPartOfZeroCrossing

  double *a = new double[2*p+2]; // polynomial coefficients for interpolant 
  double nf;                     // fractional part of zero-crossing sample-index

  // adjust (reduce) p, when we are at the borders of the signal to avoid access violation
  p = rsMin(p, n0);
  p = rsMin(p, N-n0-2);

  nf = x[n0]/(x[n0]-x[n0+1]);    // linear zero
  if( p > 0 )
  {
    // refine linear zero estimate by Newton iteration on a higher order interpolating 
    // polynomial:
    rsInterpolatingPolynomial(a, -p, 1, &x[n0-p], 2*p+2);
    nf = getRootNear(nf, a, 2*p+1, 0.0, 1.0);
  }

  delete[] a;
  return n0 + nf;
}

double RSLib::rsSinePhaseAtViaZeros(double *x, int N, int n0, int precision)
{
  // find integer zero crossings before and after n0
  //....
  int nzl, nzr;   // left and right zero crossing index
  double zl, zr;  // left and right zero crossing with subsample precision
  double p = 0;   // computed phase

  if(x[n0] == 0.0)
  {

  }
  else if(x[n0] > 0.0) 
  {
    // to the left, there's an upward zero crossing, to the right a downward one:
    nzl = rsFindZeroNear(x, N, n0, -1, true, false);
    nzr = rsFindZeroNear(x, N, n0, +1, false, true);

    zl  = rsSubSamplePrecisionZeroCrossing(x, N, nzl, precision);
    zr  = rsSubSamplePrecisionZeroCrossing(x, N, nzr, precision);

    if( zl < n0 && zr > n0 )
      p = PI*(n0-zl)/(zr-zl);
    else
    {
      rsError("not yet implemented");
    }
  }
  else
  {
    // x[n0] < 0.0
    // to the left, there's an downward zero crossing, to the right a upward one
    nzl = rsFindZeroNear(x, N, n0, -1, false, true);
    nzr = rsFindZeroNear(x, N, n0, +1, true, false);

    zl  = rsSubSamplePrecisionZeroCrossing(x, N, nzl, precision);
    zr  = rsSubSamplePrecisionZeroCrossing(x, N, nzr, precision);

    if( zl < n0 && zr > n0 )
      p = PI*(n0-zl)/(zr-zl) - PI;
    else
    {
      rsError("not yet implemented");
    }
  }

  return p;
}

double RSLib::rsSineShiftAmount(double *x, int N, int n0, double p0, double w)
{
  // compute desired phase-shift: 

  //double p  = rsSinePhaseAt(x, N, n0, w);  // instantaneous phase at n0
  double p  = rsSinePhaseAtViaZeros(x, N, n0, 3); // test - with new function

  double dp = p - p0;                      // desired phase-shift

  // map to the range -pi...pi (never shift more than half a cycle):
  if( dp < -PI )
    dp += 2*PI;
  else if( dp > PI )
    dp -= 2*PI;
    // todo: maybe provide an optional integer parameter to always shift left if its negative, always 
    // right if it's positive and use the minimum distance if it's 0

  // return time-shift in samples:
  return dp / w;
}

double RSLib::rsSineShiftAmount(double *x, int N, int n0, double p0)
{
  double w = rsSineFrequencyAt(x, N, n0);  // instantaneous normalized radian frequency at n0
  return rsSineShiftAmount(x, N, n0, p0, w);
}