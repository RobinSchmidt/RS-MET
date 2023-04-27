#include "SaturationExperiments.h"

void powRatioParametricSigmoid()
{
  // Test for a famliy of sigmoid curves that have a parameter "p" that determines, how fast the 
  // sigmoid rises into its saturating regime. The function plotted here is given by
  // y = x / (1 + |x|^p)^(1/p)

  static const int N = 1001; // number of points
  double xMin = 0.0;   
  double xMax = 4.0;
  double p = 1.0;            // parameter p for 1st curve, doubles for each subsequent curve
  int numCurves = 7;         // number of curves

  // test with large input:
  //double test = rsNormalizedSigmoids::powRatio(-1000, 1000);

  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;
  for(int i = 0; i < numCurves; i++)
  {
    for(int n = 0; n < N; n++)
      y[n] = rsNormalizedSigmoidsD::powRatio(x[n], p);
    plt.addDataArrays(N, x, y);
    p *= 2;
  }
  plt.plot();
}

void parametricSigmoid()
{
  // Test of the rsParametricSigmoid class. We plot a family of curves for different values of the
  // y1-parameter (the function value of the sigmoid at x=1).

  static const int N = 501; // number of points
  double xMin = -3.0;
  double xMax =  3.0;
  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;

  double f1 = 0.5;                  // desired value at x=1 for 1st curve
  rsParametricSigmoidD sig;
  sig.setPiecewiseBreakpoint(0.75); // value for f1, above which we use piecewise functions
  int numCurves = 11;
  for(int i = 0; i < numCurves; i++)
  {
    sig.setValueAt1(f1);
    for(int n = 0; n < N; n++)
      y[n] = sig.getValue(x[n]);
    plt.addDataArrays(N, x, y);
    f1 += 0.05;
  }
  plt.plot();
}

void parametricSigmoid2()
{
  // We compare the parametric sigmoid with some other sigmoid functions (such as tanh, atan, etc.)
  // to see, how well these other sigmoids can be approximated with the parametric sigmoid.

  static const int N = 501; // number of points
  double xMin = 0.0;
  double xMax = 3.0;
  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;

  rsParametricSigmoidD sig;
  sig.setPiecewiseBreakpoint(0.75);   
  int n;

  // test:
  sig.setThreshold(0.0);
  sig.setThreshold(0.2);
  sig.setThreshold(0.25);
  sig.setThreshold(0.4);
  sig.setThreshold(0.5);
  sig.setThreshold(0.6);
  sig.setThreshold(0.75);
  sig.setThreshold(0.8);
  sig.setThreshold(1.0); 


  // haprdclipper for reference:
  for(n = 0; n < N; n++)
    y[n] = rsNormalizedSigmoidsD::clip(x[n]);
  plt.addDataArrays(N, x, y); 

  //// tanh: 
  //sig.setValueAt1(tanh(1));
  //for(n = 0; n < N; n++)
  //  y[n] = sig.getValue(x[n]);
  //plt.addDataArrays(N, x, y);
  //for(n = 0; n < N; n++)
  //  y[n] = tanh(x[n]);
  //plt.addDataArrays(N, x, y);

  // quartic:
  sig.setValueAt1(rsPositiveSigmoidsD::quartic(1));
  for(n = 0; n < N; n++)
    y[n] = sig.getValue(x[n]);
  plt.addDataArrays(N, x, y);
  for(n = 0; n < N; n++)
    y[n] = rsPositiveSigmoidsD::quartic(x[n]);
  plt.addDataArrays(N, x, y);

  // cubicRational (is actually a special case, so the match should be perfect):
  sig.setValueAt1(0.75);
  for(n = 0; n < N; n++)
    y[n] = sig.getValue(x[n]);
  plt.addDataArrays(N, x, y);
  for(n = 0; n < N; n++)
    y[n] = rsPositiveSigmoidsD::cubicRational(x[n]);
  plt.addDataArrays(N, x, y);
  // ok - works

  plt.plot();
}


// move to RSLib...
double rsQuinticParametricSigmoid(double x, double s)
{
  // ...this is currently only the positive range function

  //if(x > s)
  //  return 1.0;

  // compute polynomial coefficients from desired saturation level s (factor out):
  double s2 = s*s;      // s^2
  double s3 = s*s2;     // s^3
  double s4 = s2*s2;    // s^4
  double s5 = s3*s2;    // s^5
  double a3 = (10-6*s)/s3;
  double a4 = (8*s-15)/s4;
  double a5 = (6-3*s) /s5;
  // ...optimize computations...use si = 1/s, si2 = si*si (=1/s^2), etc.


  // compute output:
  double x2 = x*x;   // x^2
  double x3 = x*x2;  // x^3
  double y  = x + a3*x3 + a4*x2*x2 + a5*x2*x3;
  return y;

  // todo: check, if for s=2, it reduces to our quartic sigmoid, i.e. a5=0 and a3,a4 the same as
  // for the quartic - yes, it is.
  // for s<2, we need a piecewise defined function, using the quartic internally on a 
  // scaled/shifted x' ...like in the rsSaturator class
  // problem, for s>2 (for example s=3), the polynomial is not monotonic anymore on the interval 
  // 0..s
  // maybe, we need to construct our polynomial differently, such that the derivative is a 
  // sum-of-squares on the interval 0..s
}
void quinticParametricSigmoid()
{
  // We select a saturation-level s, above which we just set f(x)=1, below that value (x < s), we 
  // use a quintic polynomial of the form f(x) = x + a3*x^3 + a4*x^4 + a5*x^5. By construction,
  // it satisfies f(0)=0, f'(0)=1, f''(0)=0. We choose the coefficients a3,a4,a5 according to the
  // requirements: f(s)=1, f'(s)=0, f''(s)=0.

  // This leads to the system of equations:
  // s + a3*s^3 + a4*s^4 + a5*s^5 = 1, 1 + 3*a3*s^2 + 4*a4*s^3 + 5*a5*s^4 = 0,
  // 3*a3 + 6*a4*s + 10*a5*s^2 = 0
  //

  // to wolfram alpha:
  // solve(s+a_3*s^3+a_4*s^4+a_5*s^5=1, 1+3*a_3*s^2+4*a_4*s^3+5*a_5*s^4=0, 3*a_3+6*a_4*s+10*a_5*s^2=0; a_3,a_4,a_5)
  // gives
  // a_3 = -(2 (3 s-5))/s^3 and a_4 = (8 s-15)/s^4 and a_5 = -(3 (s-2))/s^5 and s!=0

  static const int N = 1001; // number of points
  double xMin = 0.0;   
  double xMax = 4.0;
  double s = 3.0;            // the saturation level

  //double a3, a4, a5;         // polynomial coefficients

  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;
  for(int n = 0; n < N; n++)
    y[n] = rsQuinticParametricSigmoid(x[n], s);
  plt.addDataArrays(N, x, y);
  plt.plot();

  // Observations: with s = 3, the function is nonmonotonic - we must make sure that the polynomial
  // is monotonic on the interval 0..s
  // see here:
  // http://math.stackexchange.com/questions/60610/polynomial-fitting-where-polynomial-must-be-monotonically-increasing
  // we need theorem 6 from here:
  // http://stellar.mit.edu/S/course/6/sp10/6.256/courseMaterial/topics/topic2/lectureNotes/lecture-10/lecture-10.pdf
  // applied to the derivative of our polynomial. Our derivative is of even order, so we must use
  // the 1st line of the theorem.

  // this might be tricky - maybe just try an odd septic polynomial of the form:
  // f(x) = x + a3*x^3 + a5*x^5 + a7*x^7
  // ...maybe somehow the odd symmetry has the side effect of making the derivative monotonic?
  // a square of a polynomial is necessarily an even polynomial (the square of any function is an 
  // even function). choosing an odd polynomial, the derivative will be even, so it's somehow 
  // plausible that oddness is a necessary condition for the derivative to be even (everywhere)

  // for the quintic, try to leave out th a4*x^4 term and match only the 1st derivative
}

// move to RSLib...
double rsSepticParametricSigmoid(double x, double s)
{
  // ...this is currently only the positive range function

  if(x > s)
    return 1.0;

  // factor out: septicParametricCoeffs(s, a3, a5, a7)
  double s2, sn, a3, a5, a7;
  s2  = s*s;                  // s^2
  sn  = s2*s;                 // s^3
  a3  = (35-24*s)/(8*sn);     // (35-24s)/(8s^3)
  sn *= s2;                   // s^5
  a5  = (12*s-21)/(4*sn);     // (12s-21)/(4s^5)
  sn *= s2;                   // s^7
  a7  = (15-8*s)/(8*sn);      // (15-8s)/(8s^7)
  // may be optimized - avoid divisions by using 1/s^n

  double x2, x3, x5;
  x2 = x*x;    // x^2
  x3 = x*x2;   // x^3
  x5 = x3*x2;  // x^5
  return x + a3*x3 + a5*x5 + a7*x5*x2;
}
void septicParametricSigmoid()
{
  // We select a saturation-level s, above which we just set f(x)=1, below that value (x < s), we 
  // use a septic polynomial of the form f(x) = x + a3*x^3 + a5*x^5 + a7*x^7. By construction,
  // it satisfies f(0)=0, f'(0)=1, f''(0)=0. We choose the coefficients a3,a5,a7 according to the
  // requirements: f(s)=1, f'(s)=0, f''(s)=0.

  // This leads to the system of equations:
  // s + a3*s^3 + a5*s^5 + a7*s^7 = 1, 1 + 3*a3*s^2 + 5*a5*s^4 + 7*a7*s^6 = 0,
  // 3*a3 + 10*a5*s^2 + 21*a7*s^4 = 0
  //

  // to wolfram alpha:
  // solve(s+a_3*s^3+a_5*s^5+a_7*s^7=1, 1+3*a_3*s^2+5*a_5*s^4+7*a_7*s^6=0, 3*a_3+10*a_5*s^2+21*a_7*s^4=0; a_3,a_5,a_7)
  // gives
  // a_3 = (35-24 s)/(8 s^3) and a_5 = (3 (4 s-7))/(4 s^5) and a_7 = (15-8 s)/(8 s^7) and s!=0

  static const int N = 1001; // number of points
  double xMin = 0.0;   
  double xMax = 4.0;
  double s = 3.0;            // the saturation level

  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;
  for(int n = 0; n < N; n++)
    y[n] = rsSepticParametricSigmoid(x[n], s);
  plt.addDataArrays(N, x, y);
  plt.plot();

  // Observations: 
  // damn - it also becomes nonmonotonic - of course, it isn't enough for a polynomial to be even
  // - we also nee the coefficients to be positive.
}

void saturator()
{
  // Test of the rsSaturator class.

  int N = 1001;         // number of points
  double xMin = -2.2;   
  double xMax = +2.2;

  // create and set up the saturator:
  rsSaturatorDD sat;
  sat.setLowerThreshold(0.8);     // lower halfwave threshold
  sat.setUpperThreshold(0.5);     // upper halfwave threshold

  int n;
  vector<double> x(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);

  // hyperbolic tangent:
  vector<double> yTanh(N);
  sat.setSaturationFunctions(&rsTanh);
  for(n = 0; n < N; n++)
    yTanh[n] = sat.getSample(x[n]);

  // linear:
  vector<double> yLinear(N);
  sat.setSaturationFunctions(&rsPositiveSigmoidsD::linear);
  for(n = 0; n < N; n++)
    yLinear[n] = sat.getSample(x[n]);

  // rational:
  vector<double> yRational(N);
  sat.setSaturationFunctions(&rsPositiveSigmoidsD::rational);
  for(n = 0; n < N; n++)
    yRational[n] = sat.getSample(x[n]);

  // cubic:
  vector<double> yCubic(N);
  sat.setSaturationFunctions(&rsPositiveSigmoidsD::cubic);
  for(n = 0; n < N; n++)
    yCubic[n] = sat.getSample(x[n]);

  // quartic:
  vector<double> yQuartic(N);
  sat.setSaturationFunctions(&rsPositiveSigmoidsD::quartic);
  for(n = 0; n < N; n++)
    yQuartic[n] = sat.getSample(x[n]);

  // sixtic:
  vector<double> ySixtic(N);
  sat.setSaturationFunctions(&rsPositiveSigmoidsD::hexic);
  for(n = 0; n < N; n++)
    ySixtic[n] = sat.getSample(x[n]);

  //// identity:
  //vector<double> yIdent(N);
  //sat.setSaturationFunctions(&rsIdentity);
  //for(n = 0; n < N; n++)
  //  yIdent[n] = sat.getSample(x[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &yRational[0], &yTanh[0], &yLinear[0], &yCubic[0], &yQuartic[0], &ySixtic[0]);
  //plt.addDataArrays(N, &x[0], &yCubic[0], &yQuartic[0], &ySixtic[0]);
  plt.plot();

  // The sixtic curve looks (and sounds) more similar to the cubic one instead. The additionl
  // smoothness does does really pay off sonically - quartic seems to be the sweet spot.

  // Ideas for other waveshaping functions (not necessarily sigmoid):
  // -f(x) = x^(2*n+1) / (1 + x^(2*n))
  //  https://www.desmos.com/calculator/vcij6pvtzq
  //  n < 0: linear around 0 (with unit slope), squashes down to zero towards +-inf
  //  n = 0: linear (with slope 1/2 -> fix that!)
  //  n > 0: squash around zero, asymptotically linear (with unit slope)
  //  -Try fracctional n, maybe needs to insert abs() in denominator to avoid singularity
}

void sigmoidScaleAndShift()
{
  // Assuming we have some kind of sigmoid function f(x), that has the following properties:
  // f(0) = 0, f'(0) = 1, f(-inf) = -1, f(inf) = +1, for example: f(x) = tanh(x), we apply scaling
  // and shifting of the input and output such that the output range is mapped to arbitrary values.

  static const int N = 501;     // number of samples/values
  double xMin = -1;            // minimum x-value
  double xMax = +1;            // maximum x-value
  double lo   = -0.2;             // minimum of output range
  double hi   = +0.4;             // maximum of output range

  // compute center and y-width of the sigmoid (lo = c-w, hi = c+w):
  double c = 0.5*(lo+hi);  // center
  double w = hi-lo;        // width (in y-direction)

  // allocate arrays:
  double x[N];                  // array of x-values
  double yClip[N];              // hardclipped
  double yCubic[N];             // cubic polynomial
  double yQuartic[N];           // quartic polynomial
  double ySixtic[N];            // quartic polynomial
  double yTanh[N];              // hyperbolic tangent
  double yAtan[N];              // arcus tangens

  // compute scale and translation coefficients for input x and output y (i think: ty = c, 
  // sy = w/2):
  double sx, tx, sy, ty;
  rsRangeConversionCoefficients(  lo,   hi, -1.0, +1.0, &sx, &tx);
  rsRangeConversionCoefficients(-1.0, +1.0,   lo,   hi, &sy, &ty);

  // fill the arrays:
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  double t;  // temporary
  for(int n = 0; n < N; n++)
  {
    t = sx*x[n] + tx;  // transformed input to the normalized sigmoid function
    //t = x[n]; // just for test
    yClip[n]    = ty + sy * rsNormalizedSigmoidsD::clip(t);
    yCubic[n]   = ty + sy * rsNormalizedSigmoidsD::cubic(t);  
    yQuartic[n] = ty + sy * rsNormalizedSigmoidsD::quartic(t);
    ySixtic[n]  = ty + sy * rsNormalizedSigmoidsD::hexic(t);  
    yTanh[n]    = ty + sy * rsTanh(t);
    yAtan[n]    = ty + sy * rsNormalizedSigmoidsD::atan(t);  
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, x, yClip, yCubic, ySixtic, yQuartic, yTanh, yAtan);
  plt.plot();

  // Observations:
  // With this scale/shift strategy, the tansformed sigmoid does not satisfy f(0)=0 anymore. To
  // make this condition hold again, without affecting the asymptotic limits, we would need to 
  // shift the whole function slightly to the right. 
  // ToDo: figure out, how much we must rightshift...
  //
  // The curves are ordered in the plot th plot in such a way that the earlier ones reach their
  // saturation level more quickly. This seems to be the most important factor for the overall 
  // sound or behavior.
}

double evaluateQuartic(double x, double a1, double a2, double a3, double a4)
{
  // simple helper function to evaluate y = a1*x + a2*x^2 + a3*x^3 + a4*x^4
  double x2 = x*x;
  return a1*x + a2*x2 + a3*x*x2 + a4*x2*x2;
}
void monotonicQuarticCoeffs(double k, double *p2, double *p3, double *p4)
{
  rsAssert(k <= 6); // no solution for k > 6, 1.5<=k<=4: solution without turning points

  // intermediate variables:
  double b0, b1, b12, b0b1;
  b0   = sqrt(k);
  b1   = -2*b0 + sqrt(12-2*k);
  //b1   = -2*b0 - sqrt(12-2*k);  // negative solution gives saddle
  b12  = b1*b1;  // b1^2
  b0b1 = b0*b1;  // b0*b1

  // polynomial coefficients:
  *p2 = ( 2*b0b1 - k)   * 0.5;
  *p3 = (-2*b0b1 + b12) * (1/3.);
  *p4 = (        - b12) * 0.25;
}
double saturateQuarticMonotonic(double x, double k)
{
  if(x > k)
    return 1.0;
  double p2, p3, p4;
  monotonicQuarticCoeffs(k, &p2, &p3, &p4);
  return evaluateQuartic(x/k, k, p2, p3, p4);
}
void quarticMonotonic()
{
  // test:
  double p2, p3, p4;
  monotonicQuarticCoeffs(1.7, &p2, &p3, &p4);


  static const int N = 1001;
  double xMin = 0.0;
  double xMax = 4.0;
  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  GNUPlotter plt;
  double k  = 1.5;    // saturation level goes from 1.5 to...
  int numCurves = 6;  // ...4 (increment is 0.5)
  for(int i = 0; i < numCurves; i++)
  {
    for(int n = 0; n < N; n++)
      y[n] = saturateQuarticMonotonic(x[n], k);
    plt.addDataArrays(N, x, y);
    k += 0.5;
  }

  //// for comparison - add tanh and atan:
  //for(int n = 0; n < N; n++)
  //  y[n] = tanh(x[n]);
  //plt.addDataArrays(N, x, y);
  //for(int n = 0; n < N; n++)
  //  y[n] = rsNormalizedSigmoids::atan(x[n]);
  //plt.addDataArrays(N, x, y);

  plt.plot();
}

/** Soft clipping with threshold t */
//double softClipSextic(double x, double t)
//{
//  if(x <= t)
//    return x;
//  else
//  {
//    x = (x-t) / (1-t);
//    return t + (1-t) * rsNormalizedSigmoids::sixtic(x);
//  }
//}

void sigmoidPrototypes()
{
  // plots various normalized prototype sigmoids that can be used as basic building blocks for
  // saturators, soft-clippers, etc.

  static const int N = 1001;
  double xMin = 0.0;
  double xMax = 2.0;
  double t    = 0.5;    // threshold below which the function is linear (for the softclipper)

  double x[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);

  int n;
  double yHard[N], yCubic[N], yQuartic[N], yHexic[N], ySoft[N], yTanh[N];

  for(n = 0; n < N; n++)
    yHard[n] = rsNormalizedSigmoidsD::clip(x[n]);
  for(n = 0; n < N; n++)
    yCubic[n] = rsPositiveSigmoidsD::cubic(x[n]);
  for(n = 0; n < N; n++)
    yQuartic[n] = rsPositiveSigmoidsD::quartic(x[n]);
  for(n = 0; n < N; n++)
    yHexic[n] = rsPositiveSigmoidsD::hexic(x[n]);
  for(n = 0; n < N; n++)
    ySoft[n] = rsPositiveSigmoidsD::softClipHexic(x[n], t);
  for(n = 0; n < N; n++)
    yTanh[n] = rsNormalizedSigmoidsD::tanh(x[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, x, yHard);
  //plt.addDataArrays(N, x, yCubic);
  //plt.addDataArrays(N, x, yQuartic);
  plt.addDataArrays(N, x, yHexic);
  plt.addDataArrays(N, x, ySoft);
  //plt.addDataArrays(N, x, yTanh);
  plt.plot();
}

// value of the polynomial a1*x + a4*x^4 + a5*x^5 + a6*x^6
double sixticValue(double x, double a1, double a4, double a5, double a6)
{
  double x2 = x*x;    // x^2
  double x4 = x2*x2;  // x^4
  return a1*x + a4*x4 + a5*x4*x + a6*x4*x2;
}
// Computes the coefficients a4, a5, a6 for a 6th order polynomial 
// f(x) = a1*x + a4*x^4 + a5*x^5 + a6*x^6 from a given coefficient a1, such that the following
// properties hold: at some input value x = s = 2/a1, we have f(s)=1 and the derivatives 
// f', f'', f''' are zero at x = s
void sixticCoeffs(double a1, double &a4, double &a5, double &a6)
{
  double s  = 2/a1;   // found graphically by plotting the objective function g(s)
  double s2 = s*s;    // s^2
  double s4 = s2*s2;  // s^4
  a4 = -5 / s4;       // -5/s^4
  a5 =  6 / (s4*s);   //  6/s^5
  a6 = -2 / (s4*s2);  // -2/s^6
}
// Computes value of a sixtic polynomial f(x) = k*x + a4*x^4 + a5*x^5 + a6*x^6 which satifies the
// following conditions: f(0)=0, f'(0)=k, f''(0)=f'''(0)=0, f(s)=1, f'(s)=f''(s)=f'''(s)=0 for
// some value s which is given by s = 2/k. This is useful as nonlinear section for a piecewise 
// saturation function that should match derivatives at the junction points to the linear part and
// the constant part.
double sixticValue(double x, double k)
{
  double a4, a5, a6;
  sixticCoeffs(k, a4, a5, a6);
  return sixticValue(x, k, a4, a5, a6);
}
// ...simplified function for k=1
double sixticValue(double x)
{
  double x2 = x*x;    // x^2
  double x4 = x2*x2;  // x^4
  return x - 0.3125*x4 + 0.1875*x4*x - 0.03125*x4*x2;
}
// move these to RSLib....
void sixticPositive()
{
  // For positive x-values use a polynomial:

  // f(x) = k*x + a4*x^4 + a5*x^5 + +a6*x^6, for x <= s (s: saturation level)
  // = 1                                for x >  s

  // and use an odd symmetrized version of it for negative x-values. The 1st 3 derivatives are:

  // f0(x) = k*x +    a4*x^4 +    a5*x^5 +     a6*x^6
  // f1(x) = k   +  4*a4*x^3 +  5*a5*x^4 +   6*a6*x^5
  // f2(x) =       12*a4*x^2 + 20*a5*x^3 +  30*a6*x^4
  // f3(x) =       24*a4*x   + 60*a5*x^2 + 120*a6*x^3

  // By construction, the 1st derivative at x=0 equals k and derivatives 2,3 at x=0 equal zero. We 
  // additionally impose the following requirements: at some (yet unknown) saturation level x=s, we
  // want the function value f0(s) to be unity and derivatives 1..3 to be zero. So our equations 
  // are:

  // f0(s) = 1 = k*s +    a4*s^4 +    a5*s^5 +     a6*s^6
  // f1(s) = 0 = k   +  4*a4*s^3 +  5*a5*s^4 +   6*a6*s^5
  // f2(s) = 0 =       12*a4*s^2 + 20*a5*s^3 +  30*a6*s^4
  // f3(s) = 0 =       24*a4*s   + 60*a5*s^2 + 120*a6*s^3

  // we can simplify a bit by dividing f2, f3 by s^2 and s respectively:

  // f0(s) = 1 = k*s +    a4*s^4 +    a5*s^5 +     a6*s^6
  // f1(s) = 0 = k   +  4*a4*s^3 +  5*a5*s^4 +   6*a6*s^5
  // f2(s) = 0 =       12*a4     + 20*a5*s   +  30*a6*s^2
  // f3(s) = 0 =       24*a4     + 60*a5*s   + 120*a6*s^2

  // We have 4 unknowns s, a4, a5, a6 and the system is linear in a4, a5, a6 and nonlinear in s. 
  // We can cast then problem as a one-dimensional root-finding problem for s by defining the 
  // objective function:

  // g(s) = k*s + a4*s^4 + a5*s^5 + a6*s^6 - 1

  // which results from the 1st equation. We want to find a root, such that g(s) = 0. Let's choose 
  // a fixed value for k (which is the slope at the origin), for example set k=1. To evaluate g(s) 
  // for some given s, we first need to compute a4, a5, a6 from the linear system of equations 
  // (which results from the remaining 3 equations):

  // |4*s^3  5*s^4   6*s^5| * |a4| = |-k|
  // |12    20*s    30*s^2|   |a5|   | 0|
  // |24    60*s   120*s^2|   |a6|   | 0|

  // Having found a4, a5, a6, we evaluate g(s) = k*s + a4*s^4 + a5*s^5 + a6*s^6 - 1 and return the 
  // value. This function g(s) can be used as input for a general 1D root finding algorithm for 
  // some specified k. 

  // By inspecting the objective function, it seems that the solution is simply given by s=2/k, so 
  // we don't really need the root-finder. Just set s=2/k and solve the remaining linear system:

  //  s = 2/k,
  // -k = 4*a_4*s^3  +  5*a_5*s^4 +   6*a_6*s^5,
  //  0 = 12*a_4     + 20*a_5*s   +  30*a_6*s^2,
  //  0 = 24*a_4     + 60*a_5*s   + 120*a_6*s^2,

  // throwing at wolfram alpha:
  // solve(s=2/k, -k=4*a_4*s^3+5*a_5*s^4+6*a_6*s^5, 0=12*a_4+20*a_5*s+30*a_6*s^2, 0=24*a_4+60*a_5*s+120*a_6*s^2; a_4,a_5,a_6)
  // gives:
  // a_4 = -5/s^4 and a_5 = 6/s^5 and a_6 = -2/s^6 and s!=0 and k = 2/s

  // To get rid of having to symmetrize the function around x=0, we could also use a polynomial 
  // that already is symmetrical around x=0 by construction, for example:
  // f(x) = k*x + a5*x^5 + a7*x^7 + +a9*x^9
  // and use the same procedure with the modified polynomial. As and additional benefit, the 4th 
  // derivative at x=0 is zero here as well.

  static const int N = 500;
  double k    = 1.0;    // slope at the origin:
  double xMin = 0.0;
  double xMax = 3.0/k;

  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    y[n] = sixticValue(x[n], k);
    //y[n] = rsSaturateSixtic(x[n]); // just a test
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.plot();
}
