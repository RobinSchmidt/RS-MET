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


double invRat2(double y)
{
  // UNDER CONSTRUCTION. The output from Wolfram Alpha is a mess. We need to extract repetitive terms 
  // and simplify....



  double a  = cbrt(2.0);  // 2^(1/3)
  double b  = sqrt(3.0);
  double y2 = y*y;
  double y3 = y*y2;
  double y4 = y2*y2;

  double k1  = (128*y3 + 3*b* sqrt(256*y4 + 27*y2) + 27*y);
  double cr1 = cbrt(k1);
  double sr1 = sqrt((16*a*y)/cr1 + cr1/(a*y) + 4);


  //return 0;  // preliminary

  // x1 - Nan:
  // return sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)/(2*b) - 1./2 * sqrt(-(16*a*y)/(3*cr1) - cr1/(3*a*y) - (2*b)/(sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)*y) + 8./3.);

  // x2 - Nan:
  //return sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)/(2*b) + 1./2 * sqrt(-(16*a*y)/(3*cr1) - cr1/(3*a*y) - (2*b)/(sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)* y) + 8./3);

  // x3 - looks wrong:
  //return -sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)/(2*b) - 1./2 * sqrt(-(16*a*y)/(3*cr1) - cr1/(3*a*y) + (2*b)/(sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)*y) + 8./3);

  // x4 - seems to have wrong sign but look otherwise ok:
  return -(1./2 * sqrt(-(16*a*y)/(3*cr1) - cr1/(3*a*y) + (2*b)/(sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)*y) + 8./3) - sqrt((16*a*y)/cr1 + cr1/(a*y) + 4)/(2*b));
  // OK - I just flipped the sign - now it looks ok

  // x1 =  sqrt((16 a y)/cr1 + cr1/(a y) + 4)/(2 b) - 1/2 sqrt(-(16 a y)/(3 cr1) - cr1/(3 a y) - (2 b)/(sqrt((16 a y)/cr1 + cr1/(a y) + 4) y) + 8/3) 
  // x2 =  sqrt((16 a y)/cr1 + cr1/(a y) + 4)/(2 b) + 1/2 sqrt(-(16 a y)/(3 cr1) - cr1/(3 a y) - (2 b)/(sqrt((16 a y)/cr1 + cr1/(a y) + 4) y) + 8/3) 
  // x3 = -sqrt((16 a y)/cr1 + cr1/(a y) + 4)/(2 b) - 1/2 sqrt(-(16 a y)/(3 cr1) - cr1/(3 a y) + (2 b)/(sqrt((16 a y)/cr1 + cr1/(a y) + 4) y) + 8/3)
  // x4 =  1/2 sqrt(-(16 a y)/(3 cr1) - cr1/(3 a y) + (2 b)/(sqrt((16 a y)/cr1 + cr1/(a y) + 4) y) + 8/3) - sqrt((16 a y)/cr1 + cr1/(a y) + 4)/(2 b) 


  //double r1 = 


  


  // x1 =  sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4)/(2 sqrt(3)) - 1/2 sqrt(-(16 2^(1/3) y)/(3 (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)) - (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(3 2^(1/3) y) - (2 sqrt(3))/(sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4) y) + 8/3) and y!=0
  // x2 =  sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4)/(2 sqrt(3)) + 1/2 sqrt(-(16 2^(1/3) y)/(3 (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)) - (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(3 2^(1/3) y) - (2 sqrt(3))/(sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4) y) + 8/3) and y!=0
  // x3 = -sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4)/(2 sqrt(3)) - 1/2 sqrt(-(16 2^(1/3) y)/(3 (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)) - (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(3 2^(1/3) y) + (2 sqrt(3))/(sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4) y) + 8/3) and y!=0
  // x4 =  1/2 sqrt(-(16 2^(1/3) y)/(3 (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)) - (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(3 2^(1/3) y) + (2 sqrt(3))/(sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4) y) + 8/3) - sqrt((16 2^(1/3) y)/(128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3) + (128 y^3 + 3 sqrt(3) sqrt(256 y^4 + 27 y^2) + 27 y)^(1/3)/(2^(1/3) y) + 4)/(2 sqrt(3)) and y!=0


  // See:
  // https://www.wolframalpha.com/input?i=solve+y+%3D+-x+%2F+%28%28x-1%29%5E2*%28x%2B1%29%5E2%29+for+x

}

void sigmoidPrototypes()
{
  // We plot various normalized prototype sigmoids that can be used as basic building blocks for
  // saturators, soft-clippers, etc.

  static const int N = 1001;
  double xMin = 0.0;
  double xMax = 5.0;
  double t    = 0.5;    // threshold below which the function is linear (for the softclipper)

  double x[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);

  int n;
  double yHard[N], yCubic[N], yQuartic[N], yHexic[N], ySoft[N], yTanh[N], yInvRat[N], yInvRat2[N];

  using NS = rsNormalizedSigmoids<double>;

  for(n = 0; n < N; n++) yHard[n]    = NS::clip(x[n]);
  for(n = 0; n < N; n++) yCubic[n]   = NS::cubic(x[n]);
  for(n = 0; n < N; n++) yQuartic[n] = NS::quartic(x[n]);
  for(n = 0; n < N; n++) yHexic[n]   = NS::hexic(x[n]);
  for(n = 0; n < N; n++) ySoft[n]    = rsPositiveSigmoidsD::softClipHexic(x[n], t);
  for(n = 0; n < N; n++) yTanh[n]    = NS::tanh(x[n]);
  for(n = 0; n < N; n++) yInvRat[n]  = rsPositiveSigmoidsD::invRational(x[n]);
  for(n = 0; n < N; n++) yInvRat2[n] = invRat2(x[n]);

  // Using NS:: for ySoft doesn't compile - why? ...Ahh - the function with parameter t is only 
  // available in rsPositiveSigmoidsD but not in rsNormalizedSigmoids

  GNUPlotter plt;
  //plt.addDataArrays(N, x, yHard);
  //plt.addDataArrays(N, x, yCubic);
  //plt.addDataArrays(N, x, yQuartic);
  //plt.addDataArrays(N, x, yHexic);
  //plt.addDataArrays(N, x, ySoft);
  //plt.addDataArrays(N, x, yTanh);
  plt.addDataArrays(N, x, yInvRat);
  plt.addDataArrays(N, x, yInvRat2);
  plt.plot();

  // ToDo:
  //
  // - Add symmetrized versions for all the functions in rsPositiveSigmoids and the use them here.
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

void hilbertDistortion1()
{
  // This experiment implements an idea posted here:
  // https://www.kvraudio.com/forum/viewtopic.php?t=608320
  //
  // We apply waveshaping distortion to the magnitude of an analytic signal and are interested in
  // what that does to the original signal. Given a signal x[n], we first obtain the imaginary part
  // y[n] of the analytic complex signal z[n] = x[n] + i*y[n] by means of a Hilbert transform 
  // filter, then compute the magnitude of this complex signal and pass that into a waveshaper.
  // The output of that waveshaper represents our distorted magnitude. To impose the new magnitude
  // onto our signal z[n], we have to multiply it by the ratio of the distorted magnitude and the
  // original magnitude.

  using WT = RAPT::rsWindowFunction::WindowType;  // Shorthand for convenience

  // Setup:
  int numTaps    = 255;          // Number of taps for Hilbert filter.
  WT  window     = WT::blackman; // Window function for Hilbert filter
  int sampleRate = 44100;        // Sample rate in Hz
  int numSamples = 2000;         // Number of samples to render
  double drive   =  4.0;         // Drive for tanh-waveshaper as raw amplitude multiplier
  double comp    =  1.0;         // Compression amount. 1: normal, 0: none, -1: expand
  bool   makeUp  = false;        // true: Divide by drive again after computing tanh(drive * mag)
  bool   smooth  = false;        // Apply smoothing to Hilbert filter output

  // Processing:

  // Design the Hilbert filter:
  using WFD = rsWindowedFilterDesigner;
  using Vec = std::vector<double>;
  Vec h(numTaps);
  WFD::hilbert(&h[0], numTaps, window);

  // Generate input signal:
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 1, 441.0, double(sampleRate)); // 0: sine, 1: saw, 2: square, 3: triang
  //createModalPluck(&x[0], N, 69.0, sampleRate);  // key = 69 = A4 = 440 Hz

  // Apply a bell-shaped (Hanning window) amplitude-envelope:
  Vec env(N);
  RAPT::rsWindowFunction::createWindow(&env[0], N, WT::hanningZZ, false);
  for(int n = 0; n < N; n++)
    x[n] *= env[n];

  // Obtain Hilbert transform:
  using AT = rsArrayTools;
  Vec y(N+numTaps-1);                             // Convolution result length is M+N-1
  AT::convolve(&x[0], N, &h[0], numTaps, &y[0]);  // Convolution with Hilbert impulse response
  if(smooth)
  {
    if(rsIsOdd(numTaps))
      AT::weightedAverage3pt(&y[0], N+numTaps-1, &y[0], 0.25, 0.5, 0.25);
    else
      AT::movingAverage2ptBackward(&y[0], N+numTaps-1, &y[0]);
  }
  AT::shift(&y[0], N+numTaps-1, -numTaps/2);      // Compensate delay
  y.resize(N);                                    // Shorten y to original length of x

  // Obtain magnitude of analytic signal:
  Vec mag(N);
  for(int n = 0; n < N; n++)
    mag[n] = sqrt(x[n]*x[n] + y[n]*y[n]);

  // Apply distortion:
  Vec magD(N), scl(N), xD(N), yD(N);
  for(int n = 0; n < N; n++)
  {
    // Distort magnitude:
    magD[n] = tanh(drive * mag[n]);
    if(makeUp)
      magD[n] /= drive;              // This is what A_SN does

    // Scale x and y according to ratio of original and distorted magnitude:
    double scaler = magD[n] / mag[n];
    scaler = pow(scaler, comp);
    scl[n] = scaler;
    xD[n]  = scaler * x[n];
    yD[n]  = scaler * y[n];
  }


  // Visualization:
  rsPlotVectors(x, y);         // Input signal and its Hilbert transform
  rsPlotVectors(x, mag);       // Input and instantaneous envelope
  //rsPlotVectors(mag, magD);    // Magnitude and distorted magnitude
  //rsPlotVectors(xD, yD);       // Distorted real and imaginary part
  //rsPlotVectors(x, xD);        // Original and distorted signal
  //rsPlotVectors(scl);          // The scaling factor, i.e. the applied gain
  rsPlotVectors(x, xD, scl);   // Original, distorted signal and applied scaler


  // Observations:
  // -For a fundamental of 440 Hz, 127 taps are not enough. The Hilbert transform of a sinewave at 
  //  440 Hz will be too quiet with 127 taps. 255 seems to be enough. For lower fundamentals, we'll
  //  probably need more.
  // -For a sinewave of constant amplitude of 1 with drive of 4, the original and distorted 
  //  magnitudes look the same - except for artifacts at the ends and some wiggles. It's the
  //  original magnitude that wiggles - I guess, the distortion smoothes/saturates them out.
  // -For the sawtooth and a drive of 2, I get a result similar as shown in the KVR thread. But 
  //  there is some ripple at the Nyquist freq going on. What is this? Is this due to the non-ideal
  //  Hilbert filter? Increasing the length of the filter doesn't seem to help. Maybe it's because
  //  of the bandpass characteristic of the Hilbert filter. Or maybe it is because we need a 
  //  half-sample delay somewhere? But I don't think so - my odd length Hilbert filter should not 
  //  need that.
  // -The choice window function does not seem to have impact on these ripples. 
  // -I tried to use nicer sawtooth frequency (441) to make the aliasing align with the harmonics
  //  but that also doesn't help against the ripples.
  // -The nature of the ripples is: when the saw has no envelope, then the signal y[n] is constant
  //  over every two samples. Sample y[n+1] == y[n]. At the even samples, y jumps to a new value 
  //  and at the odd samples, the value is just repeated. Maybe it has to do with the fact that 
  //  every other sample in h[n] is zero? But that doesn't really explain it because the 
  //  convolution sum is taken over a different set of samples.
  // -Maybe it's because the saw has a component at the Nyquist freq which gets filtered out due
  //  to the bandpass nature?
  // -Results for a sinusoid with bell enevlope (makeUp = false):
  //  -When comp = 1, the drive parameter gives the maximum amount of amplitude boost, i.e. the 
  //   boost that signals close to zero get. For example, for drive = 4, the gain for a zero signal
  //   will be 4.
  //  -For drive = 1.5 and comp = 1.0, the signal gets actually attenuated in the middle. For
  //   drive = 0.5 (comp=1) we get an overall attenuation by factor 0.5 that doesn't depend much
  //   on the input amplitude.
  //  -For comp = -1, drive = 4, the quiet signals are attenuated by a factor of 0.25.
  //  -For comp =  1, drive = 1, we get unit gain for (near) zero signals and an attenuation for 
  //   louder signals.
  // -Results for a sinusoid with bell enevlope (makeUp = true):
  //  -When makeUp = true, we get regular compression rather than upward compression. The quiet 
  //   signals are multiplied by 1 and the loud (unit amplitude) signals are multiplied by 1/drive
  //   when drive > 1. When drive < 1 (like 0.5), not much is happening.
  // -Results for the sawtooth with bell envelope (makeUp = false):
  //  -The general amplitude expansion or compression curve that we see for sines has a 
  //   superimposed comb shape. It drives the amplitude down around the jumps of the saw. This is 
  //   what creates the smoothing/rounding effect. The amplitude goes down and thereby rounds off 
  //   the sawtooth shape. This combing is stronger for louder saws. When the saw is quiet, there 
  //   is not much rounding going on. Quiet saws are expanded, loud saws are smoothed.
  //  -For drive = 0.5, comp = 1.0, there is an overall attenuation of 0.5 and an additional (comb)
  //   smoothing for the louder section.
  // -The best results are obtained by using even Hilbert filter lengths with smoothing. The 
  //  smoothing in case of even length has actually the purpose to introduce the required 
  //  half-sample delay. We don't really need it to smooth out any ripple artifacts because the
  //  even lengths do not produce such artifacts. But if we leave the smoothing out, the alignment
  //  is off by half a sample.
  // -In case of odd filter lengths, we get ugly stairstep artifacts without smoothing. For 
  //  production code, an odd length may be preferable because half of its coeffs are zero which 
  //  allows for a nice optimization. However - even length might potentially give better quality
  //  because of its highpass (rather than bandpass) characteristic. The kernel used for smoothing
  //  is [0.25, 0.5, 0.25]. This produces a zero at the Nyquist freq. This zero goes on top of the
  //  one that the Hilbert filter itself already has.
  //
  // Conclusions:
  // -Values for drive < 1 seem to be not so interesting for sine inputs. We just get an overall 
  //  attenuation.
  // -For saw inputs, values of drive < 1 (like 0.5) give us attenutaion and smoothing. Maybe the
  //  attenuation should be compesated by a compensation gain to get only the smoothing effect?
  // -For drive > 1 and comp > 0, we get an upward compression effect on sinusoids. For comp < 0,
  //  we get an expansion effect.
  // 
  // ToDo:
  // -Try using an odd length filter and an MA filter with weights [0.25 0.5 0.25]
  // -Compare the effects on signals with different overall amplitudes. Use a signal with 
  //  amplitude 1, 2 and 0.5, apply the effect and plot the outputs with gain compensation, i.e.
  //  make the loud output quieter etc. Does an overall gain make a difference? I think so. Then
  //  we should perhaps include a pre-gain parameter. Maybe rename "drive" to "magDrive" to 
  //  distinguish it.
  // -Try some highpass Hilbert-filter design and see if this removes the Nyquist ripples.
  //  See wikipedia article and hilbertFilter experiment for more resources.
  // -Try scaling with the reciprocal. For a saw, that should make it look highpassish.
  // -Try lowpassing instead of waveshaping as "urosh" suggested in the thread.
  // -Try other (perhaps better) windows for the Hilbert filter
  // -For production code, we need to take care of possible division by zero in the computation of
  //  the scaler. We divide by mag[n]. That would be zero for zero signals.
  // -[Done] Use a Hilbert filter with an even number of taps. That gets rid of the staistep 
  //  artifacts. But it will require a shift of y[n] by a half-integer amount. If we do this shift 
  //  with linear interpolation, we essentially apply a 2-sample MA filter after the Hilbert filter 
  //  which will also supress the Nyquist freq. This MA filter will actually attenuate high 
  //  frequencies quite a lot - so it's like an additional smoother. This may not be such a bad 
  //  thing in the context of envelope detection, though - so try it. Maybe also try to realize the 
  //  half-integer shift with cubic interpolation (Lagrange, Hermite, etc. - maybe also try 
  //  "Elephant" interpolation).
  // -When we use some sort of smoothing anyway then maybe we can also use a 3-point MA on the 
  //  result of odd length Hilbert filter to get rid of the stairsteps. Maybe a kernel of
  //  [0.25 0.5 0.25] could be suitable. Making it causal would add another sample of delay to 
  //  which needs to be compensated for.
  // -Wrap the whole Hilbert filtering business into a class rsHilbertFilter. This should also
  //  include the smoothing (optionally - on by default). But maybe for even lengths, where the
  //  smoothing is not really needed, provide alternative ways to achive the half-sample delay. 
  //  Namely by polynomial interpolation. The function makeHilbertFilter could be turned into 
  //  static member function of that class. Maybe call it computeCoeffs().
  //
  // Notes:
  // -The instantaneous envelope of a signal is defined as env[n] = sqrt(x^2[n] + y^2[n]) which we 
  //  can easily solve for y[n] = sqrt(env^2[n] - x^2[n]). Here is a desmos plot for when x is a 
  //  sawtooth and the envelope is constantly 1 (I use f,g there in place of x,y because desmos
  //  already uses x for the independent variable):
  //  https://www.desmos.com/calculator/fun0mfarbo
  //  The resulting kind of looks like what we get for the Hilbert trafo - but not quite - it's 
  //  upside down and has a DC. I'm still doing something wrong, I guess. -> Figure out.
  //  Well - the upside-down could be explained by considering that the sqrt has actually two 
  //  solutions - so we are free to negate the result. Here, I have fudged the graph to make it 
  //  look more like the result we actually get:
  //  https://www.desmos.com/calculator/9rp8yrpjod
  //  ...what's going here? Where is my mistake? Why do I need this fudging? Maybe it's wrong to 
  //  assume that the instantaneous envelope of the saw is just 1? The instantaneous envelope is 
  //  *defined* to be the magnitude of the analytic signal but I guess that doesn't necessarily 
  //  mean that it is the same thing that we intuitively view as envelope visually.
}

void hilbertDistortion2()
{
  // Like hilbertDistortion1 but here we use the class that wraps the whole algorithm into a 
  // realtime processing friendly form.

  // Setup:
  int numTaps    = 255;      // Number of taps for Hilbert filter.
  int sampleRate = 44100;    // Sample rate in Hz
  int numSamples = 2000;     // Number of samples to render
  double drive   =  4.0;     // Drive for tanh-waveshaper as raw amplitude multiplier
  double comp    =  1.0;     // Compression amount. 1: normal, 0: none, -1: expand

  // Create and set up the Hilbert distortion object:
  rsHilbertDistortion<double, double> dist;
  dist.setHilbertFilterLength(numTaps);
  dist.setDrive(drive);
  dist.setCompression(comp);

  // Generate input signal:
  using Vec = std::vector<double>;
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 1, 441.0, double(sampleRate)); // 0: sine, 1: saw, 2: square, 3: triang

  // Produce and plot output:
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = dist.getSample(x[n]);
  rsPlotVectors(x, y); 
}

void hilbertDistortion()
{
  //hilbertDistortion1();
  hilbertDistortion2();
}