#include "MathExperiments.h"
using namespace RAPT;

void ellipseLineIntersections()
{
  // create and set up ellipse:
  rsEllipseF ellipse;
  ellipse.setParameters(1.5f, 3.f, float(PI/4), 0.5f, 1.f); 

  // create line parameters:
  float x, y, dx, dy;
  x  = 0.5f;
  y  = 0.5f;
  dx = 0.5f;
  dy = 0.2f;

  // create data for drawing the ellipse:
  static const int Ne = 500;
  float xe[Ne], ye[Ne];
  float err;
  for(int n = 0; n < Ne; n++)
  {
    float phi = float(2*PI*n) / (Ne-1);
    ellipse.getPointOnEllipse(phi, &xe[n], &ye[n]);
    err = ellipse.evaluate(xe[n], ye[n]);
  }

  // create data for drawing the line:
  float tMin = -3.5;
  float tMax = +3.0;
  float xl[2], yl[2];
  xl[0] = x + tMin*dx;
  xl[1] = x + tMax*dx;
  yl[0] = y + tMin*dy;
  yl[1] = y + tMax*dy;

  // find intersection points between line and ellipse:
  float ti1, ti2, xi1, yi1, xi2, yi2;
  ellipse.lineIntersectionParameter(x, dx, y, dy, &ti1, &ti2);
  xi1 = x + ti1*dx;
  yi1 = y + ti1*dy;
  xi2 = x + ti2*dx;
  yi2 = y + ti2*dy;

  // find tangent line to intersection point:
  float A, B, C, a, b;
  float xt1[2], yt1[2], xt2[2], yt2[2];
  ellipse.getTangentCoeffs(xi1, yi1, &A, &B, &C); // implicit  A*x + B*y + C = 0
  a = -A/B;     // explicit y = a*x + b
  b = -C/B;
  xt1[0] = -1.2f;   // use points where x=0, x=2 to draw the tangent:
  yt1[0] = a*xt1[0] + b;   
  xt1[1] = -0.5f;
  yt1[1] = a*xt1[1] + b; 

  // same for the 2nd intersection:
  ellipse.getTangentCoeffs(xi2, yi2, &A, &B, &C);
  a = -A/B;
  b = -C/B;
  xt2[0] = 0;
  yt2[0] = a*xt2[0] + b;   
  xt2[1] = 2;
  yt2[1] = a*xt2[1] + b; 

  GNUPlotter plt;
  plt.setRange(-1.5, +2.5, -1, +3);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");   // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(Ne, xe,  ye);   // ellipse
  plt.addDataArrays(2,  xl,  yl);   // line
  plt.addDataArrays(2,  xt1, yt1);  // tangent at 1st intersection
  plt.addDataArrays(2,  xt2, yt2);  // tangent at 2nd intersection
  plt.plot();
}

void linearRegression()
{
  static const int N = 500;  // number of samples
  float minDist = 0.1f;       // minimum distance between successive samples
  float maxDist = 1.0f;      // maximum ...
  float a       = 0.8f;
  float b       = 20.0f;
  float noise   = 10.0f;
  int   seed    = 0;

  // create input data:
  float x[N], y[N];
  float dx;      // delta between samples
  x[0] = 10.f;
  y[0] = a*x[0] + b + noise * (float)round(rsRandomUniform(-1, +1, seed)); 
  int n;
  for(n = 1; n < N; n++){
    dx = (float)rsRandomUniform(minDist, maxDist);
    x[n] = x[n-1] + dx;
    y[n] = a*x[n] + b + noise * (float)rsRandomUniform(-1, +1); 
  }

  // retrieve a,b from the data:
  float ar, br = 0;
  rsStatistics::linearRegression(N, x, y, ar, br);
  //ar = Statistics::proportionalRegression(N, x, y);

  // create lines for correct and retrieved line parameters:
  float yc[N], yr[N];
  for(n = 0; n < N; n++){
    yc[n] = a *x[n] + b;
    yr[n] = ar*x[n] + br; }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0], &yc[0], &yr[0]);
  plt.plot();
}




// maybe move to RAPT into the Statistics section
double variance(double *x, int N)
{
  double mx  = RAPT::rsArray::mean(x, N); // mean of x
  double sum = 0;
  for(int n = 0; n < N; n++) {
    double d = x[n] - mx;
    sum += d*d;
  }
  return sum / (N-1);
  // todo: make an optimized version that takes the mean as argument
}
double standardDeviation(double *x, int N)
{
  return sqrt(variance(x, N));
}
double covariance(double *x, double *y, int N)
{
  double mx = RAPT::rsArray::mean(x, N); // mean of x
  double my = RAPT::rsArray::mean(y, N); // mean of y
  double sum = 0;
  for(int n = 0; n < N; n++)
    sum += (x[n]-mx) * (y[n]-my);
  return sum / (N-1);
}
double correlation(double *x, double *y, int N)
{
  double vxy = covariance(x, y, N);
  double vx  = variance(x, N);
  double vy  = variance(y, N);
  return vxy / (vx*vy);
}
double conditionalProbability(double* a, double* b, int N)
{
  // conditional probability of "a given b"
  // a and b must be arrays of boolean values (0 or 1) but as double data type
  int na = 0, nb = 0;
  for(int n = 0; n < N; n++) {
    if(b[n] == 1) {
      nb++;
      if(a[n] == 1)
        na++;
    }
  }
  return (double) na / (double) nb;
}
void probabilityLogic()
{
  // consider an example: produce random number x in 0..1, 

  // let: A:  x < 0.5, B: 0.3 < x < 0.6,
  // -> P(A)=0.5, P(B)=0.3, P(A|B)=2/3, P(B|A)=2/5=0.4 -> what's the correlation C(A,B)?


  // setup:
  int N = 100000;   // number of realizations to produce
  double aL = 0.0;  // lower limit for x, such that event A is considered to have occured
  double aU = 0.5;  // upper limit ....
  double bL = 0.3;  // ...same for event B
  double bU = 0.6;


  typedef std::vector<double> Vec;
  Vec A(N), B(N);   // realizations of events A and B (true/false, represented as 0/1)
  Vec x(N);
  int n;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0, 1);
  for(n = 0; n < N; n++)
  {
    double xn = prng.getSample();
    x[n] = xn;

    if(xn > aL && xn < aU)  
      A[n] = 1;
    else         
      A[n] = 0;

    if(xn > bL && xn < bU)
      B[n] = 1;
    else
      B[n] = 0;
  }

  // compute relative frequencies of events A and B (should approximate their probabilities):
  //double fA = RAPT::rsArray::sum(&A[0], N) / N;
  //double fB = RAPT::rsArray::sum(&B[0], N) / N;
  // are actually the mean values

  // compute sample mean values for event A and B:
  double mA = RAPT::rsArray::mean(&A[0], N);
  double mB = RAPT::rsArray::mean(&B[0], N);

  // compute sample variances:
  double vA = variance(&A[0], N);
  double vB = variance(&B[0], N);

  // compute sample covariance and correlation:
  double cov = covariance( &A[0], &B[0], N);
  double cor = correlation(&A[0], &B[0], N);

  // compute empirical conditional probabilities:
  double cab = conditionalProbability(&A[0], &B[0], N); // P(A|B), empirical prob of A given B
  double cba = conditionalProbability(&B[0], &A[0], N); // P(B|A), empirical prob of B given A



  int dummy = 0;
}

double productLog(const double z) 
{
  // Evaluates the product-log function W(z) defined as the inverse of the function f(W) = W * e^W.
  // This function is also known as Lambert-W or Omega function.

  // adapted from http://keithbriggs.info/software/LambertW.c
  // here's more: http://keithbriggs.info/software.html

  if(0.0 == z) 
    return 0.0;

  double eps = 4.0e-16;
  double em1 = 0.3678794411714423215955237701614608; // 1/Euler
  double p, e, t, w;

  if(z < -em1+1e-4) // series near -em1 in sqrt(q)
  { 
    double q = z+em1, r = sqrt(q), q2 = q*q, q3 = q2*q;
    return 
      -1.0
      +2.331643981597124203363536062168*r
      -1.812187885639363490240191647568*q
      +1.936631114492359755363277457668*r*q
      -2.353551201881614516821543561516*q2
      +3.066858901050631912893148922704*r*q2
      -4.175335600258177138854984177460*q3
      +5.858023729874774148815053846119*r*q3
      -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }

  // initial approximation for Halley iteration:
  if(z < 1.0)     // series near 0
  { 
    p = sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));                 // euler-number
    w = -1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); // -1/3, 11/72
  } 
  else 
    w = log(z);   // asymptotic

  if(z > 3.0) 
    w -= log(w);  // useful?

  // Halley iteration:
  for(int i = 0; i < 10; i++) 
  { 
    e  = exp(w); 
    t  = w*e-z;
    p  = w+1.0;
    t /= e*p-0.5*(p+1.0)*t/p; 
    w -= t;
    if(fabs(t) < eps*(1.0+fabs(w)))  // rel-abs error
      return w;
  }

  //rsAssertFalse(); // no convergence
  return 0.0;
}
void productLogPlot()
{
  int N = 1000;        // number of values
  float xMin = -3.0;
  float xMax = +15.0;

  // evaluate function:
  vector<float> x(N), y(N);
  rsArray::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
    y[n] = (float) productLog(x[n]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}

void sinCosTable()
{
  // A test for the rsSinCosTable class.

  rsSinCosTableF table(8); // parameter is the table size 

  // create data:
  int N = 2000;  // number of values to plot
  float xMin = -15.0;
  //float xMin =   0.0;
  float xMax = +15.0;
  vector<float> x(N), ySin(N), yCos(N), ySinTbl(N), yCosTbl(N);
  rsArray::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    ySin[n] = sin(x[n]);
    yCos[n] = cos(x[n]);
    //table.getValuesRounded(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesTruncated(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesLinear(x[n], &ySinTbl[n], &yCosTbl[n]);
    table.getValuesCubic(x[n], &ySinTbl[n], &yCosTbl[n]);
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &ySin[0], &yCos[0], &ySinTbl[0], &yCosTbl[0]);
  plt.plot();
}

void expBipolar()
{
  // find coeffs for a parametrized exponential function y = f(x) = a * e^(b*x) + c such that 
  // f(0) = y0, f(1) = y1, f'(0) = s. We have f'(x) = a*b*e^(b*x). The equation system is:
  // y0 = a + c         
  // y1 = a * e^b + c
  // s  = a * b
  // solving 1 and 3 for c and b and plugging into 2 gives:
  // d = y1 - y0 == a e^(s/a) - a
  // which can be given to wolfram alpha as:
  // Solve[d == a e^(s/a) - a, a ]

  double y0 = 0.2;
  double y1 = 8.0;
  double s  = 1.0;

  // compute coeffs:
  //double dy  = y0 - y1;
  //double tmp = (exp(s/dy)*s)/dy;
  //tmp = productLog(tmp);
  //double a = s*dy / (-dy * tmp + s); // nope - formula must be wrong 

  double d = y1 - y0;
  double w = productLog(- (exp(-s/d)*s) / d);
  double a = -d*s / (d*w+s);
  double b = s  / a;
  double c = y0 - a;

  // compute values at the endpoints for test:
  double f0  = a * exp(b * 0) + c; // f0 is wrong
  double f1  = a * exp(b * 1) + c;
  double fp0 = a * b;

  int dummy = 0;
}

void expGaussBell()
{
  // This is still wrong. The idea is to find a,b,c,d parameters for the function
  // f(x) = a * exp(b*x) + c * exp(d*x^2)
  // and impose the constraints f(0) = 1, f'(0) = s0, f(1) = 0, f'(1) = s1.
  // we have:
  // f (x) = a * exp(b*x) +   c  *  exp(d*x^2)
  // f'(x) = a*b*exp(b*x) + 2*c*d*x*exp(d*x^2)
  // so:
  // (1) f (0) = 1  = a + c
  // (2) f'(0) = s0 = a * b
  // (3) f (1) = 0  = a*exp(b) + c*exp(d)
  // (4) f'(1) = s1 = a*b*exp(b) + 2*c*d*exp(d)
  // solve(1=a+c,s_0=a*b,0=a*exp(b)+c*exp(d),s_1=a*b*exp(b)+2*c*d*exp(d);a,b,c,d) -> nope
  // we still need to solve the system - maybe an expression that contains only one unknown 
  // - say a - can be derived and we can st up a Newton iteration to the 1D root finding problem
  // we have: (1) c = (1-a), (2) b = (s0/a), 
  // (3) d = ln(-a*exp(b)/c) = ln(-a*exp(s0/a)/(1-a))
  // (4) s1 = a*b*exp(b) + 2*c*d*exp(d)
  // s1 = s0*exp(s0/a) + 2*(1-a)*d*exp(d)
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a))*exp(ln(-a*exp(s0/a)/(1-a))) - s1
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a)) * -a*exp(s0/a)/(1-a) - s1
  // ...more math to do....


  // settings:
  int N = 300;        // number of values
  float xMin = 0.0;
  float xMax = 1.0;
  double s0   = -0.1f; // slope at 0
  double s1   = -0.2f; // slope at 1

                       // compute coefficients:
  double a, b, c, d;
  a = -s0 / productLog(0.5*s0*s1);
  b = -s0 / a;
  c = 1 - a;
  d = b; // wrong!

         // evaluate function:
  vector<float> x(N), y(N);
  rsArray::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    double t = x[n]; // temporary
    t = a*exp(b*t) + c*exp(d*t);
    y[n] = (float)t;
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}

// fun stuff:

void primeAlternatingSums()
{
  // Take the array of prime numbers with alternating signs
  // 2 -3 5 -7 11 -13 17 -19 ...
  // and take running sums of various order ...just for fun to see what happens
  // ...do they also all have alternating signs? what happens, if we add them to one of the
  // previous arrays...just mess around a little - may interesting patterns emerge...
  // -what changes, if we give the even-indexed primes a negative sign?
  // -what happens, if we flip every k-th sign instead of every 2nd?
  // -what happens, if we filter the arrays with integer-coefficient IIR filters? this is a 
  //  generalization of (iterated) running sums (i think)
  // -what about nonlinear and/or time-variant recursive filters with coeffs depending on index?
  // -what about trying to use non/linear prediction methods on the time series?
  //  -maybe the equations (correlation coeffs, etc) should use rational numbers?
  //  -test on number sequences with known structure first (like fibonacci numbers), see, if they
  //   pick it up correctly
  // -what if we sign invert not every second number but every second pair of numbers? ..and then 
  //  combine the resulting sequences with those obtained from inerting every other number ...then
  //  generalize: invert every triple, quadruple, etc...

  int N = 200; // number of primes

  //typedef std::vector<int> IntVec;
  typedef RAPT::rsArray AR;

  // create array of primes:
  std::vector<int> primes(N);
  RAPT::rsFillPrimeTable(&primes[0], N);  // make convenience function getPrimes

  // create arrays of primes with alternating signs:
  int n;
  std::vector<int> ev, od; // ev: have even-numberd signs flipped, od: odd numbered signs flipped
  ev = primes; od = primes;  
  for(n = 0; n < N; n += 2) ev[n] = -ev[n];
  for(n = 1; n < N; n += 2) od[n] = -od[n];

  // 1st order running sums:
  std::vector<int> ev1(N), od1(N);
  AR::cumulativeSum(&ev[0], &ev1[0], N);
  AR::cumulativeSum(&od[0], &od1[0], N);

  // 2nd order running sums:
  std::vector<int> ev2(N), od2(N);
  AR::cumulativeSum(&ev1[0], &ev2[0], N);
  AR::cumulativeSum(&od1[0], &od2[0], N);

  // 3rd order running sums:
  std::vector<int> ev3(N), od3(N);
  AR::cumulativeSum(&ev2[0], &ev3[0], N);
  AR::cumulativeSum(&od2[0], &od3[0], N);

  // ...at this point, we really should use a loop



  // plot:
  GNUPlotter plt;
  plt.setPixelSize(1000, 500);
  //plt.addDataArrays(N, &primes[0]);
  //plt.addDataArrays(N, &ev[0]);
  plt.addDataArrays(N, &od[0]);
  //plt.addDataArrays(N, &ev1[0]);
  plt.addDataArrays(N, &od1[0]);
  //plt.addDataArrays(N, &ev2[0]);
  plt.addDataArrays(N, &od2[0]);
  //plt.addDataArrays(N, &ev3[0]);
  plt.addDataArrays(N, &od3[0]);
  plt.plot();
  // todo: maybe plot with stems or at least, get rid of the linear interpolation (draw steps)

  // the 3rd order sum is very smooth ... subtract even and odd 3rd order sums...

  // ..but maybe it's actually faster to noodle around with that stuff in python or sage
}