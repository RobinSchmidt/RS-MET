#include "OscillatorExperiments.h"

/** Morph between saw-down (p = -1), triangle (p = 0), saw-up (p = +1). A bit care is to be taken 
for p = -1, because there's a potential division by zero in this case (it actually mostly ends
up in the branch where ther is no such division but i'm not totally sure if there can't be 
conditions where this may happen). So maybe best restrict the range to -0.99...+0.99 */
/*
double triSaw(double x, double p) // p=0: triangle, p=-0.99: saw down, p=+0.99: saw up
{
  x /= 2*PI;
  x  = fmod(x, 1);
  double r2 = 0.25*(p+1);
  if(x < r2)
    return x / r2;
  else if(x > 1-r2) 
    return (x-1)/r2;
  else
    return 2*(r2-x)/(1-2*r2) + 1;
}
*/
// was moved to RAPT now

void triSaw()
{
  double p = -0.5;           // parameter 0: triangle, -1: saw down, +1: saw up
  static const int N = 501;  // number of function values
  double t[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(t, N, 0.0, 20.0); 
  for(int n = 0; n < N; n++)
    y[n] = RAPT::rsTriSaw(t[n], p);

  // regular triangle:
  //for(int n = 0; n < N; n++)
  //{
  //  if(t[n] < 0.25)
  //    y[n] = 4*t[n];
  //  else if(t[n] > 0.75)
  //    y[n] = -4 + 4*t[n];
  //  else
  //    y[n] = 2 - 4*t[n];
  //}

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, t, y);
  plt.plot();

  int dummy = 0;
}

void phaseShapeCoeffsXY4(double x, double y, double s, double *a1, double *a2, double *a3, double *a4)
{
  // f(0)=0, f(1)=1, f'(0)=f'(1), f(x)=y, f'(x)=s
  double sx = s*x;
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x2*x2;
  double k  = 1 / (x2*(2*x4-6*x3+7*x2-4*x+1));
  *a1 = k * (x-1)*x*(-s*x2+sx+2*x4-4*x3+(4*x-2)*y);
  *a3 = k * (-sx*(  x3-2*x +1) - 3*x4+2*x2 + y*( 4*(x3-x )+1))*2;
  *a4 = k * ( sx*(2*x2-3*x +1) + 4*x3-3*x2 - y*( 6*(x2-x )+1));
  *a2 = k * ( sx*(3*x3-4*x2+1) + 9*x4-8*x3 - y*(12*(x3-x2)+1));
}
double phaseShapeXY4(double p, double x, double y, double s)
{
  double a1, a2, a3, a4;
  phaseShapeCoeffsXY4(x, y, s, &a1, &a2, &a3, &a4);
  double p2 = p*p;
  return a1*p + a2*p2 + a3*p2*p + a4*p2*p2;
}

void phaseShapeCoeffsXY4(double x, double y, double *a1, double *a2, double *a3, double *a4)
{
  // f(0)=0, f(1)=1, f'(0)=f'(1), f(x)=y, f'(x)=1 (simplified last constraint)
  double x2  = x*x;
  double x3  = x2*x;
  double x4  = x2*x2;
  double k   = (x-y) / (2*x4*x2-6*x4*x+7*x4-4*x3+x2);
  *a1 = 1 + (-4*x3+6*x2-2*x) * k;
  *a2 =   - (12*(x2-x3)-1)   * k;
  *a3 =       (8*(x-x3)-2)   * k;
  *a4 =   -   (6*(x-x2)-1)   * k;
}
double phaseShapeXY4(double p, double x, double y)
{
  double a1, a2, a3, a4;
  phaseShapeCoeffsXY4(x, y, &a1, &a2, &a3, &a4);
  double p2 = p*p;
  return a1*p + a2*p2 + a3*p2*p + a4*p2*p2;
}

void phaseShapingCurvePoly4()
{
  // Plots the 4th order polynomial phase shaping curve for a given set of parameters.

  double x = 0.4; // x for xy-shaper
  double y = 0.8; // y for xy-shaper
  double s = y/x; // slope for xy-shaper - can be assigned arbitrary but y/x kind of makes sense

  s = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)); 
    // Euclidean distance of (x,y) from (0.5,0.5)
    // satisfies f(x,y) = f(1-x,1-y) i think - provides for a mirrored curve when we apply
    // the transformations x -> 1-x, y -> 1-y to x,y
    // but for x=y=0.5 it gives zero slope at (x,y) - neutral would be unit slope

  s = y/x; // gives neutral curve when x=y=0.5
  s = 1;
    // maybe form a weighted sum of sqrt(..) and y/x the user decides, how much of which function
    // he wants with a parameter p. p=0: s=y/x, p=1: s=sqrt(..)

  static const int N = 1001;
  double p[N];    // unshaped phase
  double ps[N];   // shaped phase
  RAPT::rsArrayTools::fillWithRangeLinear(p, N, 0.0, 1.0);
  for(int n = 0; n < N; n++)
    ps[n] = phaseShapeXY4(p[n], x, y, s);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, p, p, ps);
  plt.plot();

  // Observations:
  // mirror curves: x=0.4,y=0.8,s=1 / x=0.6,y=0.2,s=1
  // -> applying the transformation x -> 1-x, y -> 1-y, s -> s gives mirror curves
  // maybe find a formula for s as function of x,y such that f(x,y) = s = f(1-x, 1-y)
}

double rsRationalMap(double x, double a)
{
  // "a" is a parameter from -1..+1 (ends excluded)
  return (x-a)/(1-a*x);
}
double phaseShapeRational(double p, double a)
{
  double y;
  if(p <= 0.5)
  {
    p = rsLinToLin(p, 0.0, 0.5, -1.0, +1.0);  // optimize...
    y = rsRationalMap(p, a);
    y = rsLinToLin(y, -1.0, +1.0, 0.0, 0.5);
  }
  else
  {
    p = rsLinToLin(p, 0.5, 1.0, -1.0, +1.0);  // optimize...
    y = rsRationalMap(p, -a);
    y = rsLinToLin(y, -1.0, +1.0, 0.5, 1.0);
  }
  return y;
}

void phaseShapingCurvesRational()
{
  // Plots a family of rational phaseshaping curves for various value of the parameter 

  static const int N = 201;        // number of function values
  static const int numCurves = 7;  // number of curves in the function family on each side

  double x[N], y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.0);

  GNUPlotter plt;
  plt.addDataArrays(N, x, x);

  double a1 = 1.0 / (numCurves+1);   // 1st parameter value
  double ai;                          // i-th parameter value
  int i, n;
  for(i = 1; i <= numCurves; i++)
  {
    ai = i * a1;
    for(n = 0; n < N; n++)
      y[n] = phaseShapeRational(x[n], ai);
    plt.addDataArrays(N, x, y);
    for(n = 0; n < N; n++)
      y[n] = phaseShapeRational(x[n], -ai);
    plt.addDataArrays(N, x, y);
  }
  plt.plot();

  // ToDo:
  // -Use shaping function y(x) = (a + b*x) / (c + dx) and require: y(0)=0, y(1)=1, y'(0)=y'(1)
  //  (1) leads to a=0, then (2) leads to b= c+d
}

vector<double> createPhase(int N, double f, double fs)
{
  // creates a phase-sawtooth of length N for frequency f and samplerate fs
  vector<double> p(N);
  createWaveform(&p[0], N, 1, f, fs, 0.0, false);  // sawtooth
  //int n;
  for(int n = 0; n < N; n++)
    p[n] = 0.5*(p[n]+1);                           // normalize to 0..1
  return p;
}

void phaseShaping()
{
  // We consider a normalized phase-function p(t) that would normally a periodic upward sawtooth 
  // function between 0..1 with a period length of 1. A sinusoidal waveform could the be created 
  // from this phase function as y(t) = sin(2*pi*p(t)). But instead of using p(t) as is, we apply a 
  // 4th order polynomial to it. The polynomial satisfies 
  // f(0)=0, f(1)=1, f'(0)=f'(1), f(x)=y, f'(x)=s where x,y,s are user parameters.
  
  // We create a few standard-waveforms (sine, triangle, saw, square), sweeping the midpoint 
  // linearly between some user defined values. ..currently only sine

  double fs = 44100;        // samplerate
  double f  = 100;          // frequency
  //int N = 80000;
  int N       = (int)(5*fs);  // number of samples
  double xMin = -1.0;         // minimum x-value
  double xMax = +2.0;         // maximum x-value
  double yMin =  0.5;         // minimum y-value
  double yMax =  0.5;         // maximum y-value
  double sMin =  1.0;         // minimum slope
  double sMax =  1.0;         // maximum slope
  double amp  =  0.5;         // amplitude

  // create the normalized phase-function:
  //vector<double> p(N);
  //createWaveform(&p[0], N, 1, f, fs, 0.0, false);  // sawtooth
  //int n;
  //for(n = 0; n < N; n++)
  //  p[n] = 0.5*(p[n]+1);                           // normalize to 0..1

  vector<double> p = createPhase(N, f, fs);

  // create the x/y sweep functions:
  vector<double> x(N), y(N), s(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  RAPT::rsArrayTools::fillWithRangeLinear(&y[0], N, yMin, yMax);
  RAPT::rsArrayTools::fillWithRangeLinear(&s[0], N, sMin, sMax);

  // apply the the phase-shaping:
  int n;
  vector<double> ps(N);
  for(n = 0; n < N; n++)
    ps[n] = phaseShapeXY4(p[n], x[n], y[n], s[n]);

  // create a phase-shaped sine wave:
  vector<double> ySin(N);
  for(n = 0; n < N; n++)
    ySin[n] = amp * sin(2*PI*ps[n]);

  // write output wavefile
  writeToMonoWaveFile("PhaseShaping.wav", &ySin[0], N, (int) fs, 16);

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, &p[0]);
  //plt.addDataArrays(N, &ps[0]);
  //plt.addDataArrays(N, &ySin[0]);
  //plt.addDataArrays(N, &m[0]);
  //plt.plot();

  int dummy = 0;

  // Ideas:
  //
  // Before the polynomial, apply a (appropriately symmetrized) Laguerre-type warping map to the 
  // phase. Depensing on the warp-parameter, this should introduce an effect that creates the
  // skewed sine (sine skewed into saw).
  //
  // Let x(t), y(t), s(t) be decaying sinusoids (scaled and shifted such that x,y center at 0.5,
  // s centers at 1.0)
  //
  // Instead of setting s directly, let it be a function of x,y for example:
  // s = f(x,y) = y/x or s = f(x,y) = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)). The second one
  // is the Euclidean distance of (x,y) from (0.5,0.5). When using this formula, it is ensured that
  // f(x,y) = f(1-x,1-y). This provides for a mirrored curve when we apply the transformations 
  // x -> 1-x, y -> 1-y to x,y, but for x=y=0.5 it gives zero slope at (x,y) - neutral would be 
  // unit slope. It's also possible to use a weighted sum of sqrt(..) and y/x where the relative
  // weights are controlled by a user parameter.
  //
  // Instead of a single 4th order polynomial, use two 2nd order polynomials for a piecewise 
  // defined function f(p) = f1(p) for p <= x and f(p) = f2(p) for p > x
  // require: f1(0)=0, f1(x)=y, f2(x)=y, f2(1)=1, f1'(0)=f2'(1), f1'(x)=f2'(x) with
  // f1(x) = a0 + a1*x + a2*x^2, f2(x) = b0 + b1*x + b2*x^2
  //
  // Maybe compare to phase-distortion synthesis

  // todo: have 2 additional control-points for the new 1/4 and 3/4 cycle-point
  // maybe require the 2nd derivative to be zero at 0 and 1
}

void phaseShapingSkew()
{
  double fs   = 44100;        // samplerate
  double f    = 100;          // frequency
  int N       = (int)(3*fs);  // number of samples
  double aMin = -0.99;        // minimum a-value
  double aMax = +0.99;        // maximum a-value
  double amp  =  0.5;         // amplitude

  vector<double> p = createPhase(N, f, fs);

  // create the parameter sweep function:
  vector<double> a(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&a[0], N, aMin, aMax);

  // apply the the phase-shaping:
  int n;
  vector<double> ps(N);
  for(n = 0; n < N; n++)
    ps[n] = phaseShapeRational(p[n], a[n]);

  // create a phase-shaped sine wave:
  vector<double> ySin(N);
  for(n = 0; n < N; n++)
    ySin[n] = amp * sin(2*PI*ps[n]);

  // write output wavefile
  writeToMonoWaveFile("PhaseShapeSkewSine.wav", &ySin[0], N, (int) fs, 16);


  // todo: maybe experiment with different shaping functions (for example, Laguerre map)

  int dummy = 0;
}

void zeroDelayFeedbackPhaseMod()
{
  // We want to generate a waveform via feedback phase-modulation without delay in the feedback 
  // path. The equation is 
  //   y[n] = sin(p[n] + a*y[n])  
  // where p[n] = w*n is the regular sine phase. We can't solve explicitly for y[n], so we consider
  // it as a root finding problem for the equation f(y) = y - sin(p + a*y) = 0. We solve it via 
  // Newton iteration. The derivative is given by f'(y) = 1 - a*cos(p + a*y)

  int N = 6000;         // number of samples to render
  double cycle = 2000;  // length of one cycle of the sinewave

  double a = 0.996;     
  // feedback amount, 0.996 works for cycle=2000 with Newton, 0.997 doesn't

  double tol = 1.e-13;
  int maxIts = 100;     // maximum number of Newton steps, if 0, we get unit-delay feedback PM

  double w = 2*PI/cycle;
  std::vector<double> y(N), z(N), q(N);

  double p;

  // Define objective function for root-finder:
  auto func = [&](double y) { return y - sin(p + a*y); };


  for(int n = 1; n < N; n++)
  {
    p = w*n;                           // sine phase
    double t = sin(p + a*y[n-1]);      // initial guess using UDF-FB-PM, may be better to use z
    y[n] = t;                          // y is the UDF signal
    
    // Newton iteration:
    //t = 0;  // test - use 0 as initial guess - nope!
    for(int i = 0; i < maxIts; i++)    // Newton iteration
    {
      double f  = t -   sin(p + a*t);  // function whose root we want to find
      double fp = 1 - a*cos(p + a*t);  // derivative of our function
      double d  = -f/fp;               // calculate Newton step
      t += d;                          // perform Newton step
      if(abs(d) <= tol) break;         // break if converged
    }
    
    // Test to use bisection instead of Newton:
    //t = rsRootFinder<double>::bisection(func, -1.0, 1.0);

    z[n] = t;                          // z is the ZDF signal
    q[n] = asin(t);                    // the new phase angle
  }

  // Observations
  // -For a = 1.5, we get some discontinuities. Is that correct? When we set maxIts=0, we get
  //  garbage, so probably the Newton iteration does not converge to the correct value.
  // -The greatest difference betwen UDF and ZDF exists around the discontinuity of the saw. There
  //  the UDF version tends to get asymmteric whereas the "cleaned up" ZDF version is still nicely
  //  symmetric.
  // -For cycle=2000, N=6000, a=0.996 it still works but for 0.997 it doesn't. But increasing cycle
  //  to 3000 and N to 9000, it works again. Conclusion: convergence of Newton iteration depends
  //  on feedback factor (higher -> unstable) and cycle-length (shorter -> unstable)
  // -We get a slightly larger stable range when using bisection instead of Newton, but it also 
  //  starts to behave weird from 1 upwards
  //
  // Ideas:
  // -How about apply a saturation like tanh to the phase where the phase in in -pi..+pi and the
  //  output of the tanh is also scaled into that range. The goal is to get a smooth transition
  //  between sine and saw. When applying tanh to the phase (with proper scaling of input and
  //  output), would we get a smooth wrap-around? I think so but I'm not sure. It would be 
  //  useful, if we could achieve the goal with an explicit feedforward algorithm. Feedback PM
  //  works nicely but is computationally expensive. Maybe implement a few variants of a 
  //  sineToSaw function taking as input the phase p and maybe an initial guess
  // -I suppose, the problem with higher feedback is that the initial guess for the iteration 
  //  becomes more and more wrong around the signal's discontinuity. Maybe we need a different
  //  strategy to come up with an initial guess. Having to pass in the previous value is 
  //  inconvenient anyway. Maybe just take a triangular approximation of the sine that squeezes
  //  into a saw with higher feedback. And/or maybe try bisection instead of Newton iteration. Or
  //  Brent's method or regula falsi. Or: when the phase is close to a multiple to 2*pi, drag the
  //  initial guess toward zero
  // -Maybe try to use peridoci cubic spline interpolation of the phase: fix points (-pi,-pi),
  //  (+pi,+pi) and for the middle point (0,y) make y adjustable. The result would be only 2nd 
  //  order smooth but maybe that's good enough. After all, we actually want to create some 
  //  harmonics with this technique anyway.

  // Try to unwrap the phase:
  //rsArrayTools::add(&q[0], +2*PI, &q[0], N);
  //rsArrayTools::unwrap(&q[0], N, 2*PI);
  // ...that doesn't seem to work yet. ah - because the unwrapping assumes discontinuities in the
  // function itself but we have only discontinuities in the derivative. Maybe try to
  // differentiate -> unwrap -> integrate (numerically)? But wrap at what value?

  //rsPlotVectors(y);  // only UDF signal
  rsPlotVectors(y, z, z-y, q);
  int dummy = 0;


  // Resources:
  // https://ristoid.net/modular/fm_variants.html
  // https://ccrma.stanford.edu/software/snd/snd/fm.html
}