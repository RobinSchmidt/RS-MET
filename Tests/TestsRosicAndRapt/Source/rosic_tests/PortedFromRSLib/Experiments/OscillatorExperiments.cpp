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
  // Oscillator that moprhs between sawDown/triangle/sawUp. With sinusoidal waveshaping we can also
  // make it morph between fatSawDown/sine/fatSawUp where "fatSaw" is a variation of the saw that
  // bulges outward.

  double p = -0.5;           // parameter 0: triangle, -1: saw down, +1: saw up
  static const int N = 501;  // number of function values
  double t[N], y[N], z[N];
  RAPT::rsArrayTools::fillWithRangeLinear(t, N, 0.0, 20.0); 
  for(int n = 0; n < N; n++)
  {
    y[n] = RAPT::rsTriSaw(t[n], p); // TriSaw
    z[n] = sin(0.5*PI * y[n]);      // SinSaw
  }

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
  plt.addDataArrays(N, t, y, z);
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

// Maps the interval x = [-1,+1] to itself monotonically where the curve is determined by the
// parameter a
double rsRationalMap(double x, double a)
{
  // "a" is a parameter from -1..+1 (ends excluded)
  return (x-a)/(1-a*x);
}
// https://www.desmos.com/calculator/imnmitfegz
// Maybe rename! Maybe interpret the a coeff as negative such that the curve gets concave for 
// positive a. This is more intuitive. Or maybe not? Currently the weight of the graph moves to 
// the left with negative values which is intuitive when using a horizontal slider for a. On the 
// other, it also moves up for negative values which is counterintuitive for a vertical slider.

double phaseShapeRational(double p, double a, double b = 1, double c = 1)
{
  double y;
  p = pow(p, b);
  if(p <= 0.5)
  {
    p = rsLinToLin(p, 0.0, 0.5, -1.0, +1.0);  // 0..0.5 -> -1..+1
    y = rsRationalMap(p, a);
    y = rsLinToLin(y, -1.0, +1.0, 0.0, 0.5);  // -1..+1 -> 0..0.5
  }
  else
  {
    p = rsLinToLin(p, 0.5, 1.0, -1.0, +1.0);  // 0.5..1 -> -1..+1
    y = rsRationalMap(p, -a);
    y = rsLinToLin(y, -1.0, +1.0, 0.5, 1.0);  // -1..+1 -> 0.5..1
  }
  y = pow(y, c);
  return y;
}
// This is the same map as now implemented in rsLinearFractionalInterpolator - ah - no - not quite.
// Here we use
// Maybe introduce an exponent applied to p before the switch


//  a = 0.5..2
// get rid
double phaseShapePower(double p, double a)
{
  return rosic::rsPhaseShaper::powerLaw(p, a);

  /*
  p = rsLinToLin(p, 0.0, 1.0, -1.0, +1.0);   // 0..1 -> -1..+1
  if(p >= 0)
    p =  pow( p, a);
  else
    p = -pow(-p, a);
  p = rsLinToLin(p, -1.0, +1.0, 0.0, 1.0);   // -1..+1 -> 0..1
  return p;
  */
}
// needs test


void phaseShapingCurvesRational()
{
  // Plots a family of rational phaseshaping curves for various value of the parameter 

  static const int N = 201;        // number of function values
  static const int numCurves = 7;  // number of curves in the function family on each side

  // some more tweakable parameters:
  double b = 2.0;   // 0.5..2, default 1
  double c = 0.5;   // ....dito........., 1/b may make sense

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
      y[n] = phaseShapeRational(x[n],  ai, b, c);
    plt.addDataArrays(N, x, y);
    for(n = 0; n < N; n++)
      y[n] = phaseShapeRational(x[n], -ai, b, c);
    plt.addDataArrays(N, x, y);
  }
  plt.plot();

  // Observations:
  // -I think for b != 1 and/or c != 1, the slopes at 0 and 1 don't match (verify!)

  // ToDo:
  // -Use shaping function y(x) = (a + b*x) / (c + dx) and require: y(0)=0, y(1)=1, y'(0)=y'(1)
  //  (1) leads to a=0, then (2) leads to b= c+d
  // -Try to apply the exponent at a different place. It doesn't have the effect that I wanted of 
  //  changing the shape symmetrically. I think, it's because the function here is based on mapping
  //  -1..+1 rather than 0..1 to itself
  // -Try using a symmetrized power function
  // -Try the symmetrized LinFrac formula
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
  plt.addDataArrays(N, &p[0]);
  //plt.addDataArrays(N, &ps[0]);
  //plt.addDataArrays(N, &ySin[0]);
  //plt.addDataArrays(N, &m[0]);    // m does not exist
  //plt.plot();
  // Plotting phase does not make sense for a full signal

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
  // Common parameters:
  double fs   = 44100;        // samplerate
  double f    = 100;          // frequency
  int N       = (int)(3*fs);  // number of samples
  double amp  =  0.5;         // amplitude

  // Declare variables:
  using Vec = std::vector<double>;
  double aMin;                      // minimum a-value (defined per shaping function)
  double aMax;                      // maximum a-value (defined per shaping function)
  Vec p = createPhase(N, f, fs);    // Raw phase (as input for computing shaped phase)
  Vec ps(N);                        // Shaped phase (for optional plotting)
  Vec y0(N), y1(N), y2(N), y3(N);   // Output signals (for wav-writing)

  // Generate the signal with rational phase-shaping map:
  aMin = -0.99;        // minimum a-value
  aMax = +0.99;        // maximum a-value
  for(int n = 0; n < N; n++)
  {
    double a = rsLinToLin(double(n), 0.0, N-1.0, aMin, aMax);
    ps[n] = phaseShapeRational(p[n], a);
    y0[n] = amp * sin(2*PI*ps[n] + 0.0*PI);
    y1[n] = amp * sin(2*PI*ps[n] + 0.5*PI);
    y2[n] = amp * sin(2*PI*ps[n] + 1.0*PI);
    y3[n] = amp * sin(2*PI*ps[n] + 1.5*PI);


    //ps[n] = phaseShapePower(p[n], pow(2.0, 2*a[n]));  // new, test
  }
  //rsPlotVector(ps);
  //rsPlotVectors(y0, y1);
  //writeToMonoWaveFile("PhaseShapeRatSine0.wav", &y0[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapeRatSine1.wav", &y1[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapeRatSine2.wav", &y2[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapeRatSine3.wav", &y3[0], N, (int) fs, 16);
  // Observations:
  // -The signal "PhaseShapeRatSine0.wav" morphs from a highpassed downward saw through sine to a
  //  highpassed upward saw. The morph looks like: highpass -> lowpass -> inverting highpass if we
  //  falsely imagine the signal as being generated by a filtered sawtooth wave. But what we want
  //  is: neutral -> lowpass -> inverted, i.e. we want: sawDown -> sine -> sawUp. The shaped phase
  //  goes from spikey (highpassed saw) through saw to sigmoid.
  // -PhaseShapeRatSine1 looks like DC+spikes -> sine -> DC+spikes
  // -Signals with a phase difference of pi are negatives of one another. So, it's actually 
  //  pointless to generate 4 signals. 0 and 1 are enough. 3 and 4 are just negatives of these.


  // Generate the signal with power-law phase-shaping map:
  aMin = -3.0;        // minimum a-value
  aMax = +3.0;        // maximum a-value
  for(int n = 0; n < N; n++)
  {
    double a = rsLinToLin(double(n), 0.0, N-1.0, aMin, aMax);
    ps[n] = phaseShapePower(p[n], pow(2.0, a));
    y0[n] = amp * sin(2*PI*ps[n] + 0.0*PI);
    y1[n] = amp * sin(2*PI*ps[n] + 0.5*PI);
    y2[n] = amp * sin(2*PI*ps[n] + 1.0*PI);
    y3[n] = amp * sin(2*PI*ps[n] + 1.5*PI);
  }
  rsPlotVector(ps);
  writeToMonoWaveFile("PhaseShapePowSine0.wav", &y0[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapePowSine1.wav", &y1[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapePowSine2.wav", &y2[0], N, (int) fs, 16);
  //writeToMonoWaveFile("PhaseShapePowSine3.wav", &y3[0], N, (int) fs, 16);
  // Observations:
  // -The signal using the power map looks nice at the beginning (similar to an unfiltered saw) but
  //  strange at the end. It's not like an inverted saw but like a squeezed one. May be useful but
  //  is not exactly what I wanted. Maybe what I want (inversion) can be achieved by some math. 
  //  Plug in the desired output signal and solve for the phase via asin?
  // -Try using:
  //  a >= 0: f(p,a) = powerLaw(p, 2^a), y = sin(f(p,a))
  //  a <  0: y = -sin(f(p,-a))  ->  asin(-y) = f(p,-a)
  //  where f(p,a) is the phase-shaping function with input phase p and parameter a. The idea is 
  //  that we first compute what the desired output signal should be (namely -sin(f(p,-a))) and 
  //  from that, compute the phase
  // -Or maybe instead of just inverting the power, we also need to reverse the direction and then
  //  flip up/down?
  // -Given f(p,a) we want to find a function g(p,a) such that sin(g(p,a)) = -sin(f(p,-a)). 
  //  Solving for g(p,a) gives g(p,a) = asin(-sin(f(p,-a))) = asin(sin(-f(p,-a))) = -f(p,-a)
  //  Note that in that equation, f(p,a) = phaseShapePower(p, 2^a)

  for(int n = 0; n < N; n++)
  {
    double a = rsLinToLin(double(n), 0.0, N-1.0, aMin, aMax);
    if(a <= 0)
      ps[n] = phaseShapePower(p[n], pow(2.0, a));
    else
    {
      //ps[n] = -phaseShapePower(p[n], pow(2.0, -a));     // (1)
      //ps[n] = 1 - phaseShapePower(p[n], pow(2.0, -a));  // (2)
      ps[n] = 1 - phaseShapePower(1-p[n], pow(2.0, -a));  // (3)
      //ps[n] =  -phaseShapePower(1-p[n], pow(2.0, -a));  // (4)

      // (5)
      //double tmp = phaseShapePower(p[n], pow(2.0, -a));
      //double y   = sin(2*PI*tmp);
      //ps[n]      = asin(-y) / (2*PI) + 0.5;
    }
    y0[n] = amp * sin(2*PI*ps[n] + 0.0*PI);
  }
  rsPlotVector(ps);
  writeToMonoWaveFile("PhaseShapePow2Sine0.wav", &y0[0], N, (int) fs, 16);
  // (1) and (2) seem to work but there's a discontinuity in the derivative. The phase reverses 
  // direction, I think. For (3), we don't have a phase discontinuity - but the saw goes up in any
  // case, i.e. at both ends. For (4), this is also true, but the phase is not in the right range.
  // So, (3) is best although it still doesn't do exactly what I want. (5) has the same behavior. 
  // Here, the problem could be due to multi-valuedness of asin. Maybe we need phase-unreflection. 
  // That would imply that the osc needs a state to remember the previous phase to achieve this 
  // disambiguation. I tend to think that we need to conditionally switch between (1) and (2) based
  // on what the phase (and maybe it's derivative) was at the previous sample. Maybe add 

  //
  // ToDo: 
  // -Experiment with different shaping functions. The goal is so morph between saws and regular 
  //  saw - not this highpassed looking one. Just like the ZDF-FM experiment. There, the phase 
  //  looks more like maximal slope is 2 rather than 1? -> figure that out!

  int dummy = 0;
}

void phaseShapingLinFrac()
{
  // We use the symmetrized linear fractional mapping function as phase-shaping function for a sine
  // wave with the goal of shaping it into a waveform that looks similar to a sawtooth wave. With
  // negative values for the slope parameter, we'll ge an upward saw and with positive ones, a
  // downward saw.

  using Real = double;
  using Vec  = std::vector<Real>;
  using LFI  = rsLinearFractionalInterpolator<Real>;

  // User parameters for the plots:
  int   N           = 257;         // Number of samples
  Real  shapePar    =  0.0;        // 0.0: symmetric (default) ...may be superfluous
  Real  minSlopePar = -4.0;
  Real  maxSlopePar = +4.0;
  int   numGraphs   =  9;


  // Helper function to produce the phase-shaping curve:
  auto createPhaseGraph = [](Real* x, Real* y, int N, 
                             Real slopeAt0, Real slopeAt1, Real shape = 0.0)
  {
    for(int n = 0; n < N; n++)
      y[n] = LFI::getNormalizedY(x[n], slopeAt0, slopeAt1, shape);
    // Try using LFI::symmetricMap. I think, it's simpler and does the same job.
  };

  // Helper function to create the phase-shaped sinusoid:
  auto createShapedSine = [&](Real* x, Real* y, int N,
    Real slopeAt0, Real slopeAt1, Real shape = 0.0)
  {
    createPhaseGraph(x, y, N, slopeAt0, slopeAt1, shape);
    for(int n = 0; n < N; n++)
      y[n] = sin(2*PI*y[n]);
  };


  // Maybe factor this plot out into its own function:
  Real slopeParInc = (maxSlopePar - minSlopePar) / (numGraphs - 1);
  Vec  x = rsLinearRangeVector(N, 0.0, 1.0);
  Vec  y(N);
  GNUPlotter plt;
  for(int i = 0; i < numGraphs; i++)
  {
    Real slopePar = minSlopePar + i * slopeParInc;
    Real slope    = pow(2.0, slopePar);

    // Uncomment one of the two:
    //createPhaseGraph(&x[0], &y[0], N, slope, slope, shapePar);
    createShapedSine(&x[0], &y[0], N, slope, slope, shapePar);

    plt.addDataArrays(N, &x[0], &y[0]);
  }
  plt.plot();

  int dummy = 0;




  // Observations:
  //
  // - We can indeed shape the sine into a waveform that looks similar to a sawUp or sawDown wave.
  //
  // - shapePar does not seem to make a visible difference (verify!)
  //
  //
  // ToDo:
  //
  // - Maybe allow moving the split point around like in phaseShapingCurvesRational().
  //
  // - To use this for an oscillator, maybe negate the slope parameter such that positive values
  //   turn a sine into an upward saw. Maybe write a class SinSawOscLinFrac for this. Maybe use it
  //   in the bassdrum synth. Or make it a more general PhaseShapeOsc with various options for the
  //   shaping function.
  //
  // - Make a morphable SawDown/Sin/SawUp oscillator based on the TriSaw oscillator by applying a
  //   sinusoidal waveshaper to its output
  //
  // - Implement a class rsPhaseShaper. There's a stub for this in rosic_MiscUnfinished.h
  //
  //
  // See also: 
  //
  // - linearFractionalInterpolation()
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

  a = 0.8;  // test

  double tol = 1.e-13;
  int maxIts = 100;     // maximum number of Newton steps, if 0, we get unit-delay feedback PM

  double w = 2*PI/cycle;
  std::vector<double> y(N), z(N), q(N);

  double p;

  // Define objective function for root-finder:
  auto func = [&](double y) { return y - sin(p + a*y); };


  for(int n = 1; n < N; n++)           // starts from 1 to enable UDF, y[0]=z[0]=q[0]=0 anyway
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
  // differentiate -> unwrap -> integrate (numerically)? But wrap at what value? Or use the more 
  // advanced phase-unwrapping algorithms in rsSingleSineModeler. It has these unreflectPhase...
  // functions. The phase shouldn't go backward.

  using SSM = rsSingleSineModeler<double>;
  SSM::unreflectPhaseFromSig(&z[0], &q[0], N);


  //rsPlotVectors(y);  // only UDF signal
  rsPlotVectors(q);    // only the phase
  rsPlotVectors(y, z, z-y, q);
  int dummy = 0;

  // Observations:
  // -The phase looks qualitatively similar to what is produced in the phaseShapingCurvesRational
  //  experiment. That means that with such a phase-shaping algo, we could produce a similar signal
  //  with an explicit algorithm that doesn't need Newton iteration. phaseShapingSkew seems to try 
  //  that


  // Resources:
  // https://ristoid.net/modular/fm_variants.html
  // https://ccrma.stanford.edu/software/snd/snd/fm.html
}