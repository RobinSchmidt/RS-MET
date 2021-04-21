#include "FilterExperiments.h"

//using namespace RAPT; // remove when possible

/*

ToDo:

make a test for intermodulation distortion:
feed 2 sinewaves  with frequencies f1 and f2, for
example f1 = 1000 Hz, f2 = 1082 Hz -> the frequencies should be chosen such that multiples of the
difference frequency and multiples of the sum-frequency are easily distinguished (if several
saturators are applied in series, we'll see harmonics of sum- and difference frequncies as well9
...maybe find an even better frequency ratio

make a test for time-varying behavior:
feed a square-wave (naive, with aliasing) at some frequency around 100 Hz (pick one where aliasing
is very apparent) into a filter, switch the cutoff of the filter back and forth between 2000 and
500 Hz, at a rate of, say, 5 Hz -> look at time-domain output


Ladder filter:
introduce a nonlinearity of the form f(x) = ((1+a)*x)/(1+a*x^2)
a determines the amount of distortion, the slope at the origin will be 1+a (verify this)
...try to find an analytic solution for the instantaneous feedback problem (this will depend on
where the nonlinerity is placed), if this is not possible, at least, with rational functions, it
should be possible to find reasonable expressions for the derivative to be used in a newton
iteration, also nice: f(x) = (x+a*x^3)/(1+a*x^4) - very linear in the range -1...+1, for low a
maybe in the feedback-path, we need another nonlinearity, we'll see
-introduce a pole-spread: scale poles such that the frequencies are:
f1 = fc/b, f2 = fc/sqrt(b), f3 = fc*sqrt(b), f4 = fc*b, where
b = sqrt(2^poleSpread) -> poleSpread is in octaves

Biquad:
implement a biquad based on a rotating complex phasor - this should produce the least amount of
  artifacts when the cutoff is switched. the rotating phasor implements the two poles, 2 zeros may
  be placed afterwards or before - see what gives best results
-maybe apply distortion by using the rational functions above in the complex domain

*/


// Example code for a TPT/ZDF-SVF (obviously, this could be optimized for production code),
// written by Robin Schmidt, licensing/copyright: none (public domain)
// Parameters:
// fc: cutoff frequency
// fs: sample-rate
// R:  damping coefficient (todo: maybe replace this later by a resonance gain, maybe pass an array)
// N:  number of samples
// x:  input signal
// yL: lowpass output
// yB: bandpass output
// yH: highpass output
void svf(double fc, double fs, double R, int N, double x[], double yL[], double yB[], double yH[])
{
  double wd = 2*PI*fc;           // target radian frequency
  double T  = 1/fs;              // sampling period
  double wa = (2/T)*tan(wd*T/2); // prewarped radian frequency for analog filter (Eq. 3.7)
  double g  = wa*T/2;            // embedded integrator gain (Fig 3.11), wc == wa

  // states of the 2 integrators, static so we can call the function block-wise while maintaining 
  // the states from block to block:
  static double s1 = 0.0;
  static double s2 = 0.0;

  // loop over the samples:
  for(int n = 0; n < N; n++)
  {
    // compute highpass output via Eq. 5.1:
    yH[n] = (x[n] - 2*R*s1 - g*s1 - s2) / (1 + 2*R*g + g*g);

    // compute bandpass output by applying 1st integrator to highpass output:
    yB[n] = g*yH[n] + s1;
    s1    = g*yH[n] + yB[n]; // state update in 1st integrator

    // compute lowpass output by applying 2nd integrator to bandpass output:
    yL[n] = g*yB[n] + s2;
    s2    = g*yB[n] + yL[n]; // state update in 2nd integrator

    // Remark: we have used two TDF2 integrators (Fig. 3.11) where one of them would be in code:
    // y = g*x + s; // output computation
    // s = g*x + y; // state update

    // as a cheap trick to introduce nonlinear behavior, we apply a nonlinearity to state of 2nd 
    // integrator:
    //s1 = tanh(s1); // maybe don't do this
    //s2 = tanh(s2);
    //int dummy = 0;
  }
}


void rsAnalogBandpassCoeffsConstantPeak(double *B0, double *B1, double *B2, double *A0, double *A1,
  double *A2, double Q, double wc = 1.0)
{
  // give the function another parameter to switch between constant peak and constant skirt
  // move to RSLib

  *B0 = 0;
  *B1 = 1/(wc*Q);
  *B2 = 0;
  *A0 = 1;
  *A1 = *B1;
  *A2 = 1/(wc*wc);
}

/*
void bandwidthScaling()
{
  // We plot magnitude responses of analog bandpass filters which are to be applied multiple times
  // (1, 2, 4, 8, 16 times) to the same signal. If the filter is applied multiple times, the
  // bandwidth of a single-pass filter must be scaled (broadened) such that the cutoff frequency of
  // the multipas filter matches that of a nominal single pass filter.

  // user parameters:
  static const int N = 1000;   // number of frequencies
  double wMin  =  0.001;       // minimum radian frequency for plot
  double wMax  = 1000.0;       // maximum radian frequency for plot
  double wc    = 1.0;          // radian center frequency
  double Q1    = 0.1;          // quality factor of 1-pass filter
  double gSq   = 0.5;          // target squared gain where curves should cross
  int    order = 3;            // order of the single-pass filter

  // internal variables:
  double B0, B1, B2, A0, A1, A2; // analog biquad coeffs
  double w[N], p[5][N];
  rsFillWithRangeExponential(w, N, wMin, wMax);

  // compute magnitude responses, convert to dB and plot:
  double Q;  // scaled Q
  int np;    // np number of passes
  for(int k = 0; k < 5; k++)
  {
    np = rsPowInt(2, k);
    Q = Q1 / rsBandwidthConverter::multipassScalerButterworth(np, order, gSq);
    rsAnalogBandpassCoeffsConstantPeak(&B0, &B1, &B2, &A0, &A1, &A2, Q, wc);
    for(int n = 0; n < N; n++)
    {
      p[k][n] = analogBiquadMagnitudeSquaredAt(B0, B1, B2, A0, A1, A2, w[n]);

      p[k][n] = pow(p[k][n], 0.5*np*order);
        // nope -this is wrong - this is not a Butterworth response with given order applied
        // np times

      p[k][n] = rsMax(rsAmp2dB(p[k][n]), -200.0);
    }
  }
  plotDataLogX(N, w, p[0], p[1], p[2], p[3], p[4]);
}
*/

void bandwidthScaling()
{
  // user parameters:
  static const int N = 1000;  // number of frequencies
  double wMax = 2.0;          // maximum radian frequency for plot
  double wc   = 1.0;          // radian center frequency
  double Q1   = 2.0;          // quality factor of 1-pass filter
  double g    = 0.5;          // target gain where curves should cross

  // internal variables:
  double B0, B1, B2, A0, A1, A2; // analog biquad coeffs
  double w[N], p[5][N];
  RAPT::rsArrayTools::fillWithRangeLinear(w, N, 0.0, wMax);

  // compute and plot magnitude-squared responses:
  double Q;
  for(int k = 0; k < 5; k++)
  {
    Q = Q1 / RAPT::rsBandwidthConverter::multipassScalerButterworth(k+1, 1, g);
    rsAnalogBandpassCoeffsConstantPeak(&B0, &B1, &B2, &A0, &A1, &A2, Q, wc);
    for(int n = 0; n < N; n++)
    {
      p[k][n] = analogBiquadMagnitudeSquaredAt(B0, B1, B2, A0, A1, A2, w[n]);
      p[k][n] = pow(p[k][n], k+1);
    }
  }
  plotData(N, w, p[0], p[1], p[2], p[3], p[4]);
}

void butterworthEnergy()
{
  int N = 5;  // Butterworth order
  int M = 8;  // number of passes
  double E = RAPT::rsPrototypeDesigner<double>::butterworthEnergy(N, M);

  // todo: maybe make plots of both normalizing functions (energy and cutoff based

  int dummy = 0;
}


// these tow are obsolete - the code has been moved into rsOnePoleFilter::prepareForBackwardPass
template<class T>
T rsOnePoleInitialStateForBackwardPass(T a, T b, T y)
{
  return (a*y*b) / (T(1)-b*b);
  /*
  T k = T(1) / b; 
  T p = a*y  / (b*b-T(1));
  T c = -p*k; 
  return c - p*(b-k);
  */
}
// -rename a,b to b0,a1, maybe change sign convention for feedback coeff
// -dividing by 1-b^2 may be numerically inaccurate, if |b| is close to 1, which is typical for 
//  filters - can the formula be rewritten in a numerically more accurate way?

template<class T>
void rsOnePoleInitialStateForBackwardPass(T a, T b, T* y1, T d, T* x1)
{
  T q = b * *y1 + d * *x1;
  T p = q*(a+b*d) / (b*(b*b-T(1)));
  T k = T(1) / b;
  T c = -p*k;
  *x1 = q;            // == t[1]
  *y1 = c - p*(b-k);  // == s[1]
}
// -try to simplify
// -move to rsOnePoleFilter
// -have a member prepareForBackwardPass that does that to the internal state variables
// -rename a,b,d to b0,a1,b1
//   k = 1/b; 
//   p = a * (b*y + d*x) / (b*(b^2-1));
//   c = -p*k;
//   yNew = c - p*(b-k); // == s(1)
//   xNew = t[1] = b*y + d*x

void biDirectionalStateInit1()
{
  // Implementation of and experiments with closed form expressions to set the internal state of a 
  // 1st order bidirectional filter directly before the backward run. We compare the results from 
  // the direct setting of the state and the ring-out/warm-up strategy.
  // see BiDirectionalInitialStates.txt for more details

  // filter coeffs:
  double a = 4.0;
  double b = 0.5;
  double d = -0.25;

  // input signal:
  static const int N = 7;
  double x[N] = { 1,2,3,4,7,4,3 };

  // tail buffers:
  static const int Nt = 200;  // tail buffer length
  double t[Nt], s[Nt];

  // output after forward pass and after both passes:
  double yf[N], y[N];

  // create and set up filter object:
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setCoefficients(a, d, b);  // we use different conventions here

  // compute forward pass output:
  int n;
  for(n = 0; n < N; n++)
    yf[n] = flt.getSample(x[n]);

  // Compute filter state analytically - the computed value should match t[0] (not yet computed) up 
  // to roundoff (if Nt is too short, maybe worse). t[0] can be considered as having computed the 
  // desired inital state numerically:
  double xOld =  x[N-1];                 // last input into filter, equals its x[n-1] state
  double yOld = yf[N-1];                 // last output of filter, equals its y[n-1] state
  double q = b*yOld + d*xOld;
  double p = q*(a+b*d) / (b*(b*b-1));
  double k = 1 / b;
  double c = -p*k;                       // constant determined by boundary condition s[inf] == 0
  double xNew = q;                       // == t[1], first sample of forward tail
  double yNew = c - p*(b-k);             // == s[1], first sample of bidirectional tail

  // test, if the function computes the same:
  double xTmp = xOld, yTmp = yOld; // should become equal to xNew, yNew:
  rsOnePoleInitialStateForBackwardPass(a, b, &yTmp, d, &xTmp); // ok, looks good


  // ring-out/warm-up, using tail buffers:
  for(n = 0;    n <  Nt; n++) t[n] = flt.getSample(0.0);   // fill forward tail buffer, ring out
  for(n = Nt-1; n >= 0;  n--) s[n] = flt.getSample(t[n]);  // fill backward tail buffer, warm up

  // compare numerically and analytically computed states:
  double err = yNew - s[0]; // should be zero (up to roundoff)
  rsAssert(rsAbs(err) < 1.e-12); 
  err = xNew - t[0];
  rsAssert(rsAbs(err) < 1.e-12); 


  // compute backward pass output:
  for(n = N-1; n >= 0; n--)
    y[n] = flt.getSample(yf[n]);

  // Compute the tail generated by the closed form expression:
  double sa[Nt];    // analytic tail
  for(n = 0; n < Nt; n++)
  {
    //sa[n] = c*pow(k, n-1) - (a * yOld * (pow(b, n) - pow(k, n))) / (b*b-1); // only correct if d==0
    sa[n] = c*pow(k, n-1) - (a * (b*yOld+d*xOld) * (pow(b, n) - pow(k, n))) / (b*(b*b-1));
  }
  // Looks good, if b is negative power of 2 and d==0, otherwise numerical error messes up the
  // analytic solution. Result ta is equal to t but one sample shifted - when plugging n=1 into the 
  // formula, we can compute our desired state variable yNew.


  // Compare results from first doing a forward pass, then a backward pass and the other way 
  // around: backward first, then forward - the results should be the same:
  double xL = -7, xR = -3;                                     // left and right boundary conditions
  double yfb[N]; flt.applyForwardBackward(x, yfb, N, xL, xR);  // yfb: forward/backward output
  double ybf[N]; flt.applyBackwardForward(x, ybf, N, xL, xR);  // ybf: backward/forward output
  double diff[N]; rsArrayTools::subtract(yfb, ybf, diff, N);   // diff = yfb-ybf, should be zero


  // plot results:
  //rsPlotArrays(Nt, s, sa);     // numerically and analytically computed tail
  //rsPlotArrays(Nt, t, s);   // tail buffers
  //rsPlotArrays(N, x, yf, y); // input and outputs
  rsPlotArrays(N, x, yfb, ybf, diff);  // compare forward-first vs backward-first


  // solving recursions:
  // https://ask.sagemath.org/question/35626/how-to-solve-non_homogenous-recurrence-in-sage-using-rsolve-functions/
  // https://trac.sagemath.org/ticket/1291
  // https://groups.google.com/forum/#!topic/sage-support/pYvjN7da9LY
  // Computational Mathematics with SageMath, page 229
}

template<class T>
void rsOnePoleInitialStateForBackwardPass(T a, T b, T d, T r, T* x1, T* y1)
{
  // compute some intermediate values:
  T b2 = b*b;                               // b^2
  T b3 = b*b2;                              // b^3
  T cx =    a*b + (b2-b)*d  - a;            // coeff for r in update of y1
  T cr = a*(a*b + (b +1)*d) + d*d;          // coeff for x1 in update of y1

  // compute new state variables:
  *x1  = a*r + b * *y1 + d * *x1;           // update x1
  *y1  = (cr*r - cx * *x1) / (b3-b2-b+1);   // update y1 using updated x1
}

void biDirectionalStateInit2()
{
  // We try the generalized formula where we don't assume the input to go to zero but to some 
  // arbitrary constant value r.

  // filter coeffs:
  double a =  4.0;
  double b =  0.5;
  double d = -0.25;

  // final state:
  double x1 = 0.2;
  double y1 = 0.3;

  // target constant:
  double r = -0.2;

  static const int N = 40;  // tail length
  int n;

  // compute forward tail numerically, i.e. with the filter:
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setCoefficients(a, d, b);
  flt.setInternalState(x1, y1);
  double tn[N];
  tn[0] = y1;
  for(n = 1; n < N; n++)
    tn[n] = flt.getSample(r);

  // compute forward tail analytically, i.e. via the formula t[n] = v + (q-v)*b^(n-1):
  double ta[N];
  double q = a*r + b*y1 + d*x1;
  double u = (a+d)*r;
  double v = u/(1-b);          // == t[inf]
  ta[0] = y1;
  ta[1] = q;
  for(n = 2; n < N; n++)
    ta[n] = v + (q-v) * pow(b,n-1);

  // compute forward/backward tail numerically:
  flt.setStateForConstInput(v);  // x1 == t[inf] == v, y1 == s[inf] == (a+d)*v / (1-b)
  double sn[N];
  double xNew, yNew;             // desired new states - will be recorded at n==1
  for(n = N-1; n >= 0; n--)
  {
    sn[n] = flt.getSample(tn[n]);
    if(n == 1) {
      xNew = flt.getStateX();
      yNew = flt.getStateY(); }
  }

  // compute forward/backward tail analytically:
  double sa[N];
  double k = 1/b;
  double k2 = k*k;   // k^2
  double k3 = k2*k;  // k^3
  double c  = (k2*(a*(k2*(q-v)+k*v+v)+d*(k*q+v)))/(k2-1); // try to simplify
  for(n = 2; n < N; n++)
  {
    double K = pow(k, n-1);

    // try to simplify this:
    sa[n] = c*K - ( ((K*k - 1)/K) * (a* (K*k3* (q-v) 
                    + v * K*k2 + v * K*k + k2 * (q-v)) + d * (K*k2 * (q-v) 
                    + v * K*k2 + v * K*k + k  * (q-v))))/(k2 - 1);
  }
  sa[1] = a*ta[1] + b*sa[2] + d*ta[2];
  sa[0] = a*ta[0] + b*sa[1] + d*ta[1];

  // test the state-setter function - after this call, our local variables xx1,yy1 should match the
  // corresponding recorded state variables xNew,yNew at n == 1:
  double xx1 = x1, yy1 = y1;
  rsOnePoleInitialStateForBackwardPass(a, b, d, r, &xx1, &yy1);
  rsAssert(xx1 == xNew && yy1 == yNew);  // ...yes! they do match indeed! :-)

  // test - after these two calls, the filter's state should be xNew, yNew:
  flt.setInternalState(x1, y1);
  flt.prepareForBackwardPass(r);
  rsAssert(flt.getStateX() == xNew && flt.getStateY() == yNew);

  // plot forward tails and forward/backward tails, numerically and analytically evaluated:
  rsPlotArrays(N, tn, ta, sn, sa);  // OK - that looks good. 

  // todo: try if it makes a difference if we apply the filter forward first or backward first - it
  // doesn't make any with the simpler formulas where r=0 - but verify that it still makes no 
  // difference for arbitrary choices of xL,xR
}

void biDirectionalStateInit()
{
  biDirectionalStateInit1();    // assumes the input to be zero outside n = 0,...,N-1
  biDirectionalStateInit2();    // assumes the input to go to some user selectable constants for
                                // n < 0 and n >= N (they may be different for each side)
}

void biquadTail()
{
  // We check the explicit formula for a biquad tail against the tail computed by the actual 
  // filter.

  // coeffs:
  double b0 =  8;
  double b1 =  4;
  double b2 =  2;
  double a1 = -0.5;
  double a2 = -0.25;

  // state variables:
  double x1 = -4;
  double x2 =  2;
  double y1 =  8;
  double y2 = -4;

  // number of samples:
  int N = 20;  

  // Notation change to have single letter variables in the formulas:
  //   y[n] = u*x[n] + v*x[n-1] + w*x[n-2] + a*y[n-1] + b*y[n-2]
  double u = b0, v = b1, w = b2;
  double a = -a1, b = -a2;

  // compute intermediate variables (needs to be done before computing the tail using the filter)
  double t0 = y1;                         // 0th tail sample is last actual output sample y[N-1]
  double t1 = v*x1 + w*x2 + a*y1 + b*y2;  // computed directly
  double t2 = w*x1 + a*t1 + b*t0;         // computed directly
  //double t3 = a*t2 + b*t1;                // computed recursively
  //double t4 = a*t3 + b*t2;                // computed recursively
  double g = t1, h = t2;
  double S = sqrt(a*a + 4*b);             // may have to be complex when a*a + 4*b < 0
  double T = a - S, U = a + S;            // ...then these all have to be complex as well
  double V = a*a*g + 2*b*g - a*h;
  double W = a*g*S - h*S;
  double A = W+V;
  double B = W-V;
  // g = -9, h = -10.5, a = 0.5, b = 0.25

  // compute tail using the explicit formula:
  //   t(n) = -(A*T^n + B*U^n) / (b*S*2^(n+1))
  // which has been found by wolfram alpha using
  //   RSolve[{t[n] == a*t[n-1] + b*t[n-2], t[1] == g, t[2] == h}, t[n], n]
  // and doing a lot of simplification by hand
  using Vec = std::vector<double>;
  Vec tFrm(N);
  for(int n = 0; n < N; n++)
    tFrm[n] = -(A*pow(T,n) + B*pow(U,n)) / (b*S*pow(2,n+1));

  // compute tail using the filter:
  Vec tFlt(N);
  for(int n = 0; n < N; n++) {
    double x0 = 0;
    tFlt[n] = u*x0 + v*x1 + w*x2 + a*y1 + b*y2;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = tFlt[n]; 
  }

  // plot both results:
  rsPopFront(tFrm); // because they agree only up to a shift by one sample
  rsPlotVectors(tFlt, tFrm); 
  int dummy = 0;
}

void biquadModulation()
{
  // Compares the responses of different biquad implementation structures (DF1, DF2, SVF, etc.) to
  // modulation of the parameters...
  
  // not yet finsihed

  int    N  = 10000;         // number of samples
  double fs = 44100;         // sample rate
  double f  = 50;            // input signal frequency
  double f1 = 5000;          // first filter frequency
  double f2 = 500;           // second filter frequency
  double q  = 5.0;           // Q value
  //double g  = 2;             // gain factor (for shelf or bell filters)


  double b0, b1, b2, a1, a2; // biquad coeffs
  typedef rosic::BiquadDesigner DSN;
  rosic::BiquadMonoDF1 bqd;
  rsStateVectorFilter<double, double> stVecFlt;

  // has to be put into a loop:
  //DSN::calculateCookbookPeakFilterCoeffsViaQ(b0, b1, b2, a1, a2, 1/fs, f, q, g);
  DSN::calculateCookbookLowpassCoeffs(b0, b1, b2, a1, a2, 1/fs, f, q);
  // uses + convention for a-coeffs
  //bqd.setCoefficients(b0, b1, b2, a1, a2);
  //stVecFlt.setupFromBiquad(b0, b1, b2, -a1, -a2);



  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, &yBqd[0]);
  //plt.addDataArrays(N, &yStVec[0]);
  plt.plot();
}

void stateVariableFilter()
{
  double fs = 44100;  // samplerate in Hz
  double fc = 5000;   // cutoff/center frequency in Hz
  double G  = 4.0;    // cutoff gain or peak/shelf gain
  double B  = 0.25;   // bandwidth (for peak, const-peak-bandpass, bandreject)

  static const int N = 2048;  // number of sample to generate

  rsStateVariableFilterDD svf;
  svf.setSampleRate(fs);
  svf.setFrequency(fc);
  svf.setGain(G);
  svf.setBandwidth(B);  //later

  //svf.setBandwidth(0.1);  // test

  // time axis and impulse responses:
  double t[N], yL[N], yB[N], yH[N], yPK[N], yLS[N], yHS[N], yAP[N], yBR[N], yBPCP[N];

  // frequency axis and magnitude responses:
  double f[N/2], mL[N/2], mB[N/2], mH[N/2], mPK[N/2], mLS[N/2], mHS[N/2], mAP[N/2], pAP[N/2];
  double mBR[N/2], mBPCP[N/2];

  // create time-axis and frequency-axis scaled in Hz (put into function):
  createTimeAxis(N, t, fs);
  for(int k = 0; k < N/2; k++)
    f[k] = k*fs/N;

  // obtain impulse-responses for various modes:
  svf.setMode(rsStateVariableFilterDD::LOWPASS);
  svf.reset();
  getImpulseResponse(svf, yL, N);

  svf.setMode(rsStateVariableFilterDD::HIGHPASS);
  svf.reset();
  getImpulseResponse(svf, yH, N);

  svf.setMode(rsStateVariableFilterDD::BANDPASS_SKIRT);
  svf.reset();
  getImpulseResponse(svf, yB, N);

  svf.setMode(rsStateVariableFilterDD::BANDPASS_PEAK);
  svf.reset();
  getImpulseResponse(svf, yBPCP, N);

  svf.setMode(rsStateVariableFilterDD::BANDREJECT);
  svf.reset();
  getImpulseResponse(svf, yBR, N);

  svf.setMode(rsStateVariableFilterDD::BELL);
  svf.reset();
  getImpulseResponse(svf, yPK, N);

  svf.setMode(rsStateVariableFilterDD::LOWSHELF);
  svf.reset();
  getImpulseResponse(svf, yLS, N);

  svf.setMode(rsStateVariableFilterDD::HIGHSHELF);
  svf.reset();
  getImpulseResponse(svf, yHS, N);

  svf.setMode(rsStateVariableFilterDD::ALLPASS);
  svf.reset();
  getImpulseResponse(svf, yAP, N);

  // get magnitude responses from impulse responses:
  rsMagnitudeAndPhase(yL, N, mL);
  rsMagnitudeAndPhase(yB, N, mB);
  rsMagnitudeAndPhase(yH, N, mH);
  rsMagnitudeAndPhase(yPK, N, mPK);
  rsMagnitudeAndPhase(yLS, N, mLS);
  rsMagnitudeAndPhase(yHS, N, mHS);
  rsMagnitudeAndPhase(yAP, N, mAP, pAP); RAPT::rsArrayTools::scale(pAP, N/2, 180.0/PI);
  rsMagnitudeAndPhase(yBR, N, mBR);
  rsMagnitudeAndPhase(yBPCP, N, mBPCP);

  // magnitude response plots:
  //plotData(N/2, f, mL, mB, mH);    // LP, BP, HP - the bandwidth doesn't have any effect - why?
  plotData(N/2, f, mBPCP, mBR);    // BR and BP (const peak gain)
  //plotData(N/2, f, mPK, mLS, mHS); // bell, low- and highshelf
  //plotData(N/2, f, mAP, pAP);      // allpass, plot also phase
}

void stateVariableFilterMorph()
{
  double fs = 44100;  // samplerate in Hz
  double fc = 10000;  // cutoff/center frequency in Hz
  double G  = 2.0;    // cutoff gain or peak/shelf gain
  double B  = 2.0;    // bandwidth (for peak, bandpass, bandreject)

  static const int N = 2048;  // number of sample to generate

  rsStateVariableFilterDD svf;
  svf.setSampleRate(fs);
  svf.setFrequency(fc);
  svf.setGain(G);
  svf.setBandwidth(B);  //later

  // impules responses:
  double y0[N], y1[N], y2[N], y3[N], y4[N]/*, y5[5]*/;

  // frequency axis and magnitude responses:
  double f[N/2], m0[N/2], m1[N/2], m2[N/2], m3[N/2], m4[N/2]/*, m5[N/2]*/;

  // frequency-axis scaled in Hz (put into function):
  for(int k = 0; k < N/2; k++)
    f[k] = k*fs/N;

  // obtain impulse-responses for various modes:
  svf.setMode(rsStateVariableFilterDD::MORPH_LP_BP_HP);
  svf.setMorph(0.0);
  svf.reset();
  getImpulseResponse(svf, y0, N);

  svf.setMorph(0.125);
  svf.reset();
  getImpulseResponse(svf, y1, N);

  svf.setMorph(0.25);
  svf.reset();
  getImpulseResponse(svf, y2, N);

  svf.setMorph(0.375);
  svf.reset();
  getImpulseResponse(svf, y3, N);

  svf.setMorph(0.5);
  svf.reset();
  getImpulseResponse(svf, y4, N);

  /*
  svf.setMorph(0.5);
  svf.reset();
  getImpulseResponse(svf, y5, N);
  */


  // get magnitude responses from impulse responses:
  rsMagnitudeAndPhase(y0, N, m0);
  rsMagnitudeAndPhase(y1, N, m1);
  rsMagnitudeAndPhase(y2, N, m2);
  rsMagnitudeAndPhase(y3, N, m3);
  rsMagnitudeAndPhase(y4, N, m4);
  //rsMagnitudeAndPhase(y5, N, m5);



  // magnitude response plot:
  plotData(N/2, f, m0, m1, m2, m3, m4);
}

void stateVectorFilter()
{
  // Idea: create a filter similar to the complex-phasor filter, but more general, using a 2D
  // (x,y) vector as state and update equations:
  //
  // t = x;                       // temporary
  // x = xx*x + xy*y + in;        // update x
  // y = yx*t + yy*y + in;        // update y
  // return cx*x + cy*y * ci*in;  // output
  //
  // It should be able to have a pair of complex conjugate poles, if the (xx,xy;yx,yy) matrix is
  // a spiraling matrix (xx = yy = r*cos(w), -xy = yx = r*sin(w)) or a parallel connection of two
  // real poles (xy = yx = 0), the output coeffs cx,cy,ci determine the zeros.
  //
  // To find the coefficients, take a biquad reference filter, compute the first 5 samples of its
  // impulse response, find formulas of the state-vector-filter's first 5 impulse response samples
  // in terms of its coeffs and solve the resulting system of equations for the coeffs (maybe use
  // the additional constraints on the matrix coeffs above, if necessary - maybe we need to choose
  // which constraint to use based on the discriminant of the a1,a2 biquad coeffs). To make a 
  // general modulatable filter, boil it down to a parallel connection of biquads and implement 
  // each biquad as state-vector filter - make comparison experiments of modulation properties of
  // DF1, DF2, TDF1, TDF2, state-variable and state-vector implementations. Maybe if the 
  // equation system based on the impulse response is hard or impossible to solve analytically, try
  // something based on the transfer function - maybe it makes sense to derive formulas for the
  // coeffs in terms of poles and zeros - how? ..or maybe first derive formulas for poles and zeros
  // in terms of the coeffs? maybe compute matrix coefficients from the two poles and weights
  // cx,cy,ci from the 1st 3 samples of the impulse repsonse when xx,xy,yx,yy are already known
  // ...this should be tractable
  //
  // 0: x = 0, y = 0, in = 1
  //    x = xx*0 + xy*0 + 1 = 1
  //    y = yx*0 + yy*0 + 1 = 1
  //    h[0] = cx + cy + ci
  // 1: x = 1, y = 1, in = 0
  //    x = xx*1 + xy*1 + 0 = xx + xy
  //    y = yx*1 + yy*1 + 0 = yx + yy
  //    h[1] = cx*(xx + xy) + cy*(yx + yy)
  // 2: x = xx+xy, y = yx+yy, in = 0
  //    x = xx*(xx+xy) + xy*(yx+yy)
  //    y = yx*(xx+xy) + yy*(yx+yy)
  //    h[2] = cx*(xx*(xx+xy) + xy*(yx+yy)) + cy*(yx*(xx+xy) + yy*(yx+yy))
  // uh oh - expression complexitiy blows up exponentially and we get lots on products, i.e. a 
  // nonlinear system...but it's a 3x3 linear system, if we compute xx,xy,yx,yy from the poles
  // beforehand and treat them as knowns - we need a branch for computing these for the two cases
  // of two real poles and a pair of complex conjugate poles
  // h[3] = 
  // h[4] = 

  // transfer functions (verify):
  // X(z) = 1 / (1 - xx * z^-1 * X(z) - xy * z^-1 * Y(z))
  // Y(z) = 1 / (1 - yx * z^-1 * X(z) - yy * z^-1 * Y(z))
  // H(z) = ci + cx*X(z) + cy*Y(z)
  // maybe solve 2nd equation for Y(z) and plug into 1st - try to equate to biquad transfer
  // function and compare coeffs
  // H(z) =
  // try to manipulate with sage:
  // var("H z X Y ci cx cy xx xy yx yy")
  // H = ci + cx/(1-(xx*X + xy*Y)/z) + cy/(1-(yx*X+yy*Y)/z)
  // simplify(H) # doesn't help
  // i want to write down H(z) in terms of xx,xy,yx,yy,ci,cx,cy and z without using X(z),Y(z)
  // H.simplify_full()
  // maybe eliminate(X)? or maybe i need a matrix z-transform

  int    N  = 100;           // number of samples
  double fs = 44100;         // sample rate
  double f  = 1000;          // filter frequency
  double q  = 0.5;           // Q value
  double g  = 2;             // gain factor (for shelf or bell filters)

  typedef rosic::BiquadDesigner DSN;
  double b0, b1, b2, a1, a2; // biquad coeffs
  //DSN::calculateCookbookPeakFilterCoeffsViaQ(b0, b1, b2, a1, a2, 1/fs, f, q, g);
  DSN::calculateCookbookLowpassCoeffs(b0, b1, b2, a1, a2, 1/fs, f, q);
   // uses + convention for a-coeffs

  //// test:
  //b0 = 1;
  //b1 = 0;
  //b2 = 0;
  //a1 = 0.5;
  //a2 = 0.2;

  rosic::BiquadMonoDF1 bqd;
  rsStateVectorFilter<double, double> stVecFlt;

  bqd.setCoefficients(b0, b1, b2, a1, a2);
  stVecFlt.setupFromBiquad(b0, b1, b2, -a1, -a2);

  std::vector<double> yBqd   = impulseResponse(bqd,      N, 1.0);
  std::vector<double> yStVec = impulseResponse(stVecFlt, N, 1.0);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &yBqd[0]);
  plt.addDataArrays(N, &yStVec[0]);
  plt.plot();

  // Observations:
  // -cookbook lowpass with q = 0.5 leads to a singularity (a matrix becomes non-invertible)
  //  -in this case, we fudge with the poles (todo: make more tests for this)
}

void transistorLadder()
{

  double fs     = 44100;  // samplerate in Hz
  double length = 2.1;    // length of signal in seconds


  int N = (int)(length*fs); // number of sample to generate


  double *t = new double[N];
  double *x = new double[N];
  createTimeAxis(N, t, fs);

  RAPT::rsArrayTools::fillWithRandomValues(x, N, -0.01, +0.01, 1);
  //fillWithZeros(x, N);
  //x[0] = 1.0;


  rsLadderMystran<double, double> ldr;
  ldr.setMode(rsLadderMystran<double, double>::LP_24);
  ldr.setCutoff(1000.0);
  ldr.setResonance(1.0);
  for(int n = 0; n < N; n++)
  {
    // apply envelope (exponential sweepdown to create a 'zap' sound)
    //double f = 55.0 + 20000.0 * exp(-n / (0.01*fs));
    //tl.setCutoff(f);

    x[n] = ldr.getSample(x[n]);
  }


  // todo: compare against the non-ZDF ladder, try multimode    
  plotData(N, t, x);

  RAPT::rsArrayTools::normalize(x, N, 1.0);
  writeToMonoWaveFile("d:\\TmpData\\MystranLadderZap.wav", x, N, (int)fs, 16);

  delete[] t;
  delete[] x;
}

void phonoFilterPrototypePlot()
{
  // create an exponentially scaled frequency axis:
  static const int N = 500; // number of frequencies for the plot
  double fL = 1;            // low frequency
  double fH = 1.e+6;        // high frequency
  double f[N], dB[N];
  for(int n = 0; n < N; n++)
  {
    f[n]  = rsLinToExp((double)n, 0.0, (double)N-1, fL, fH);
    dB[n] = rsAmp2dB(rsPhonoFilterDD::prototypeMagnitudeAt(f[n]));
  }
  plotDataLogX(N, f, dB);
}
void magnitudeMatchedOnePoleFilter()
{
  double fs = 44100;  // sample-rate
  int i;

  // frequency/magnitude pairs:
  double f[3], m[3];
  f[0] = 1000.0; m[0] = 2.0;
  f[1] = 2000.0; m[1] = 1.5;
  f[2] = 4000.0; m[2] = 1.0;

  // normalized radian frequencies:
  double w[3];
  for(i = 0; i < 3; i++)
    w[i] = 2*PI*f[i]/fs;

  // compute filter-coeffs:
  double b0, b1, a1;
  magnitudeMatchedOnePoleCoeffs(b0, b1, a1, w, m);

  // compute and plot the magnitude response
  static const int N = 500;
  double fp[N], mp[N];
  for(i = 0; i < N; i++)
  {
    fp[i] = i * fs/(2*N);
    mp[i] = onePoleMagnitudeAt(b0, b1, a1, 2*PI*fp[i]/fs);
  }
  plotData(N/5, fp, mp);
}

void phonoFilterModelPlot()
{
  //double fs = 192000;  // sample-rate
  double fs = 44100;  // sample-rate
  //double fs = 8000;  // sample-rate

  rsPhonoFilterDD preEmphasisFilter;
  preEmphasisFilter.setSampleRate(fs);
  preEmphasisFilter.setMode(rsPhonoFilterDD::PRE_EMPHASIS);

  rsPhonoFilterDD deEmphasisFilter;
  deEmphasisFilter.setSampleRate(fs);
  deEmphasisFilter.setMode(rsPhonoFilterDD::DE_EMPHASIS);

  // compute and plot the magnitude response
  static const int N = 1000;
  double fL = 1.0;
  double fH = 1000000.0;
  double fp[N], dBp[N], dBm[N], dBmi[N]; // frequencies for the plot, dB-gains of prototype, model
                                         // and inverse model
  for(int n = 0; n < N; n++)
  {
    // frequency in Hz:
    fp[n]  = rsLinToExp((double)n, 0.0, (double)N-1, fL, fH);

    // analog prototype gain:
    dBp[n] = rsAmp2dB(rsPhonoFilterDD::prototypeMagnitudeAt(fp[n]));

    // digital model gain:
    dBm[n] = rsAmp2dB(preEmphasisFilter.getMagnitudeAt(fp[n]));
    if(fp[n] > fs/2)
      dBm[n] = 0.0; // artificially cut off the response plot

    // inverse model gain:
    dBmi[n] = rsAmp2dB(deEmphasisFilter.getMagnitudeAt(fp[n]));  // dB-gain of the model
    if(fp[n] > fs/2)
      dBmi[n] = 0.0;
  }

  plotDataLogX(N, fp, dBp, dBm, dBmi);
}

void phonoFilterSimulation()
{
  double fs = 44100;           // samplerate in Hz
  static const int N = 8192;  // number of sample to generate

  rsPhonoFilterDD preEmphasisFilter;
  preEmphasisFilter.setSampleRate(fs);
  preEmphasisFilter.setMode(rsPhonoFilterDD::PRE_EMPHASIS);

  rsPhonoFilterDD deEmphasisFilter;
  deEmphasisFilter.setSampleRate(fs);
  deEmphasisFilter.setMode(rsPhonoFilterDD::DE_EMPHASIS);

  // impules responses (pre- and de-emphasis):
  double hp[N], hd[N];

  // frequency axis and magnitude responses:
  double f[N/2], mp[N/2], md[N/2];
  for(int k = 0; k < N/2; k++)
    f[k] = k*fs/N;

  getImpulseResponse(preEmphasisFilter, hp, N);
  getImpulseResponse(deEmphasisFilter, hd, N);
  rsMagnitudeAndPhase(hp, N, mp);
  rsMagnitudeAndPhase(hd, N, md);
  RAPT::rsArrayTools::applyFunction(mp, mp, N/2, &rsAmp2dB);
  RAPT::rsArrayTools::applyFunction(md, md, N/2, &rsAmp2dB);

  plotDataLogX(N/2, f, mp, md);
}


void serialParallelBlend()
{
  // tests a continuous blend between a serial and a parallel connection of some number of filters
  // for testing, we use bell filters implemented as TPT/ZDF SVF

  double fs = 44100.0;
  double bw = 0.2;       // bandwidths
  double g  = 2.0;       // linear gains
  double f1 = 1000.0;    // frequency of 1st filter
  double gs = 0.2;       // gain for serial path
  double gp = 1.0 - gs;  // gain for parallel path
  static const int numFilters = 3;
  static const int numSamples = 2048;
  rsStateVariableFilterDD f[numFilters];

  int n, m;
  for(m = 0; m < numFilters; m++)
  {
    f[m].setSampleRate(fs);
    f[m].setMode(rsStateVariableFilterDD::BELL);
    f[m].setBandwidth(bw);
    f[m].setGain(g);
    f[m].setFrequency((m+1)*f1);
  }

  // obtain impulse-response:
  double x[numSamples], y[numSamples];
  RAPT::rsArrayTools::fillWithZeros(x, numSamples);
  x[0] = 1;
  double w, z;
  double accu;
  for(n = 0; n < numSamples; n++)
  {
    accu  = 0.0;
    z     = gp*x[n];
    w     = f[0].getSample(x[n]);
    accu += w;
    for(m = 1; m < numFilters; m++)
    {
      w     = f[m].getSample(gs*w + z);
      accu += w;
    }
    y[n] = gs*w + gp*accu/numFilters;
  }

  // obtain and plot magnitude response:
  double frq[numSamples/2], mag[numSamples/2], dB[numSamples/2];
  for(int k = 0; k < numSamples/2; k++)
    frq[k] = k*fs/numSamples;
  rsMagnitudeAndPhase(y, numSamples, mag);
  RAPT::rsArrayTools::applyFunction(mag, dB, numSamples/2, &rsAmp2dB);
  plotDataLogX(numSamples/2, frq, dB);

}


void averager()
{
  // averager based on a linear combination of exponential decays, that has an approximately 
  // constant (weight) value over some time and then deacys away

  static const int N = 5000;

  double tc = 1000.0; // cutoff time

  double t1 = tc;
  double t2 = tc * 0.25;
  double w1 = +1.0;
  double w2 = -1.0;

  double t3 = tc * 0.14;
  double w3 = 0.5;

  double t[N], y1[N], y2[N], y3[N], ySum[N];

  RAPT::rsArrayTools::fillWithRangeLinear(t, N, 0.0, (double)(N-1));
  RAPT::rsArrayTools::fillWithZeros(ySum, N);

  // accumulate the decay functions
  int n;
  for(n = 0; n < N; n++)
  {
    y1[n]    = w1 * exp(-n/t1);
    ySum[n] += y1[n];
  }
  for(n = 0; n < N; n++)
  {
    y2[n]    = w2 * exp(-n/t2);
    ySum[n] += y2[n];
  }
  for(n = 0; n < N; n++)
  {
    y3[n]    = w3 * exp(-n/t3);
    ySum[n] += y3[n];
  }


  plotData(N, t, ySum, y3);
  // maybe the remaining bump can be compesnated by another attach/decay function


  int dummy = 0;
}


void movingAverage()
{
  static const int L = 1000;
  static const int N = 2*L;

  double t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsMovingAverageDD ma;
  //ma.setLengthInSamples(L);
  ma.setLengthInSeconds(0.02);
  //ma.setDeviation(0.001); // good for practical use
  ma.setDeviation(0.01);
  //ma.setDeviation(0.1);   // good for showing the leakage in a plot
  //ma.setDeviation(0.0); 


  getImpulseResponse(ma, h, N);
  plotData(N, t, h);
  double sum = RAPT::rsArrayTools::sum(h, L);
}

void trapezAverager()
{
  static const int L1 = 10;
  static const int L2 = 30;
  static const int N = L1+L2+10;

  double t[N], h[N];
  RAPT::rsArrayTools::fillWithIndex(t, N);
  rsMovingAverageDD ma1, ma2;
  ma1.setLengthInSamples(L1);
  ma2.setLengthInSamples(L2);

  getImpulseResponse(ma1, h, N);
  for(int n = 0; n < N; n++)
    h[n] = ma2.getSample(h[n]);
  plotData(N, t, h);

  //double sum = rsSum(h, N);
}

void compareApproximationMethods()
{
  // filter parameters:
  //double fs    = 96000.0;  // samplerate
  double fs    = 44100.0;  // samplerate
  double fc    =  1000.0;  // cutoff frequency
  double Ap    =     1.0;  // passband ripple in dB
  double As    =    50.0;  // stopband rejection in dB
  int    order =    5;    // prototype filter order

  typedef RAPT::rsInfiniteImpulseResponseDesigner<double> IIRD;
  typedef RAPT::rsPrototypeDesigner<double> PTD;



  // create and set up the filter:
  rsEngineersFilterMono flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  //flt.setMode(IIRD::LOWPASS);
  //flt.setMode(IIRD::HIGHPASS);
  flt.setMode(IIRD::BANDPASS);
  //flt.setMode(IIRD::BANDREJECT);
  flt.setRipple(Ap);
  flt.setStopbandRejection(As);
  flt.setPrototypeOrder(order);

  // create impulse responses of different types of filters:
  static const int N = 7000;    // number of samples in impulse- and magnitude response
  double hBes[N], hBut[N], hCheb1[N], hCheb2[N], hEll[N], hPap[N], hHalp[N], hGaus[N];  // impulse responses
  double mBes[N], mBut[N], mCheb1[N], mCheb2[N], mEll[N], mPap[N], mHalp[N], mGaus[N];  // magnitude responses
  double f[N];
  //rsFillWithRangeLinear(f, N, 0.0, fs/2);
  RAPT::rsArrayTools::fillWithRangeExponential(f, N, 1.0, fs/2);

  flt.setApproximationMethod(PTD::BESSEL);
  getImpulseResponse(flt, hBes, N);
  flt.getMagnitudeResponse(f, mBes, N, true);

  flt.setApproximationMethod(PTD::BUTTERWORTH);
  getImpulseResponse(flt, hBut, N);
  flt.getMagnitudeResponse(f, mBut, N, true);

  flt.setApproximationMethod(PTD::CHEBYCHEV);
  getImpulseResponse(flt, hCheb1, N);
  flt.getMagnitudeResponse(f, mCheb1, N, true);

  flt.setApproximationMethod(PTD::INVERSE_CHEBYCHEV);
  getImpulseResponse(flt, hCheb2, N);
  flt.getMagnitudeResponse(f, mCheb2, N, true);

  flt.setApproximationMethod(PTD::ELLIPTIC);
  getImpulseResponse(flt, hEll, N);
  flt.getMagnitudeResponse(f, mEll, N, true);

  flt.setApproximationMethod(PTD::PAPOULIS);
  getImpulseResponse(flt, hPap, N);
  flt.getMagnitudeResponse(f, mPap, N, true);

  flt.setApproximationMethod(PTD::HALPERN);
  getImpulseResponse(flt, hHalp, N);
  flt.getMagnitudeResponse(f, mHalp, N, true);

  flt.setApproximationMethod(PTD::GAUSSIAN);
  getImpulseResponse(flt, hGaus, N);
  flt.getMagnitudeResponse(f, mGaus, N, true);


  // plots:
  //plotData(N, 0, 1/fs, hBes, hBut, hCheb2, hCheb1, hEll);
  //plotData(N/10, 0, 1/fs, hBut, hCheb2);
  //plotData(N, 0, 1/fs, hBut, hPap, hCheb1);
  //plotDataLogX(N, f, mBes, mBut, mCheb2, mCheb1, mEll);
  plotDataLogX(N, f, mBut, mPap, mHalp);   // Butterworth, Papoulis, Halpern
  plotDataLogX(N, f, mBut, mBes, mGaus);     // Butterworth, Bessel, Gaussian - gauss too wide - use asymptotic normalization
  //plotDataLogX(N, f, mBut, mPap, mCheb1);  // Butterworth, Papoulis, Chebychev1
  //plotDataLogX(N, f, mBut, mPap, mEll);
  //plotDataLogX(N, f, mBut, mCheb2);

  // for the asymptotic normalization (for Gaussian), we take the Butterworth frequency response as
  // reference - the mag-squared response goes like: 
  //   |H(w)|^2 = 1 / (1 + w^(2*N))
  // and that's what we want the asymptotic response of all the other allpole filters also look 
  // like. In general, for a transfer functions given by 
  //   H(s) = (sum_i=0^M b_i s^i) / (sum_j=0^N a_i s^i)
  // we see an aymptotic behavior of the mag-squared response of
  //   |A(w)|^2 = (b_M * w^M)^2 / (a_N * w^N)^2 ...verify this
  // and in case of allpole filters just:
  //   |A(w)|^2 = 1 / (a_N * w^N)^2
  // so, to get the same asymptotic behavior, we need to scale the poles by a_N 
  // ....or a_N^2 or sqrt(a_N)
  //
  // i think, the normalization should be implemented in rsPrototypeDesigner<T>::zpkFromTransferCoeffsLP
  // maybe we should just multiply all poles (and zeros, if any), by a factor b_M / a_M and scale 
  // the gain by the reciprocal of that value? maybe intorduce functions 
  // normalizeAsymptote(poles, numPoles, zeros, numZeros, b, a) ...maybe that should call a 
  // getAsymptoticNormalizer(b, Nb, a, Na)

  // ...actually, that normalization may make sense for Papoulis and Halpern too - they are 
  // asymtotocially too narrow compared to Butterworth - the same is true for Chebychev1

  /*
  // normalize impulse repsonses and write to wavefiles (maybe remove)
  double s = 1.0;
  //s = 1.0 / rsMaxAbs(hBut, N); // Butterworth impulse response has highest peak
  RAPT::rsArrayTools::scale(hBes, N, s);
  RAPT::rsArrayTools::scale(hBut, N, s);
  RAPT::rsArrayTools::scale(hCheb1, N, s);
  RAPT::rsArrayTools::scale(hCheb2, N, s);
  RAPT::rsArrayTools::scale(hEll, N, s);
  RAPT::rsArrayTools::scale(hPap, N, s);
  writeToMonoWaveFile("Bessel.wav", hBes, N, (int)fs, 16);
  writeToMonoWaveFile("Butterworth.wav", hBut, N, (int)fs, 16);
  writeToMonoWaveFile("Chebychev1.wav", hCheb1, N, (int)fs, 16);
  writeToMonoWaveFile("Chebychev2.wav", hCheb2, N, (int)fs, 16);
  writeToMonoWaveFile("Elliptic.wav", hEll, N, (int)fs, 16);
  writeToMonoWaveFile("Papoulis.wav", hPap, N, (int)fs, 16);
  */

  // BUG: bandpass and bandreject don't work for chebychev2/elliptic (types with finite zeros)
  // ...must be in LP->BP, LP->BR transform ...check old code, write a test in the old codebase
  // with the same parameters, looks like some zeros are missing(?)

  // Observations: 
  // -Butterworth has highest peak value
  // -Chebychev 2 impulse response is similar to Butterworth
  //  -increasing stopband rejection...
  // -Chebychev 1 impulse response is similar to elliptic (for low orders)
  //  -increasing passband ripple increases the ringing time (for both)
  // -Papoulis impulse response is in between Butterworth and Chebychev1
  //  -more or less halfway, if the passband-ripple is around 3dB, in this case, the magnitude
  //   in the stopband is also halfway between
  // -the ringing time seems to correlate well with the slope at the cutoff frequency (check this)

  // todo: see, what happens, if we increase the order

  // the impulse reponses can sound interesting as transients for synthetic drums - highpass 
  // responses may make for a good "click" at the start of a bassdrum
  // maybe "zaps" are also good for that (maybe they can be created with allpasses)
  // make an electronic drum-synthesizer "FilterDrum" or "Impulse"

  // interesting read:
  // http://en.wikipedia.org/wiki/Elliptic_filter
  // http://en.wikipedia.org/wiki/Elliptic_rational_functions
}

void compareOldAndNewEngineersFilter()
{
  // obsolete - the old version has been removed

  double fs    = 44100.0;  // samplerate
  double fc    =  1000.0;  // cutoff frequency
  double Ap    =     1.0;  // passband ripple in dB
  double As    =    50.0;  // stopband rejection in dB
  int    order =    6;     // prototype filter order

  typedef RAPT::rsInfiniteImpulseResponseDesigner<double> IIRD;
  typedef RAPT::rsPrototypeDesigner<double> PTD;

  rosic::rsEngineersFilterMono efOld;
  efOld.setPrototypeOrder(order);
  efOld.setSampleRate(fs);
  efOld.setFrequency(fc);
  efOld.setRipple(Ap);
  efOld.setStopbandRejection(As);
  efOld.setApproximationMethod(PTD::ELLIPTIC);
  efOld.setMode(IIRD::BANDPASS);

  rosic::rsEngineersFilterMono  efNew;
  efNew.setPrototypeOrder(order);
  efNew.setSampleRate(fs);
  efNew.setFrequency(fc);
  efNew.setRipple(Ap);
  efNew.setStopbandRejection(As);
  efNew.setApproximationMethod(PTD::ELLIPTIC);
  efNew.setMode(IIRD::BANDPASS);

  // in rsInfiniteImpulseResponseDesigner<T>::getBiquadCascadeCoefficients, the zeros are not 
  // properly sorted

  int dummy = 0;
}

void testPoleZeroMapper()
{
  // some debug/throwaway code to figure out why the new RAPT elliptic bandpass zeros end up wrong

  //// new version:
  //typedef rsPoleZeroMapperD PZM;     
  //typedef std::complex<double> cmplx;

  // old version:
  typedef rosic::rsPoleZeroMapperD PZM; 
  typedef std::complex<double> cmplx;

  // create some pole/zero arrays:
  cmplx c1(1.,1.), c2(1.,-1.), c3(2.,2.), c4(2.,-2.), z(0.,0.);
  cmplx poles[8] = { c1, c2, c3, c4, z, z, z, z }; // 4 additional zeros because of order doubling
  cmplx zeros[8] = { c1, c2, c3, c4, z, z, z, z };
  cmplx poles2[8], zeros2[8]; // targets for non-in-place processing
  double g = 1;

  //PZM pzm;
  //PZM::sLowpassToBandpass(zeros, poles, &g, zeros2, poles2, &g, 4, 0.5, 2.0); // old/new identical
  PZM::sPlanePrototypeToBandpass(poles, zeros, poles2, zeros2, 4, 0.5, 2.0); // old/new identical but different from sLowpassToBandpass
  //PZM::prototypeToAnalogBandpass(poles, 4, zeros, 4, &g, 0.5, 2.0); // old/new identical


  PZM::bilinearAnalogToDigital(poles2, 8, zeros2, 8, 5., &g); // old/new identical

  // what's the difference between prototypeToAnalogBandpass and sLowpassToBandpass?  maybe the 
  // in-place processing sounds like they should do the same thing? are they redundant?
  // ...if so, get rid of redundant code
  // and what about sPlanePrototypeToBandpass - this seems to be actually used in InfiniteImp..

  int dummy = 0;
}

void ringingTime()
{
  typedef RAPT::rsInfiniteImpulseResponseDesigner<double> IIRD;
  typedef RAPT::rsPrototypeDesigner<double> PTD;


  // filter parameters:
  double fs        = 44100.0;  // samplerate
  double fc        = fs/16;
  double Ap        =     1.0;  // passband ripple in dB
  double As        =    50.0;  // stopband rejection in dB
  double threshold =     0.01; // threshold amplitude
  int    order     =     5;    // filter order
  //int    method    = rsPrototypeDesigner::BUTTERWORTH; // approximation method
  int    method    = PTD::ELLIPTIC;

  // create and set up the filter:
  rsEngineersFilterMono flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setMode(IIRD::LOWPASS);
  flt.setRipple(Ap);
  flt.setStopbandRejection(As);
  //flt.setOrder(order);
  flt.setPrototypeOrder(order);
  flt.setApproximationMethod(method);

  // compute ringing time:
  double rt = flt.getRingingTimeEstimate(threshold);

  // create impulse and normalize response:
  static const int N = 1000;    // number of samples
  double h[N];
  getImpulseResponse(flt, h, N);
  double s = 1.0 / RAPT::rsArrayTools::maxAbs(h, N);
  RAPT::rsArrayTools::scale(h, N, s);

  // plot:
  plotData(N, 0, 1, h);
}

void butterworthSquaredLowHighSum()
{
  // We plot the Butterworth squared magnitude response (that may result from applying a 
  // Butterworth filter bidirectionally) for the lowpass and highpass prototype and the sum of
  // both. 

  static const int P = 1000;  // number of datapoints for the plot
  double wMin = 0.01;
  double wMax = 100;

  int N = 5;        // order of Butterworth filter
  int M = 3;        // number of (bidirectional) passes

  double eps = pow(0.5, -1.0/M) - 1; // follows from fixing the sum to unity at w=1

  double w[P];      // radian frequencies
  double H2L[P];    // |H(w)|^2 for the lowpass
  double H2H[P];    // |H(w)|^2 for the highpass
  double H2S[P];    // |H(w)|^2 for the sum of low- and highpass

  RAPT::rsArrayTools::fillWithRangeExponential(w, P, wMin, wMax);
  for(int i = 0; i < P; i++)
  {
    H2L[i] = pow(1.0 / (1.0 + eps*pow(w[i]*w[i], N)), M);
    H2H[i] = pow(1.0 / (1.0 + eps*pow(1.0/(w[i]*w[i]), N)), M);
    H2S[i] = H2L[i] + H2H[i];
  }

  plotDataLogX(P, w, H2L, H2H, H2S);

  // Observations:
  // for M=1, N arbitrary: leads to a unity sum
  // for other M, we see a wiggle around the cutoff frequency which increases moderately with M

  // Conclusion: 
  // unity sum filtering works only for a single bidirectional pass of the Butterworth filter
}

void gaussianPrototype()
{
  int N = 5;  // 1...10
  double c[21];

  double wc = 0.5;                 // cutoff frequency
  double g  = log(2.0) / (wc*wc);  // "gamma"


  typedef RAPT::rsPrototypeDesigner<double> PTD;
  PTD::gaussianPolynomial(c, N, wc);

  static const int numBins = 1000;
  double wMin = 0.0;
  double wMax = 3.0;
  double w[numBins], m[numBins], mt[numBins];
  RAPT::rsArrayTools::fillWithRangeLinear(w, numBins, wMin, wMax);
  for(int i = 0; i < numBins; i++)
  {
    m[i]  = 1.0 / (RAPT::rsPolynomial<double>::evaluate(w[i], c, 2*N));
    mt[i] = exp(-g*w[i]*w[i]);  // ideal Gaussian response
  }

  plotData(numBins, w, mt, m);
}


void halpernPrototype()
{
  // compare Halpern and Papoulis prototype magnitude responses

  double eps = 1.0;

  double aH[50], aP[50];

  int N = 10;

  typedef RAPT::rsPrototypeDesigner<double> PTD;
  PTD::halpernPolynomial(aH, N);
  PTD::papoulisPolynomial(aP, N);


  //maximumSlopeMonotonicPolynomial(aP, N);
    // this function uses a different convention for the coeffs - they act on powers of w^2 instead
    // of powers of w - maybe use the same convention later for the Halpern polynomial, too

  // we have a sign alternation in our numbers which is not present in Paarmann's MatLab code
  // ...why? ...because we compute coeffs of H(w^2), the Paarmann code for H(s)*H(-s) - setting
  // s = j*w in H(s)*H(-s) leads to this sign alternation


  // evaluate magnitude-squared function and plot:
  double wMin = 0.01;
  double wMax = 100;
  static const int numBins = 1000;
  double w[numBins], mH[numBins], mP[numBins];
  RAPT::rsArrayTools::fillWithRangeExponential(w, numBins, wMin, wMax);
  for(int i = 0; i < numBins; i++)
  {
    mH[i] = 1.0 / (1.0 + eps*eps * RAPT::rsPolynomial<double>::evaluate(w[i], aH, 2*N));
    mP[i] = 1.0 / (1.0 + eps*eps * RAPT::rsPolynomial<double>::evaluate(w[i], aP, 2*N));

    //mP[i] = 1.0 / (1.0 + eps*eps * evaluatePolynomialAt(w[i]*w[i], aP, N));

    mH[i] = rsAmp2dB(sqrt(mH[i]));
    mP[i] = rsAmp2dB(sqrt(mP[i]));

  }
  plotDataLogX(numBins, w, mH, mP);

  int dummy = 0;
}


void plotMaxSteepResponse(double *p, int M, double *q, int N, double k)
{
  // a function to plot a magnitude squared reponse together with a Butterworth response as
  // reference for evaluating out max-steep designs

  // evaluate magnitude-squared function and plot:
  double wMin = 0.01;
  double wMax = 100;
  static const int numBins = 1000;
  double w[numBins], r[numBins], rB[numBins];
  RAPT::rsArrayTools::fillWithRangeExponential(w, numBins, wMin, wMax);
  int n;
  double num, den, x;
  for(n = 0; n < numBins; n++)
  {
    x = w[n]*w[n];

    // evaluate our filter's response:
    num   = RAPT::rsPolynomial<double>::evaluate(x, p, M);
    den   = RAPT::rsPolynomial<double>::evaluate(x, q, N);
    r[n]  = k * num / den;

    // evaluate the reference Butterworth response:
    rB[n] = 1 / (1 + pow(x, N));
  }

  // plot:
  plotDataLogX(numBins, w, rB, r);
}
void maxFlatMaxSteepPrototypeM1N2()
{
  // Consider an analog (s-domain) magnitude squared frequency response function:
  //
  //             p0 + p1*w^2 + p2*w^4 + ... + pM*w^2M
  // |H(w)|^2 = --------------------------------------
  //             q0 + q1*w^2 + q2*w^4 + ... + qN*w^2N
  //
  // where N >= M. It is desired to arrange the zeros in a way so as to have a maximally flat
  // response in the passband, a maximum slope at w = 1 and a monotonically decreasing response.
  // The Butterworth filter is maximally flat but has no finite zeros. Maybe it's possible to
  // choose poles and zeros that retain the maximum flatness and steepen the slope at the cutoff
  // frequency. For a monotonic response, these zeros should not be on the imaginary axis (as is
  // the case in Chebychev2 and elliptic filters).
  // To simplify matters, we replace w^2 by x and call the function r(x) and normalize p0, q0 to
  // unity and introduce an overall gain factor k instead, so:
  //
  //             1 + p1*x + p2*x^2 + ... + pM*x^M
  // r(x) = k * ----------------------------------
  //             1 + q1*x + q2*x^2 + ... + qN*x^N
  //
  // we choose M = N-1 and as a first example let M = 1, N = 2 such that:
  //
  //             1 + p1*x              p(x) 
  // r(x) = k * ------------------- = ------
  //             1 + q1*x + q2*x^2     q(x)
  //
  // we have 4 parameters: k, p1, q1, q2 and impose 4 constraints:
  //
  // r(x=0) = 1, r(x=1) = g, r'(x=0) = 0, r'(x=1) = max -> r''(x=1) = 0
  //
  // r(0) = 1 fixes k = 1
  // r(1) = g: fixes q2: g = (1+p1)/(1+q1+q2) -> q2 = (1+p1)/g - (1+q1)
  // the first derivative is:
  //
  //          p'(x)*q(x) - q'(x)*p(x)     p1*(1+q1*x+q2*x^2) - (q1+2*q2*x)*(1+p1*x)
  // r'(x) = ------------------------- = -------------------------------------------
  //                 p(x)*p(x)                      (1+p1*x) * (1+p1*x)
  //
  // r'(x=0) = 0 gives: (p1*(1+0+0) - (q1+0)*(1+0)) / ((1+0)*(1+0)) -> (p1-q1)/1 = 0
  // so p1 = q1
  //
  // We have now only one free parameter left: p1 - all others can be computed from the conditions.
  // Trying various values, we see that for p1 < 0, the magnitude response undershoots the x-axis
  // and negative magnitudes are not permitted. For p1 > 0, the response is less flat in the 
  // passband and less steep at the cutoff than the Butterworth response. Conclusion: p1 = 0 must
  // be chosen and this reduces to the Butterworth response again. So, a real zero seems not to be 
  // able to improve a 2nd order Butterworth response. Maybe we need a pair of complex conjugate
  // zeros a + j*b, a - jb where b is slightly above unity and a is a small negative number?
  // This should be able steepen the slope at unity at the expense of having the response not
  // approach 0 as x approaches infinity. If we want this, we would need to make the denominator
  // 3rd order.

  static const int M = 1;       // numerator order
  static const int N = 2;       // denominator order
  double p[M+1];                // numerator coeffs
  double q[N+1];                // denominator coeffs
  double g = 0.5;               // magnitude squared at x = 1
  double k = 1;                 // from r(x=0) = 1
  p[0] = 1;                     // normalized
  p[1] = 0.0;                   // trial and error
  q[0] = 1;                     // normalized
  q[1] = p[1];                  // from r'(x=0) = 0
  q[2] = (1+p[1])/g - (1+q[1]); // from r(x=1) = g

  plotMaxSteepResponse(p, M, q, N, k);
}
void maxFlatMaxSteepPrototypeM2N2()
{
  //         1 + p1*x + p2*x^2     p(x) 
  // r(x) = ------------------- = ------   // the gain k is always equal to 1, so it's left out
  //         1 + q1*x + q2*x^2     q(x)

  // r'(x) = (p'(x)*q(x) - q'(x)*p(x)) / (p(x)*p(x))  quotient rule
  //
  // p'(x) = p1 + 2*p2*x, q'(x) = q1 + 2*q2*x
  //
  //          (p1+2*p2*x)*(1+q1*x+q2*x^2) - (q1+2*q2*x)*(1+p1*x+p2*x^2)      s(x)
  // r'(x) = ----------------------------------------------------------- =: ------
  //                      (1+p1*x+p2*x^2) * (1+p1*x+p2*x^2)                  t(x)
  //
  // r'(0) = 0 = p1 - q1 -> p1 = q1
  //
  // s(x) := (p1+2*p2*x)*(1+q1*x+q2*x^2) - (q1+2*q2*x)*(1+p1*x+p2*x^2)
  // t(x) := (1+p1*x+p2*x^2) * (1+p1*x+p2*x^2)
  //
  // r''(x) = (s'(x)*t(x) - t'(x)*s(x)) / (s(x)*s(x))  quotient rule
  //
  // s'(x) = 
  // t'(x) =
  // maybe use maxima - it's messy to do by hand

  static const int M = 2;
  static const int N = 2;
  double p[M+1];                // numerator coeffs
  double q[N+1];                // denominator coeffs
  double g = 0.5;
  double k = 1.0;

  p[0] = 1;                          // normalized
  p[1] = 0.2;                        // trial and error
  p[2] = 0.4;
  q[0] = 1;                          // normalized
  q[1] = p[1];                       // from r'(x=0) = 0
  q[2] = (1+p[1]+p[2])/g - (1+q[1]); // from r(x=1) = g

  // in general qN = (sum_{m=0}^M p_m / g) - sum_{n=0}^{N-1} q_n, in code:
  // q[N] = rsSum(p, M+1) / g - rsSum(q, N);

  plotMaxSteepResponse(p, M, q, N, k);

  // or - maybe alternatively to impose conditions on derivatives at x=1 (that they should be 
  // maximal), use the additional degrees of freedom to make the slope at DC even more flat
  // maybe, for the filter to be lowpass, we need N > M, then (N-M)*6 gives the slope in dB/oct
  // furthermore, it seems necessarry that M is even, because real zeros in r(x) don't seem
  // to be of value (we need the to come in complex conjugate pairs)
  // maybe the responses can be seen as to decouple the flatness from the asymptotic slope
  // with Butterworth filters, choosing higher order makes the passband flatter and the slope
  // steeper. here, L=N-M controls the slope and M controls the flatness. by keeping L constant
  // and increasing M, the slope remains the same, but the response approximates the straight-line
  // idealization (flat until cutoff, straight downward line from cutoff) better

  // In the derivation of the Butterworth response, a polynomial long division was made on the 
  // magintude-squred function (this is a deconvolution of the p, q coefficient arrays) (this
  // seems to give rise to an infinite sequence - is this right?). Then a Taylor expansion was 
  // of the magnitude-squared function was equated to this and derivatives at 0 were set to zero.
  // see rsDeConvolve and PowerSeries.pdf

  // maybe these filters can be designed by using rational interpolation with a number of
  // prescribed derivatives at x=0 and x=1
}

void envelopeFilter1(double *x, double *yEnv, double *yMod, int N, double fc, double fs,
  int numPasses = 8)
{
  double mean = RAPT::rsArrayTools::mean(x, N);     // compute mean
  RAPT::rsArrayTools::add(x, -mean, yEnv, N);       // subtract it from x
  rsBiDirectionalFilter::applyButterworthLowpass(yEnv, yEnv, N, fc, fs, 1, numPasses);
  RAPT::rsArrayTools::add(yEnv, mean, yEnv, N);     // undo mean subtraction
  RAPT::rsArrayTools::subtract(x, yEnv, yMod, N);   // compute residual
}

void envelopeFilter2(double *x, double *yEnv, double *yMod, int N, double fc, double fs,
  int numPasses = 8)
{
  double mean = RAPT::rsArrayTools::mean(x, N);      // compute mean
  RAPT::rsArrayTools::add(x, -mean, yEnv, N);        // subtract it from x
  RAPT::rsArrayTools::difference(yEnv, N);           // take difference signal
  rsBiDirectionalFilter::applyButterworthLowpass(yEnv, yEnv, N, fc, fs, 1, numPasses);
  RAPT::rsArrayTools::cumulativeSum(yEnv, yEnv, N);  // undo difference
  RAPT::rsArrayTools::add(yEnv, mean, yEnv, N);      // undo mean subtraction
  RAPT::rsArrayTools::subtract(x, yEnv, yMod, N);    // compute residual

  // transfer mean of residual signal into envelope signal:
  mean = RAPT::rsArrayTools::mean(yMod, N);
  RAPT::rsArrayTools::add(yEnv, +mean, yEnv, N);
  RAPT::rsArrayTools::add(yMod, -mean, yMod, N);
}

void envelopeFilter3(double *x, double *yEnv, double *yMod, int N, double fc, double fs,
  int numPasses = 8)
{
  // use the filter on the normal signal
  double *tmp = new double[N];
  envelopeFilter2(x, tmp, yMod, N, fc, fs, numPasses);

  // use the filter on the reversed signal (then undo reversal):
  RAPT::rsArrayTools::reverse(x, N);
  envelopeFilter2(x, yEnv, yMod, N, fc, fs, numPasses);
  RAPT::rsArrayTools::reverse(x, N);
  RAPT::rsArrayTools::reverse(yEnv, N);

  // crossfade:
  for(int n = 0; n < N; n++)
  {
    double k = n / double(N-1);
    yEnv[n] = k*tmp[n] + (1-k)*yEnv[n];
  }

  RAPT::rsArrayTools::subtract(x, yEnv, yMod, N);    // compute residual
  delete[] tmp;
}

void splitLowFreqFromDC()
{
  // We try to split a low frequency oscillation from a DC component using a bidirectional lowpass
  // filter.

  static const int N = 1000;     // number of samples
  double fs          = 44100.0;  // samplerate
  double f           =   400.0;  // sine frequency
  double a           =     0.2;  // sine amplitude
  double dc1         =     1.0;  // amount of DC at start
  double dc2         =     1.5;  // DC at end
  double fc          =   100.0;  // filter cutoff frequency
  double order       =     1;    // order of the single-pass filter
  int    numPass     =     8;    // number of passes


  double *x = new double[N];
  double *yl = new double[N];
  double *yh = new double[N];

  // create the input signal:
  createSineWave(x, N, f, a, fs);
  RAPT::rsArrayTools::fillWithRangeLinear(yl, N, 1.0, 1.5);
  for(int n = 0; n < N; n++)
    x[n] += yl[n];

  //envelopeFilter1(x, yl, yh, N, fc, fs, numPass);
  //envelopeFilter2(x, yl, yh, N, fc, fs, numPass);
  envelopeFilter3(x, yl, yh, N, fc, fs, numPass);

  // plot:
  GNUPlotter plt;
  //plt.setLegends("Input", "Lowpass", "Highpass");
  plt.plotArrays(N, x, yl, yh);

  // todo: make a version of plotFunctionTables that don't need an explicit x-axis
  delete[] x;
  delete[] yl;
  delete[] yh;

  int dummy = 0;
}

template<class T>
void plotFreqRespZ(const std::vector<T>& b, const std::vector<T>& a, int numSamples = 501)
{
  FilterPlotter<T> plt;
  plt.addFilterSpecificationBA((int)b.size(), &b[0], (int)a.size(), &a[0], 1.0);
  plt.plotFrequencyResponses(numSamples, 0.0, 0.5, false);
  //plt.plotMagnitude(numSamples, 0.0, 0.5, false, false);

  // todo: plotMagnitudeAndPhase, take optional parameter for number of samples - should behave
  // similar to matlab's freqz, move to plotting: 
  // https://www.mathworks.com/help/signal/ref/freqz.html
}

void directFormFreqResp()
{
  // We plot magnitude- and phase responses of a direct form filter

  using Real = double;
  using Vec  = std::vector<Real>;

  Vec b({0.1});        // feedforward coeffs
  Vec a({1.0, -0.9});  // feedback coeffs
  plotFreqRespZ(b, a);  // We get a DC gain of 10 - why? Shouldn't it be 1?
}
// todo: plot freq-responses for enigineers filter vor various orders


void ladderResonanceGain()
{
  // An experiment to figure out the appropriate s-values (scalers for k in the compensation gain
  // formula - this scaler should depend on the mode). The goal is that, at high resonance, all
  // modes produce the same resonance gain. The 4-pole lowpass is taken as reference. It uses s=1.
  // Currently, i'm using some rough values that "look nice" i.e. are simple, round numbers and
  // work as rough approximation - but we want a better match...

  using Real = double;
  using LDR  = RAPT::rsLadderFilter<Real, Real>;
  using Vec  = std::vector<Real>;
  using Mode = LDR::Mode;

  Real fs       = 44100;   // sample rate
  Real fc       = 8000;    // cutoff frequency
  Real r        = 0.99;    // resonance (should be high)
  Real b1       = 0.23;    // 0.0: zeroless, 0.5: bilinear, 0.23: good compromise
  bool withGain = false;   // switch from plotting with or without compensation gain applied

  LDR ldr;
  ldr.setSampleRate(fs);
  ldr.setCutoff(fc);
  ldr.setResonance(r);
  ldr.setB1(b1);

  int numFreqs = 1001;
  Real fLo   = 10;
  Real fHi   = fs/2;
  Vec  freqs = rsRangeExponential(fLo, fHi, numFreqs);
  Vec  gains(numFreqs);

  std::vector<Mode> modes({ Mode::LP_24, Mode::LP_18, Mode::LP_12, Mode::LP_6, Mode::FLAT });
  //std::vector<Mode> modes({ Mode::HP_24, Mode::HP_18, Mode::HP_12, Mode::HP_6, Mode::FLAT });

  GNUPlotter plt;
  for(size_t i = 0; i < modes.size(); i++) 
  {
    ldr.setMode(modes[i]);

    // Resonance gains as function of cutoff freq:
    for(int k = 0; k < numFreqs; k++) 
    {
      ldr.setCutoff(freqs[k]);
      gains[k] = ldr.getMagnitudeResponseAt(freqs[k], withGain);
      gains[k] = rsAmp2dB(gains[k]); 
    }
    plt.addDataArrays(numFreqs, &freqs[0], & gains[0]); 

    /*
    // DC gains as function of cutoff freq:
    for(int k = 0; k < numFreqs; k++) 
    {
      ldr.setCutoff(freqs[k]);
      gains[k] = ldr.getMagnitudeResponseAt(0.0, withGain);
      gains[k] = rsAmp2dB(gains[k]); 
    }
    plt.addDataArrays(numFreqs, &freqs[0], & gains[0]); 
    */

  }
  plt.setLogScale("x");
  plt.plot();
  int dummy = 0;
  

  // Observations:
  // -When withGain = false, the curves of the 4 lowpasses have qualitatively the same shape,
  //  the LP_24 curve is the lowest and they seem to be rough equidistant for low frequencies
  //  on a dB scale. The gain goes down towards high frequencies and the effect is most severe
  //  for LP_24. For FLAT, it doesn't go down at all. For the HP types, it's the same.
  // -The DC gains have a similar shape, but they are all equal for all types. They seem to start
  //  around -13.9 dB for low cutoffs and go down to around -25.4 dB for high frequencies but 
  //  there's a spike at the Nyquist freq: it goes back up to -22.3. sooo - does that mean, k 
  //  (and therefor s) should depend on the cutoff setting and the current k is correct only
  //  for low cutoffs? ..with bilinear, it actually goes up: from -13.9 to -6
  //  -> does that perhaps indicate that we should use something in between regular and bilinear?
  //  maybe: G1(z) = b*((1-c) + c/z) / (1 + a/z) for some pre-selected c? the regular filter
  //  has c=0, the bilinear has c=0.5. To tune this correctly, assume we have chosen b0, b1. To 
  //  find the 1-pole coeff a1, we need to solve: 
  //    t = ((a1*b0-b1)*s) / (b0 + a1*b1 + (b1+a1*b0)*c)
  //  where
  //    t = tan((wc-pi)/4), s = sin(wc), c = cos(wc)
  //  as always. This gives:
  //    a1 = (b0 t + b1 c t + b1 s)/(-b0 c t + b0 s - b1 t) 
  //  regular:  a1 = t / (s-c*t)                   b0 = 1, b1 = 0
  //  bilinear: a1 = (c*t+s+t) / (s-(c+1)*t)       b0 = b1 = 1 (or 0.5 or - doesnt matter)
  //  general:  a1 = (b0 t + b1 c t + b1 s) / (-b0 c t + b0 s - b1 t) 
  // -comparing the plots for regular and bilinear, we see that they deviate from the ideal 
  //  constant in opposite ways, so it may indeed make sense to strike a compromise between the
  //  two designs. We can tune it to some optimal default value and give it as additional control
  //  to the user


  // ToDo:
  // -maybe plot not only the resonance gains but also the DC gains
  // -implement the difference equation for the compromise design and make listening tests
  //  -> provide a user parameter for b1 in the APE script
  //  ...then clean up the code
  // -try to work out the transfer function canonical from for the general case and plot the 
  //  magnitude- and pahse responses for various b1 values
  // -get rid of the setBilinear function

  // maybe make the resulting s-values static const members of the class...or just static local
  // constants in the function that needs them
}

void ladderTransferFunction()
{
  // We test the functionality of class rsLadder to return its transfer function in canoncial
  // BA form

  using Real = double;
  using LDR  = RAPT::rsLadderFilter<Real, Real>;
  using RF   = RAPT::rsRationalFunction<Real>;
  using BA   = RAPT::rsFilterSpecificationBA<Real>;
  using Mode = LDR::Mode;

  Real fs  = 44100;   // sample rate
  Real fc  = 1000;    // cutoff frequency
  bool plotPhase = true;

  LDR ldr;
  ldr.setSampleRate(fs);
  ldr.setCutoff(fc);
  //ldr.setResonance(res);
  //ldr.setMode(Mode::LP_24);  // the basic "Moog" configuration
  ldr.setMode(Mode::LP_18); 
  //ldr.setMode(Mode::HP_24);
  //ldr.setMode(Mode::FLAT);  // still problematic when reso is zero
  //ldr.setMode(Mode::BP_6_18);
  //ldr.setMode(Mode::LP_6);
  //ldr.setBilinear(true);
  ldr.setB1(0.23);


  // Plot frequency responses for resonances from 0.0 to 0.9 in 0.1 steps:
  FilterPlotter<Real> plt;
  auto addPlotWithReso = [&](Real r)
  { 
    ldr.setResonance(r); 
    plt.addTransferFunction(ldr.getTransferFunction(), fs); 
  };
  for(int i = 0; i <= 9; i++)
    addPlotWithReso(0.1 * i);
  plt.setPixelSize(800, 400);
  plt.plotFrequencyResponses(501, 31.25, 32000, true, true, true, plotPhase);
  //plt.plotPolesAndZeros(400);  // multiplicities not shown

  // Plot variuos response types with the same resonane:
  FilterPlotter<Real> plt2;

  auto addPlotWithMode = [&](Mode m)
  {
    ldr.setMode(m); 
    plt2.addTransferFunction(ldr.getTransferFunction(), fs); 
  };
  ldr.setResonance(0.99);

  //addPlotWithMode(Mode::FLAT);
  //addPlotWithMode(Mode::LP_6);
  //addPlotWithMode(Mode::LP_12);
  addPlotWithMode(Mode::LP_18);
  //addPlotWithMode(Mode::LP_24);

  /*
  addPlotWithMode(Mode::LP_24);
  addPlotWithMode(Mode::LP_18);
  addPlotWithMode(Mode::LP_12);
  //addPlotWithMode(Mode::LP_6);
  //addPlotWithMode(Mode::BP_6_18);
  addPlotWithMode(Mode::BP_6_12);
  //addPlotWithMode(Mode::BP_6_6);
  addPlotWithMode(Mode::BP_12_12);
  //addPlotWithMode(Mode::BP_18_6);
  addPlotWithMode(Mode::BP_12_6);
  //addPlotWithMode(Mode::FLAT);
  //addPlotWithMode(Mode::HP_6);
  addPlotWithMode(Mode::HP_12);
  addPlotWithMode(Mode::HP_18);
  addPlotWithMode(Mode::HP_24);
  plt2.setPixelSize(800, 400);
  */
  plt2.plotFrequencyResponses(501, 31.25, 32000, true, true, true, false);

  /*
  // test - uses getTransferFunctionAt:
  ldr.setMode(Mode::LP_24);
  ldr.setResonance(0.99); 
  plotFrequencyResponse(ldr, 501, 31.25, 32000.0, fs, true);
  ldr.setBilinear(true);
  plotFrequencyResponse(ldr, 501, 31.25, 32000.0, fs, true);
  int dummy = 0;
  // ok - that works. next: implement the transfer function computation
  */


  // Observations:
  // -Without resonance, the phase is around -172° at the cutoff frequency. This makes sense 
  //  because it does not yet contain the contribution from the unit delay. (ToDo: figure out 
  //  the formula for exact expected phase at 1 kHz and compare to the plot)
  // -Above around 3.75 kHz, the phase goes up again and goes back to 0° at the Nyquist frequency.
  // -With resonance, we actually formally get an 8-pole filter but it has pole-zero cancellations.
  // -With increasing resonance, the phase response gets steeper around fc and more flat away from
  //  fc with increasing resonance. The value at fc seems to remain fixed, though.
  // -The LP_12 mode has a phase of -86 at fc - with or without resonance.
  //  -with resonance = 0.6, the root finder does not converge - why? with LP_6, even reso = 0.5
  //   doesn't work - try to find expressions the coeffs of the transfer function via sage - maybe
  //   we run into numerical precision issues
  // -highpasses do not seem to have unit gain at high frequencies - it's a bit less, even without 
  //  resonance
  // -HP_12/6 hass all positive phase response, phase response of HP_18 starts at -90°
  //-Gain adjustment:
  // -with BP_12_12 and s=0, the response with and without resonance looks the same for high 
  //  frequencies and low frequencies are attenuated, with s=1, it looks the same for low 
  //  frequencies and high frequencies are boosted -> s=0.5 seems the most natural compromise
  //  ...same for BP_6_6
  // -The bilinear phase response is nicer overall (monotonic)
  // -the gains at resonance do not really match exactly - they differ by about 5 dB. set up a 
  //  numerical matching algorithm to find more exact values for s. ..but this should be done for 
  //  various cutoffs and resoance values - we need to see, if the values are consistent, otherwise
  //  maybe use some sort of average. the matching can be done by just measuring the resonance gain
  //  (or the peak-freq gain?) with gain compensation turned off. from this, the s-values can be 
  //  calculated: g = 1 + s*k ->  s = (g-1)/k wher we must set g := desiredGain/measuredGain
  //  and the desired gain is the gain of the 4-pole lowpass, maybe try fc = 100,1000,10000
  //  ...maybe it should be done using a really high resonance like 0.999
  //  -figure out, if we need different factors for no-zero and bilinear design - i think so.

  // ToDo: 
  // -set up the y-ranges correctly - currently, one must use the mousewheel to get the correct
  //  range
  // -FLAT still makes problems (access violations due to empty pole/zero arrays)
  // -Compare the plot of the computed frequency response with a measured frequency response to see
  //  if they match. OK - i checked against code that uses getTransferFunctionAt and the plot looks
  //  the same. Maybe to be sure, also just record the impulse response and compute the frequency
  //  reponse from that (maybe via FFT or via z-transform of the sequence)
  // -apply a nonlinear mapping to the resonance to make it more even - the graphs are crammed
  //  towards higher resonance, so we need a function that stretches the high resonances and
  //  squishes the low ones - maybe r' = r^2...although, it's mostly the rolloff slope that looks
  //  crammed - the development of the peak looks actually quite reasonable
  // -maybe offset the phase curves by -10 dB such that they don't get in the way so much
  //  but maybe not - that value is only good for lowpass
  // -plot various response types with high resonance in one plot to compare their peak gains.
  //  They should be roughly equal.
  // -maybe try a high-shelf that attenuates abobe the cutoff freq by around 10 dB...maybe the
  //  gain should be 1/(1+s*k)? because high frequencies get a boost with resonance
  // -maybe we could make the "bilinear" parameter continuous from 0 to 1. Then we would indeed
  //  have a b0 and b1 coeff
  // -getTransferFunction is really useful - this feature should be integrated into other filters
  //  too (like rsBiquadCascade, rsOnePoleFilter, rsStateVariableFilter etc.)
}

void ladderMultipole()
{
  // We try to generalize the resonance tuning formulas from here:
  // http://www.rs-met.com/documents/dsp/ResonanceTuningForTheDigitalMoogFilter.pdf
  // to a general N stage ladder. I think, we need to replace the computations of t and k in Eq 11 
  // and 15 by: t = tan((w_c - pi) / N), k = r/g_1^N. We plot the magnitude responsed of filters
  // with a given resonance settings and different cutoff frequencies

  using Real = double;
  using Vec  = std::vector<Real>;

  int numStages     = 3;               // number of ladder stages
  int numPlotFreqs  = 5001;            // number of frequency samples in the plot
  Real reso = 0.9;                                    // resonance
  Real fs = 44100;                                    // sample rate
  Vec fc({250, 500, 1000, 2000, 4000, 8000, 16000});  // cutoff frequencies
  //Vec fc({2000});

  // Computes the filter coeffs b0, a1, k for given order N, normalized radian cutoff frequency wc 
  // and resonance r:
  auto computeCoeffs = [](int N, Real wc, Real r, Real& b0, Real& a1, Real& k)
  {
    rsAssert(N >= 2, "Filter order must be at least 2 for getting a resonance"); 
    // ...because otherwise the phase response does not cross 180° but only approaches it

    Real s = sin(wc);
    Real c = cos(wc);
    Real t = tan((wc-PI)/N);

    a1 = t / (s-c*t);
    b0 = 1 + a1;

    // test:
    //Real d  = s-c*t;   // denominator
    //Real x = s/t - c;
    //Real a = 1/x;

    Real g1 = sqrt((b0*b0) / (1+a1*a1+2*a1*c));

    k = r / pow(g1, N);
  };

  // Computes the magnitude response at the given normalized radian frequency w for given filter 
  // order and coefficients:
  auto getMagnitudeAt = [](Real w, int N, Real b0, Real a1, Real k)
  {
    using Cmp = std::complex<Real>;
    Cmp j(  0, 1);                   // imaginary unit
    Cmp one(1, 0);                   // 1 + 0j
    Cmp z  = exp(j*w);               // z at which we want to evaluate |H(z)|
    Cmp G1 = b0 / (one + a1/z);      // transfer function of single stage 
    Cmp GN = pow(G1, N);             // transfer function of N stages
    Cmp H  = GN / (one + k * GN/z);  // transfer function H(z) with resonance
    return rsAbs(H);                 // magnitude of H(z)
  };

  // Compute and plot the magnitude responses:
  GNUPlotter plt;
  Vec f = rsRangeExponential(20., 20000., numPlotFreqs); 
  Vec M(numPlotFreqs);                               // magnitudes
  for(size_t i = 0; i < fc.size(); i++)
  {
    Real wc = 2 * PI * fc[i] / fs;
    Real b0, a1, k;
    computeCoeffs(numStages, wc, reso, b0, a1, k);
    for(int j = 0; j < numPlotFreqs; j++)
    {
      Real w = 2 * PI * f[j] /fs;
      M[j] = getMagnitudeAt(w, numStages, b0, a1, k);
      M[j] = rsAmpToDb(M[j]);
    }
    plt.addDataArrays(numPlotFreqs, &f[0], &M[0]);
  }
  plt.setLogScale("x", 2.0, true);
  plt.plot();

  // Observations:
  // -for N >= 3, it seems to work
  // -for N <= 2, it cannot possibly work, because the phase does not cross -180°
  // -maybe try it with a 1-pole/1-zero filter with its zero at z = -1, such that we have 
  //  b := b1 = b0, so we still have only 2 degrees of freedom: G1(z) = b(1 + 1/z) / (1 + a1/z),
  //  maybe for computing a1, we can set b=1 - it should not really matter for the phase response
  //    phs(G1) = tan(wc-pi)/N) = t = -((a0*b1-a1*b0)*s)/(a0*b0 + a1*b1 + (a0*b1+a1*b0)*c)
  //  set a0=1, b0=b1=b=1:
  //    t = -((1-a1)*s)/(1+a1 + (1+a1)*c)
  //      =  ((a1-1)*s)/(1+a1 + (1+a1)*c)
  //  wolfram: solve t=((a1-1)*s)/(1+a1 + (1+a1)*c) for a1
  //    a1 = (c t + s + t)/(s - (c + 1) t) and s != c t + t and c s + s != 0
  //
  // ToDo: 
  // -Try to derive a formula for the poles in the z-plane. They occur at location where
  //  1 + k * GN/z = 0  ->  GN + z/k = 0  ->  (b0 / (1+a1/z))^N + z/k = 0
  //  wolfram alpha can do it only when we use particular values for N, like:
  //    solve (b / (1+a/z))^2 + z/k = 0 for z
  //  ...and only up to N=4 - seems like there is no general formula. we have also
  //  b/(1+a/z) = (b*z)/(a+z) - but that doesn't seem to help
  // -figure out what happens when positive instead of negative feedback is used
  // -plot phase responses of 1-pole filters for various values of the pole coefficient - maybe the
  //  problem is that the phase response does not really pass through -180° when there are only 
  //  2 poles? ...but i actually don't think so - we'll see....
}

void ladderResonanceModeling()
{
  // rename to ladderResonanceModeling

  // We investigate the resonance of a Moog-style ladder filter by creating a step-response of
  // a filter with and without resonance, subtracting the nonresonant from the resonant response
  // to get a pure resonance signal. of particular interest is the transient of resonance signal.
  // ...later, we try to model the ladder's pure resonance signal by the impulse-response of
  // an attack/decay-sinusoid filter plus the output of some transient correction filter.

  int    N   = 2000;        // number of samples
  double fs  = 44100;       // samplerate
  double fc  = 1000;         // cutoff frequency
  double d   = 0.02;        // decay time in seconds

  // create and set up the ladder filter:
  rsLadderResoShapedDD ldr;
  ldr.setSampleRate(fs);
  ldr.setCutoff(fc);
  ldr.setResonanceDecay(d);
  ldr.setDecayByFrequency(0.0);  

  // create filter's step response:
  vector<double> x(N), yr(N);
  RAPT::rsArrayTools::fillWithValue(&x[0], N, 1.0);
  double dummy;
  for(int n = 0; n < N; n++)
    ldr.getSignalParts(x[n], &dummy, &yr[n]);

  // create a decaying sinusoid that is supposed to model the pure resonance signal:
  rsModalFilterDD dsf;
  double a  = 0.87;     // sine amplitude
  double p  = -PI/1.7;  // start phase
  double fr = 0.985*fc;  // resonance frequency (all values were set ad-hoc - we need formulas)
  dsf.setModalParameters(fr, a, 0.7*d, p, fs);
  //dsf.setModalParameters(fc, 1.0, d, 0.0, fs);
  vector<double> ym(N);
  getImpulseResponse(dsf, &ym[0], N);

  // obtain the error-signal which is the desired impulse response of the transient correction
  // filter:
  vector<double> err(N);
  RAPT::rsArrayTools::subtract(&yr[0], &ym[0], &err[0], N);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &yr[0]);
  plt.addDataArrays(N, &ym[0]);
  //plt.addDataArrays(N, &yrf[0]);
  plt.addDataArrays(N, &err[0]);
  plt.plot();

  // Observations:
  // It seems, the ladder's resonance signal is not really well modeled as a decaying sinusoid. It
  // seems to have a frequency sweep (check this) and decay-dependent center-frequency - the higher
  // the decay time, the closer the actual frequency approaches the cutoff frequency - maybe this
  // can be explained by the peak-frequency in the magnitude response. But the amplitude, phase and
  // decay for the decaying sinusoid needs also adjustments in order to find a best match. It would
  // seem that the optimal parameters for the decaying sine are all functions of the ladder's 
  // decay-time - maybe 4 tables should be created that map the normalized resonance parameter
  // (0..1) to multipliers or direct values for the decaying sine parameters:
  // r: | 0.0  0.25  0.5  0.75  1.0
  // -------------------------------
  // f: | 0.0   ??    ??   ??   1.0    (multiplier)
  // a: | 0.0                          (direct value)
  // d: |                              (multiplier)
  // p: |                              (direct value)
  // ...but it seems, the tables would have to be 2-dimensional - the values would be different
  // for each fc

  // Maybe the transient in the resonance signal can be brought down a bit by scaling the cutoff
  // frequency of the nonresonant filter ...yes - that seems to work - good values for the 
  // multiplier for the cutoff frequency of the nonresonant filter seem to be between 1.5..1.6
  // choose it according to the criterion that the resonance-signal's bandwidth should be as small
  // as possible. A value of 1.52 seems to make be around the limit such that the attack of the
  // resonance doesn't show appreciable undershoot. ...but that seems to depend of the 
  // cutoff-frequency, too - and maybe also on the resonance decay time - or maybe on the ratio
  // or product of both values? maybe we should derive a formula in terms of cutoff and decaytime

  // ToDo:
  // -figure out a formula, that computes the actual frequency at which the ladder resonates (it 
  //  seems to deviate from the cutoff frequency - possibly depending on the decay-time, because
  //  the peak in the mgnitude response of the ladder seems to be below the resonance frequency
  //  at lower resonance settings)
  //  ...but maybe, we should switch to a ZDF ladder implementation first
  // -figure out a formula for the ladder resonance gain and phase

  // hmm...modeling the resonance signal with a decaying sinusoid seems not very promising. 
  // Theoretically, the resonance signal should asymptotically look like a decaying sinusoid, but
  // as it seems, the envelope seems to deviate considerably from the exponential decay in the 
  // earlier sections of the signal. Perhaps, it's a better idea to just use the actual resonance
  // signal as is instead of trying to model it.

  // at some point, when the rsLadderResoShaped is complete, this function may be deleted.
}

void ladderResoShape()
{
  // We investigate the step responses of the lowpass and resonance output of the 
  // rsLadderResoShaped filter class.

  int    N   = 2000;        // number of samples
  double fs  = 44100;       // samplerate
  double fc  = 500;         // cutoff frequency
  double dec = 0.01;        // decay time in seconds
  double att = 0.1;         // attack-time as multiplier for decay-time
  double phs = 0.0;         // resonance phase
  double drv = 1.0;         // waveshaper drive

  // create and set up the filter:
  rsLadderResoShaped2DD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(fc);
  flt.setResonanceDecay(dec);
  flt.setDecayByFrequency(0.0);
  flt.setResonanceAttack(att);
  flt.setResonancePhase(phs);
  flt.setFeedbackSaturationGainAt1(0.75);
  //flt.setSaturationFunction(&tanh);
  flt.setDrive(drv);


  // create step-response output signals:
  vector<double> yf(N);  // pure filter signal
  vector<double> yr(N);  // pure resonance signal
  vector<double> yrp(N); // post-processed resonance signal
  vector<double> y(N);   // final output = yf + yrp
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(1.0, &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
  }

  // plot outputs:
  GNUPlotter plt;
  plt.addDataArrays(N, &y[0]);
  plt.addDataArrays(N, &yf[0]);
  plt.addDataArrays(N, &yr[0]);
  //plt.addDataArrays(N, &yrp[0]);
  plt.plot();

  // Observations: 
  // -when the phase of the modal-filter is set to 0.5*PI, the phase of the resonance signal is
  //  unchanged -> in the resShaper filter, we need to add 0.5*PI to the "phase" user parameter
  // -when the attack time is nonzero (i.e. the decay-time of the modal filter is nonzero), the
  //  resonance is also boosted - we need an amplitude scaler that depends on the decay-time of
  //  the modal filter
  //  i think, we need the sum of the geometric series S = sum_{k=0}^{\infty} r^k = 1 / (1-r)
  //  where r = e^(-dn), where dn is the normalized decay-time ( = dec*fs)
}

void ladderThresholds()
{
  int    N  = 2000;        // number of samples
  double fs = 44100;       // samplerate
  double lo = -0.4;        // low threshold
  double hi = +0.6;        // high threshold
  double drive = 8;

  // create and set up filter:
  rsLadderResoShaped2DD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1000);
  flt.setResonanceDecay(0.01);
  flt.setDecayByFrequency(0.0);
  flt.setResonanceAttack(0.0);
  //flt.setSaturationFunction(&tanh);
  //flt.setDrive(8);

  // create input signal:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, 50., fs, 0., false);

  // create output signal:
  vector<double> y(N);
  vector<double> yf(N), yr(N);   // variables for filter- and resonance signal
  vector<double> a(N);           // amplitude multiplier
  double t;                      // temporary variable
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);

    // apply saturation to resonance:
    yr[n] = tanh(drive * yr[n]);


    //t  = yf[n] + yr[n];                  // temporary signal used for thresholding
    t  = yf[n];                            // temporary signal used for thresholding
    //t  = x[n];                           // temporary signal used for thresholding
    t  = rsLinToLin(t, lo, hi, -1.0, 1.0); // map lo..hi to -1..+1
    t *= t;                                // ^2
    t *= t;                                // ^4
    t *= t;                                // ^8
    t *= t;                                // ^16
    t *= t;                                // ^32
    a[n] = 1 / (1+t);                      // amplitude scaler for resonance

    //y[n] = yf[n] + a[n] * yr[n];  // preliminary
    y[n] = yf[n] + a[n] * (yr[n]-yf[n]);  // preliminary
    //y[n] = yf[n] + a[n] * ( 0.5 *yr[n]-yf[n]);  // preliminary

    //// or try something like:
    //// if yf+yr > hi: yr -= yf+yr-hi
    //t = yf[n] + yr[n];
    //if(t > hi)
    //  yr[n] -= t-hi;
    //if(t < lo)
    //  yr[n] -= t-lo;
    //y[n] = yf[n] + yr[n];

    //// preliminary:
    //if(t < lo)
    //  y[n] = lo;
    //else if(t > hi)
    //  y[n] = hi;
    //else
    //  y[n] = t;

    //// apply thresholding to resonance signal:

    ////a[n] = exp(-t*t);    // test
    //if(t > loThresh && t < hiThresh)
    //{
    //  //a[n] = 1; // preliminary
    //  a[n] = min(hiThresh-t, t-loThresh);
    //}
    //else
    //{
    //  a[n] = 0;
    //}

    //yr[n] *= a[n];
    //y[n] = yf[n] + yr[n];  // preliminary
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &yf[0]);
  //plt.addDataArrays(N, &yr[0]);
  plt.addDataArrays(N, &y[0]);
  plt.addDataArrays(N, &a[0]);
  plt.plot();

  // ...hmmm...this seems to be a dead-end-road
  // -> try putting a saturator in the feedback path instead
}

//double feedbackSaturation(double x, double drive, double offset)
//{
//  if(drive == 0.0)
//    return x;
//  else
//    return (tanh(drive*x+offset)-tanh(offset)) / drive;
//  // in production code tanh(offset) and 1/drive can be precomputed for efficiency
//  // maybe get rid of this..
//}
void ladderFeedbackSaturation()
{
  // We investigate the effects of putting a saturator into the feedback path of a ladder filter.
  // The saturation should affect *only* the feedback signal, the whole feedforward path is free
  // of saturation. The expected effect is that tuning up the input signal will have an effect
  // of supressing the resonance in a signal dependent way.

  int    N    = 8000;        // number of samples
  double fs   = 44100;       // samplerate
  double fc   = 1500;        // cutoff frequency
  double r    = 4.0;         // normalized resonance gain (1: self-osc limit)
  double fIn  = 100;         // frequency of input sawtooth
  double aMin = 0.5;         // minimum input amplitude
  double aMax = 2.0;         // maximum output amplitude

  // create input signal - we use a sawtooth wave with amplitude increasing from aMin to aMax:
  vector<double> x(N);
  vector<double> a(N); // amplitude 
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::fillWithRangeLinear(&a[0], N, aMin, aMax);
  RAPT::rsArrayTools::multiply(&x[0], &a[0], &x[0], N);

  // compute ladder coeffs:
  double a1, b0, k, g;
  //fc = 0.007, r = 1;  // just for test
  rsLadderFilter2DD::computeCoeffs(2*PI*fc/fs, r, &a1, &b0, &k, &g);


  // compute ladder output:
  double y0, y1, y2, y3, y4;    // outputs of individual ladder stages
  y0 = y1 = y2 = y3 = y4 = 0;   // ...are initially zero
  vector<double> y(N);          // the output signal
  for(int n = 0; n < N; n++)
  {
    //y0   = x[n]  - tanh(k*y4);  // a linear filter would just have y0 = x[n] - k*y4
    y0   = tanh(x[n]) - tanh(k*y4); // 
    //y0   = tanh(x[n]) - k*tanh(y4); // 
    //y0   = tanh(x[n] - k*y4); 
    //y0   = x[n]  - k*tanh(y4);  // gain after saturation
    //y0   = x[n]  - sqrt(k)*tanh(sqrt(k)*y4); // gain equally distributed
    //y0   = x[n]  - pow(k, 0.1)*tanh(pow(k, 0.9)*y4);
    //y0   = x[n]  - 2*tanh(k*y4/2);    // less of sawtooth bleeds through
    //y0   = x[n]  - 0.5*tanh(2*k*y4);
    //y0   = x[n]  - 3*tanh(k*y4/3);

    y1   = b0*y0 - a1*y1;
    y2   = b0*y1 - a1*y2;
    y3   = b0*y2 - a1*y3;
    y4   = b0*y3 - a1*y4;
    y[n] = y4;

    //y[n] = g*y4; // maybe we should restrict the compensation gain to the range between 1 and the
                   // value that would arise from feedback-gain == 1

    //// test:
    //double drive = 0.2;
    //y[n] = tanh(drive*g*y4)/drive; 
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();

  // Observations:
  // The initial spikes at the edges look rather exaggerated compared to the analog filter. I 
  // suppose, this is due to some saturation taking place in the feedforward path of the analog 
  // filter, too. It could be the case that the feedforward path just somewhat more headroom than
  // the feedback path. ..No - it seems like the analog filter doesn't do DC-gain compensation. If
  // we leave the gain-factor out, it looks more like the analog output.

  // ToDo: incorporate an offset "o" into the saturator input, such that the saturation becomes:
  // y0   = x[n]  - (tanh(k*y4+o)-to)  // to = tanh(o)

  // maybe this can at some point be deleted when the ladderFeedbackSaturation2 is complete
}

void ladderFeedbackSaturation2()
{
  // We test the rsLadderFilterFeedbackSaturated class. As input, we use a sawtooth wave with 
  // increasing amplitude.

  // create input signal (sawtooth with amplitude increasing from aMin to aMax):
  int    N    = 10000;       // number of samples
  double fs   = 44100;       // samplerate
  double fIn  = 100;         // frequency of input sawtooth
  double aMin = 0.5;         // minimum input amplitude
  double aMax = 3.0;         // maximum output amplitude
  vector<double> x(N);
  vector<double> a(N);       // amplitude 
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  //createWaveform(&x[0], N, 1, fIn, fs, 0, true); // with anti-aliasing
  RAPT::rsArrayTools::fillWithRangeLinear(&a[0], N, aMin, aMax);
  RAPT::rsArrayTools::multiply(&x[0], &a[0], &x[0], N);

  // create the filter object and set up its parameters:
  rsLadderFeedbackSaturatedDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1500);
  //flt.setResonanceDecay(0.1);
  flt.setResonance(1.0);
  flt.setFeedbackDrive(2.0);
  flt.setLowerFeedbackLimit(-0.3);
  flt.setUpperFeedbackLimit(+0.5);
  flt.setSaturationMode(0);
  flt.setSaturationGainAt1(0.75);

  //flt.setFeedbackRange(1.0);
  //flt.setFeedbackOffset(+0.0);

  // get the filtered output signal:
  vector<double> y(N);
  for(int n = 0; n < N; n++)
  {
    //y[n] = flt.getSample(x[n]);
    y[n] = flt.getSampleNoGain(x[n]);
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();

  // Observations:

  // using a tanh as saturation function, it looks like the feedback is actually limited to
  // a range of .25 of what is passed to setFeedbackRange. it doesn't really seems to depend 
  // strongly on the feedback drive - with a range of 1000 and drive of 2.0 it's actually around
  // +-240, with drive=20 around +-270
  // maybe it's because a sinusoid at the resonance frequency that goes through the 4 filter stages
  // will have its amplitude multiplied by 1/4

  // If the edge of the saw occurs before the resonance has been diminished, the amount of 
  // excursion depends on the exact phase of the resonance at the instant of the edge - this is 
  // strange - is this some kind of sync-effect? Elan says, this effect occurs also in analog 
  // filters, so it's probably not an artifact of the digital implementation.

  // ToDo: check, how the different saturation functions respond to zero input - i think, there 
  // could be a DC be produced for certain settings (especially asymmetric settings) due to the 
  // fact that the unnormalized saturator's outtput may not exactly be zero, even for zero input
  // and that builds up in the feedback loop - maybe compensate by subtracting the output that
  // is produced for zero input...or something
}

void ladderFeedbackSaturation3()
{
  // We investigate the stage signal gain that occurs at certain settings. 
  
  // As input, we use a swatooth of fixed amplitude and we sweep the cutoff frequency 
  // exponentially.

  int    N    = 15000;       // number of samples
  double fs   = 44100;       // samplerate
  double fIn  = 50;          // frequency of input sawtooth
  double aIn  = 0.5;         // amplitude of input
  double fc1  = 1000;         // cutoff frequency at start
  double fc2  = 4000;        // cutoff frequency at end

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create cutoff frequency sweep:
  vector<double> fc(N);
  RAPT::rsArrayTools::fillWithRangeExponential(&fc[0], N, fc1, fc2);

  // create and set up filter:
  rsLadderFeedbackSaturatedDD flt;
  flt.setSampleRate(fs);
  //flt.setResonanceDecay(0.01);
  flt.setResonance(1.0);
  flt.setFeedbackDrive(2.0);
  flt.setLowerFeedbackLimit(-0.2);
  flt.setUpperFeedbackLimit(+0.2);

  // get the filtered output signal:
  vector<double> y(N);
  for(int n = 0; n < N; n++)
  {
    flt.setCutoff(fc[n]);
    y[n] = flt.getSampleNoGain(x[n]);
  }

  // plot outputs:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}

void ladderFeedbackSatDCGain()
{
  // We investigate, hwo the DC gain of the feedback-saturated ladder filter depends on the 
  // settings of the cutoff frequency, resonance-decay and resonance drive (and possibly other
  // parameters) and whether it can be reliably and accurately computed - this is important when
  // using the feedback saturated ladder inside the resohape filter.

  int    N      = 5000;       // number of input samples
  double fs     = 44100;      // samplerate
  double inGain = 1.0;        // gain of input DC signal

  // create and set up filter:
  rsLadderFeedbackSaturatedDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1000);
  //flt.setResonanceDecay(0.01);
  flt.setResonance(1.0);
  flt.setFeedbackDrive(2.0);
  flt.setLowerFeedbackLimit(-0.5);
  flt.setUpperFeedbackLimit(+0.5);
  flt.setSaturationMode(rsLadderFeedbackSaturatedDD::PRE_FB_GAIN);
  //flt.setSaturationMode(rsLadderFilterFeedbackSaturated::POST_FB_GAIN);
  //flt.setSaturationMode(rsLadderFilterFeedbackSaturated::POST_INPUT_ADD);

  // create DC input signal:
  vector<double> x(N);
  RAPT::rsArrayTools::fillWithValue(&x[0], N, inGain);

  // compute filter output signal:
  vector<double> y(N);
  int n;
  for(n = 0; n < N; n++)
    y[n] = flt.getSampleNoGain(x[n]);

  // apply a 1st order lowpass to the output to extract the DC component:
  vector<double> yDC(N);
  rsOnePoleFilterDD lpf;
  lpf.setMode(rsOnePoleFilterDD::LOWPASS_IIT);
  lpf.setSampleRate(fs);
  lpf.setCutoff(10.0);
  for(n = 0; n < N; n++)
    yDC[n] = lpf.getSample(y[n]);

  // compute the DC gain that is expected from the settings:
  double dcGain = 1.0 / flt.getCompensationGain();

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, &y[0]);
  plt.addDataArrays(N, &yDC[0]);
  plt.plot();
}

void ladderFeedbackSatReso()
{
  // We test investigate the shape of the resonance of the feedback saturated ladder filter.
  // There are these strange ramp-like artifacts around the main resonance body - we want to
  // figure out why they occur and if they can be avoided, for example by placing the saturation
  // somewhere else...

  int    N        = 2000;       // number of samples
  double fs       = 44100;      // samplerate
  double fIn      = 50;         // frequency of input sawtooth
  double aIn      = -1.0;       // amplitude of input
  double drive    = 2.0;        // feedback drive
  //double y1       = 1.0;        // feedback saturator gain at x=1 (parameter for simoid)
  //double lo       = -0.5;       // feedback saturation low-limit
  //double hi       = +0.5;       // feedback saturation high-limit

  // create and set up the filter:
  rsLadderFeedbackSaturatedDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1000);                   // cutoff frequency
  flt.setResonance(1.0);
  flt.setSaturationGainAt1(0.9);
  flt.setFeedbackDrive(2.0);
  flt.setLowerFeedbackLimit(-0.5); 
  flt.setUpperFeedbackLimit(+0.5);

  // set the position for the saturation function(s):
  //flt.setSaturationMode(rsLadderFilterFeedbackSaturated::PRE_FB_GAIN);
  flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);     // use as reference
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_1ST_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_2ND_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_3RD_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_4TH_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_EACH_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::IN_1ST_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::IN_2ND_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::IN_3RD_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::IN_4TH_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::IN_EACH_STAGE);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::TEST1);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create output signals of the individual filter stages:
  double y[5];
  vector<double> y0(N), y1(N), y2(N), y3(N), y4(N);
  int n;
  for(n = 0; n < N; n++)
  {
    flt.getSample(x[n]);
    flt.getState(y);
    y0[n] = y[0];
    y1[n] = y[1];
    y2[n] = y[2];
    y3[n] = y[3];
    y4[n] = y[4];
  }

  // plot outputs:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays(N, &y0[0]);
  //plt.addDataArrays(N, &y1[0]);
  //plt.addDataArrays(N, &y2[0]);
  //plt.addDataArrays(N, &y3[0]);
  plt.addDataArrays(N, &y4[0]);
  plt.plot();

}

void ladderFeedbackSatGrowl()
{
  // We investigate the "growling" that occurs when a sawtooth is passed through a ladder filter 
  // with high resonance.

  int    N     = 8000;       // number of samples
  double fs    = 44100;      // samplerate
  double fc    = 1025;       // cutoff frequency
  double drive = 4.0;        // feedback drive
  double fIn   = 100;        // frequency of input sawtooth
  double aIn   = 1.0;        // amplitude of input sawtooth

  // create and set up the filter:
  rsLadderResoShapedDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(fc);                             // cutoff frequency
  flt.setResonanceDecay(0.1);
  flt.setFeedbackSaturationGainAt1(1.0);         // hardclip
  flt.setFeedbackDrive(drive);
  flt.setFeedbackLowerLimit(-1.0);
  flt.setFeedbackUpperLimit(+1.0);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_GAIN_AND_ADD);


  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create filter output:
  vector<double> yf(N), yr(N), y(N);
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
    //y[n] = flt.getSample(x[n]);
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  //plt.addDataArrays(N, &yf[0]);
  plt.addDataArrays(N, &yr[0]);
  plt.plot();

  //// scale, write to file:
  //rsScale(&x[0], N, 0.3);
  //rsScale(&y[0], N, 0.3);
  //writeToStereoWaveFile("Growl_fc=1000_f=107_PostAdd.wav", &x[0], &y[0], N, (int)fs, 16);

  // Observations:
  // -fc=1kHz, POST_FB_GAIN, drive=2:
  //  -fIn=105: no growl, excursion always the same
  //  -fIn=109: strong 2-level (alternating peak excursion) growl
  //  -fIn=111: complex growl (multiple peak excursion levels)
  // -fc=1kHz, POST_INPUT_ADD, drive=1:
  //  -fIn=115: kind of squarewave amplitude modulation
  //  -we get the highest excursions, when reso-sine is at maximum at the saw-edge 
  //   -try fIn=107, look at n=3100,3500,3900

  // ToDo:
  // -what's the relation between the edge of the saw and the phase of the resonance?
  // -or is the phase being clamped when input + resonance exceeds the clipping limit?
  // -is presence or absence of growl dependent on input amplitude (if not, it may have to do with
  //  the signal edge, if, it may be related to phase-clamping due to saturation)
}

void ladderFeedbackSatGrowl2()
{
  // I try to simulate the waveshape seen in Elan's demos with the feedback saturated ladder
  // filter.

  int    N     = 8000;       // number of samples
  double fs    = 44100;      // samplerate
  double fc    = 1000;       // cutoff frequency
  double drive = 4.0;        // feedback drive
  double fIn   = 100;        // frequency of input sawtooth
  double aIn   = 1.0;        // amplitude of input sawtooth


  // create and set up the filter:
  rsLadderFeedbackSaturatedDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(fc);                 // cutoff frequency
  flt.setResonance(1.0);
  flt.setSaturationGainAt1(1.0);
  flt.setFeedbackDrive(drive);
  flt.setLowerFeedbackLimit(-1.0);   // irrelevant, if we use the experimental saturation mode
  flt.setUpperFeedbackLimit(+1.0);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  //flt.setSaturationMode(rsLadderFeedbackSaturatedDD::POST_GAIN_AND_ADD);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create filter output:
  vector<double> y(N);
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(x[n]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}

void ladderZDF()
{
  // We implement a zero-delay-feedback version of the ladder filter. Let z1,z2,.. denote the old 
  // valuees of y1,y2..., then the equations are:
  // y0 =   x  - k*y4
  // y1 = b*y0 - a*z1 =          b*(x-k*y4) - a*z1
  // y2 = b*y1 - a*z2 =       b*(b*(x-k*y4) - a*z1) - a*z2
  // y3 = b*y2 - a*z3 =    b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3
  // y4 = b*y3 - a*z4 = b*(b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3) - a*z4
  //
  // solve(y4 = b*(b*(b*(b*(x-k*y4) - a*z1) - a*z2) - a*z3) - a*z4, y4)
  // gives:
  // y4 = (-a b^3 z1-a b^2 z2-a b z3-a z4+b^4 x) / (b^4 k+1)
  // we can write this as:
  // y4 = d0*x + d1*z1 + d2*z2 + d3*z3 + d4*z4
  // where:
  // q = 1/(b^4*k + 1), d0 = q*b^4, d1 = -q*a*b^3, d2 = -q*a*b^2, d3 = -q*a*b, d4 = -q*a

  // G1(z) = b/(1+a/z), G4(z) = (G1(z))^4, H(z) = G4(z) / (1 + k*G4(z))

  // G4(z) = (b/(1+a/z))^4
  // H(z) = (b/(1+a/z))^4 / (1 + k*(b/(1+a/z))^4)

  // throw at wolfram:
  // (b/(1+a*z))^4  / (1+k (b/(1+a*z))^4)

  static const int N = 1000;

  double fs = 40000;
  double fc = 1000;    // cutoff frequency
  double fb = 0.9;    // total feedback gain

  // lowpass coefficient formula is simpler than in UDF case:
  double a, b, k;  // filter lowpass coeffs and feedback
  double wc = 2*PI*fc/fs;
  double s  = sin(wc);
  double c  = cos(wc);
  a = -1 / (s+c);
  b =  1 + a;

  // feedback gain formula is the same as in UDF case:
  double g2 = b*b / (1.0 + a*a + 2*a*c);
  k = fb / (g2*g2);

  // compute compensation gain (formula the same as in UDF):
  double g;
  g = rsLadderFilter2DD::computeCompensationGain(a, b, k);

  // compute feedback solution coefficients:
  double d[5];
  double b4 = b*b*b*b;         // b^4
  double q  = 1 / (b4*k + 1);  // scaler
  d[0] =  q*b4;
  d[4] = -q*a;
  d[3] =  b*d[4];
  d[2] =  b*d[3];
  d[1] =  b*d[2];

  // create input signal:
  vector<double> x(N);
  RAPT::rsArrayTools::fillWithZeros(&x[0], N);
  x[0] = 1;

  // compute output:
  vector<double> output(N);
  double y[5];    // filter state
  RAPT::rsArrayTools::fillWithZeros(y, 5);
  //double tmp;
  for(int n = 0; n < N; n++)
  {
    // compute/predict output of final stage:
    y[4] = d[0]*x[n] + d[1]*y[1] + d[2]*y[2] + d[3]*y[3] + d[4]*y[4];

    // update other filter states:
    y[0] =   x[n] - k*y[4];
    y[1] = b*y[0] - a*y[1];
    y[2] = b*y[1] - a*y[2];
    y[3] = b*y[2] - a*y[3];

    // write output signal:
    output[n] = g * y[4];  // more generally, a linear combination of the y[] can be used
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &output[0]);
  plt.plot();

  // Observations: 
  // It seems to work - but we need to figure out formulas for the filter coefficients and the 
  // feedback gain. The formulas will probably be simpler than in the UDF case.
}
// may be deleted when ZDF vs UDF experiemnt is complete

void ladderZDFvsUDF()
{
  // We compare the behaviors of a zero-delay-feedback (ZDF) and unit-delay-feedback (UDF) 
  // implementation of a ladder filter fo special interest is the response to audio-rate
  // modulation of the cutoff frequency.

  int    N     = 50000;          // number of samples
  double fc    = 1000;          // cutoff frequency
  double fs    = 44100;         // samplerate
  double r     = 0.99;           // normalized resonance
  double fIn   = 120;           // input wave frequency
  double aIn   = 0.125;          // input amplitude
  double fMod  = 55;            // modulation frequency
  double depth = 0.8;           // modulation depth (as scaler for fc)

  // create and set up the unit delay feedback filter:
  //rsLadderFilter2DD udf;
  rsLadderFilter<double, double> udf;
  udf.setCutoff(fc);
  udf.setSampleRate(fs);
  udf.setResonance(r);
  udf.setB1(0.23);

  // create and set up the zero delay feedback filter:
  rsLadderFilterZDFDD zdf;
  zdf.setCutoff(fc);
  zdf.setSampleRate(fs);
  zdf.setResonance(r);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., true);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create modulation signal:
  vector<double> m(N);
  createWaveform(&m[0], N, 2, fMod, fs, 0., true);

  // create filter outputs
  vector<double> yu(N), yz(N);
  double tmp;
  for(int n = 0; n < N; n++)
  {
    tmp = fc + depth*fc*m[n];     // modulated cutoff
    udf.setCutoff(tmp);
    zdf.setCutoff(tmp);
    yu[n] = udf.getSample(x[n]);
    yz[n] = zdf.getSample(x[n]);
  }

  // write output files:
  writeToMonoWaveFile("LadderModulationUDF.wav", &yu[0], N, (int)fs, 16);
  writeToMonoWaveFile("LadderModulationZDF.wav", &yz[0], N, (int)fs, 16);

  //// plot:
  GNUPlotter plt;
  ////plt.addDataArrays(N, &x[0]);
  ////plt.addDataArrays(N, &m[0]);
  plt.addDataArrays(N, &yu[0]);
  plt.addDataArrays(N, &yz[0]);
  //plt.plot();

  // Observations:
  // -Whenever the cutoff frequency switches, there is a little discontinuity in the output sinal.
  //  This gives a kind of sizzling effect. The sizzle is stronger in the ZDF version than in the 
  //  UDF version. Moreover, the ZDF version can get unstable (fc = 5000, r = 0.8, depth = 0.8)
  // -Actually, it makes qa big difference, if we use rsLadderFilter<double, double> or 
  //  rsLadderFilter2DD ...aren't they supposed to be the same? the former sounds better but is 
  //  also more quiet by a factor of about 2
  // -The bilinear UDF version sound a bit brighter and louder than the regular version. After
  //  normalizing both, they still sound different. The regular version seems to have its 
  //  "rumble" or "lawnmower" a bit more present and at a lower tone. I'd say, bilinear sounds 
  //  a bit better. Bilinear seems have a higher resonance amplitude. ToDo: write both into a 
  //  stereo wavefile (maybe with channels individually normalized) for better comparison.
  //  In some segments, the bilinear version has higher resonance (0.218) and in others, the 
  //  regular version has higher resonance (0.254). It actually makes for a nice stereo effect
  //  to have one channel regular and the other bilinear. Maybe use bilinear for mid and regular
  //  for side channel. But maybe that effect only happens with this aggressive modulation?


  // ToDo: try a cascade of 4 1st order filters with 1-pole/1-zero each, instead of just one pole
  // the additional zero can be used to enforce the gain at the cutoff frequency at each filter
  // stage to be 1/sqrt(2). This way, we do not have to compensate for a "wrong" gain via the
  // feedback factor and can use k = 4*r where r is the normalized resonance.

  // The equations are: 
  // t   = -((b1-a1*b0)*s) / (b0 + a1*b1 + (b1+a1*b0)*c)        phase of feedback at wc is -180°
  // 1/2 =  (b0^2 + b1^2 + 2*b0*b1*c) / (1 + a1^2 + 2*a1*c)     magnitude-squared at wc is 1/2
  // 1   =  (b0^2 + b1^2 + 2*b0*b1)   / (1 + a1^2 + 2*a1)       magnitude-squared at DC is 1
  // where: s = sin(wc), c = cos(wc), t = tan((wc-pi)/4)
  // the 3rd equation can be replaced by the simpler equation b0+b1-a1 = 1 -> a1 = b0+b1-1
  // i think, we also have: a1 <= 0, b1 >= 0, b0 > 0
  // we arrive at the system of 2 conic-section equations for b0, b1:
  // 0 = b0^2 - 2(1-c)b0b1 + b1^2 + 2(1-c)b0 + 2(1-c)b1 - 2(1-c)
  // 0 = (t*c-s)b0^2 + (t+t*c-s)b0b1 + t*b1^2 + (t-t*c+s)b0 + (t*c-t+s)b1
  // ...so we would need an algorithm to find the intersection of 2 conic sections to solve this 
  // problem

}

void ladderResoModulation()
{

  int dummy = 0;
}

void resoShapeFeedbackSat()
{
  // We test the combination of feedback saturation and resonance shaping with the ResoShape 
  // filter. We use a sawtooth with linearly ascending amplitude as input.

  int    N    = 8000;        // number of samples
  double fs   = 44100;       // samplerate
  double fIn  = 100;         // frequency of input sawtooth
  double aMin = 0.5;         // minimum input amplitude
  double aMax = 3.0;         // maximum output amplitude

  // create and set up the filter object:
  //rsLadderFilterFeedbackSaturated flt;
  rsLadderResoShapedDD flt;  // later: rsLadderResoShaped2
  flt.setSampleRate(fs);
  flt.setCutoff(1500);
  flt.setResonanceDecay(0.1);

  flt.setFeedbackDrive(4.0);
  flt.setFeedbackLowerLimit(-0.1);
  flt.setFeedbackUpperLimit(+0.2);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::PRE_FB_GAIN);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  flt.setFeedbackSaturationGainAt1(0.75);

  // create input signal (sawtooth with amplitude increasing from aMin to aMax):
  vector<double> x(N);
  vector<double> a(N);       // amplitude 
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::fillWithRangeLinear(&a[0], N, aMin, aMax);
  RAPT::rsArrayTools::multiply(&x[0], &a[0], &x[0], N);

  // compute output signals:
  vector<double> y(N), yr(N), yf(N);  // full, resonant and nonresonant-filtered output
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);

    // preliminary: clip the resonance signal:
    double drive  = 4.0;
    double offset = 0.0;
    //yr[n] = tanh(drive*yr[n]+offset);
    yr[n] = rsClip(drive*yr[n]+offset, -1.0, +1.0);

    y[n] = yf[n] + yr[n];
  }


  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  //plt.addDataArrays(N, &yf[0]);
  //plt.addDataArrays(N, &yr[0]);
  plt.plot();

  // things to watch out for: is the DC-gain of the resonant filter matched to the non-resonant?
  // ...if so, the pure resonance-signal (= resonant - nonresonant) should be free of DC -> this
  // can be checked.
}

void resoSaturationModes()
{
  // We investigate different approaches to appyling saturation to the resonance signal in the
  // rsResoShape filter. The main idea is to re-express the application of the stauration function
  // y = f(x) as a signal-dependent multiplication y = x * a(x) where a(x) = f(x)/x or a 
  // signal-dependent addition y = x + a(x) where a(x) = f(x)-x. This allows to use a sidechain 
  // signal as input into the respective function a(x), such that eventually, we have
  // y = x * a(z) or y = x + a(z) using the sidechain signal z, which may be the resonance-signal
  // plus some amount of the filtered output signal and/or some amount of the filtered input 
  // signal.


  int    N      = 2000;        // number of samples
  double fs     = 44100;       // samplerate
  double fc     = 880;         // cutoff frequency
  double decay  = 0.005;       // decay time in seconds
  double fIn    = 55;          // frequency of input saw-wave
  double drive  = 8.0;
  double addIn  = 0.0;         // amount by which filter-input is added to saturator-input
  double addOut = 0.5;
  double offset = 0.0;

  // create and set up the filter:
  rsLadderResoShapedDD filter;
  filter.setSampleRate(fs);
  filter.setCutoff(fc);
  filter.setResonanceDecay(decay);
  filter.setDecayByFrequency(0.0);
  //filter.setResonanceAttack(att);
  //filter.setResonancePhase(phs);

  // create input signal:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);

  // create output signals with different saturation modes:
  vector<double> yd(N), ym(N), ya(N); // direct, multiplicative, additive
  double flt, res;                    // filter-output and resonance signals
  double z;                           // sidechain signal 
  double fz;                          // f(z) - saturated sidechain signal
  double t, ft;
  double a;
  //double tmp;
  for(int n = 0; n < N; n++)
  {
    // get filtered output and resonance signal:
    filter.getSignalParts(x[n], &flt, &res);

    // create the sidechain signal z:
    t = addIn*x[n] + addOut*flt + offset;
    z = drive*res + t;

    // create the saturated sidechain-signal f(z):
    fz = tanh(z);
    ft = tanh(t);   // would be nice, if we don't need that

    // direct application of saturation:
    //yd[n] = flt + fz;
    yd[n] = fz;  // preliminary - we just want to see the saturated resonance

    // multiplicative:
    if(z == 0.0)
      a = 1.0;
    else
      a = fz/z; 
    ym[n] = flt + drive*res * a;

    // additive:
    a = fz-z;
    ya[n] = flt + (drive*res + a);
    //ya[n] = drive*res + a;
    ya[n] = drive*res + a - 0.5*t;
    //ya[n] = flt + a;
    //ya[n] = x[n] + (drive*res + a);
    //ya[n] = flt + (drive * res + a) - ft;
  }

  // plot outputs:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &yd[0]);
  //plt.addDataArrays(N, &ym[0]);
  plt.addDataArrays(N, &ya[0]);
  plt.plot();

  // Observations:
  // The multiplicative seems useless but the additive might be promising - i need to figure out
  // how to get rid of the bump at the edge
}

void resoShapeGate()
{
  // We test the gating of the resonance signal.
  // hmmm - the gate seems a dead end road - maybe delete it.

  int    N        = 4000;       // number of samples
  double fs       = 44100;      // samplerate
  double fIn      = 50;         // frequency of input sawtooth
  double aIn      = -1.0;       // amplitude of input
  double drive    = 2.0;        // feedback drive
  double y1       = 1.0;        // feedback saturator gain at x=1 (parameter for simoid)
  double lo       = -0.5;       // feedback saturation low-limit
  double hi       = +0.5;       // feedback saturation high-limit

  // create and set up the filter:
  rsLadderResoShaped2DD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1000);                   // cutoff frequency
  flt.setResonanceDecay(1.0);
  flt.setFeedbackSaturationGainAt1(y1);
  flt.setFeedbackDrive(drive);
  //flt.setFeedbackDrive(drive / (y1*y1*y1*y1)); // test - compensated drive
  flt.setFeedbackLowerLimit(lo); 
  flt.setFeedbackUpperLimit(hi);
  flt.setGateSensitivity(1.0);
  flt.setGateMix(0.0);

  // set the position for the saturation function(s):
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::PRE_FB_GAIN);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_1ST_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_2ND_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_3RD_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_4TH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_EACH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_1ST_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_2ND_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_3RD_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_4TH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_EACH_STAGE);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create output signals (filtered, resonance, complete)
  vector<double> y(N), yr(N), yf(N);
  int n;
  for(n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
  }

  // plot outputs:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  //plt.addDataArrays(N, &yf[0]);
  //plt.addDataArrays(N, &yr[0]);
  plt.plot();

  //// normalize filter output and write to wavefile:
  //rsNormalize(&y[0], N, 1.0, true);
  //writeToMonoWaveFile("ResoShapeTest.wav", &y[0], N, (int)fs, 16);

  // Observations:
  // In addition to the expected restricted resonance in the middle section of the saw, there 
  // are some strange ramp-like signal components in the resonance signal (in POST_FB_GAIN mode)
  // These do not get samller when increasing the sample rate, so they are unlikely to be
  // artifacts of the unit delay in the feedback path.
  // ToDo:
  // -we could try to punch out the middle section of the signal by using a gate, sidechained on
  //  the input or lowpass signal
  // -more elegant would be to get at the bottom why they occur in the first place and maybe try 
  //  to avoid them ...maybe a highpass in the feedback path could help?
  // -try some more saturation places - maybe there is one setting where it works as desired, maybe
  //  we can also try applying saturation in multiple places (between or inside all stages)

  // -PRE_FB_GAIN and POST_FB_GAIN look similar when in the PRE_FB_GAIN case, the lower and upper
  //  limits are scaled by 1/4.
  // -feedback drive should be scaled by 1/y1^4
  // -in POST_INPUT_ADD mode, the full filter output signal is actually a much better candidate
  //  for the "pure resonance", than resonant-minus-nonresoant output
  // -it seems like, as soon as we enter the domain of self-oscillation (implying a saturating
  //  feedback loop), the model resonance = resonant - nonresonant is not really good anymore

  // ToDo:
  // build gate into the filter - two parameters: 
}

void resoShapePseudoSync()
{
  // Create a pseudo oscillator sync sound using the ResoShape filter. The self-oscillating filter
  // takes the role of the sync-slave, the sawtooth input is the sync master. We use a fixed 
  // sawtooth (master) frequency and let the filter (slave) frequency do an exponential sweep.

  // user parameters:
  int    N        = 100000;      // number of samples
  double fs       = 44100;      // samplerate
  double fIn      = 50;         // frequency of input sawtooth
  double fc1      = 4000;       // cutoff frequency at start
  double fc2      = 700;        // cutoff frequency at end
  double aIn      = -1.0;       // amplitude of input
  double drive    = 2.0;        // feedback drive
  double y1       = 1.0;        // feedback saturator gain at x=1 (parameter for simoid)
  double lo       = -0.5;       // feedback saturation low-limit
  double hi       = +0.5;       // feedback saturation high-limit

  // create and set up the filter:
  rsLadderResoShaped2DD flt;
  flt.setSampleRate(fs);
  flt.setResonanceDecay(0.005);
  flt.setFeedbackSaturationGainAt1(y1);
  //flt.setFeedbackDrive(drive);
  flt.setFeedbackDrive(drive / (y1*y1*y1*y1)); // test - compensated drive
  flt.setFeedbackLowerLimit(lo); 
  flt.setFeedbackUpperLimit(hi);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_1ST_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_2ND_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_3RD_STAGE);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_4TH_STAGE);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create cutoff sweep:
  vector<double> fc(N);
  RAPT::rsArrayTools::fillWithRangeExponential(&fc[0], N, fc1, fc2);

  // get filter output:
  vector<double> y(N);
  for(int n = 0; n < N; n++)
  {
    flt.setCutoff(fc[n]);
    y[n] = flt.getSample(x[n]);
  }

  // normalize filter output and write to wavefile:
  RAPT::rsArrayTools::normalize(&y[0], N, 1.0, true);
  writeToMonoWaveFile("ResoShapePseudoSync.wav", &y[0], N, (int)fs, 16);
}

void resoSeparationNonlinear()
{
  // We experiment with separating the resonance signal is the case of a self-oscillating nonlinear
  // setting. The goal is to figure out, whether and how a meaningful resonance separation is 
  // possible. This may depend on the point at which the feedback saturation is applied.

  // maybe we should use the ResoReplace filter here already...


  int    N        =  4000;      // number of samples
  double fs       =  44100;     // samplerate
  double fc       =  1000;      // cutoff frequency
  double fIn      =  100;       // frequency of input sawtooth
  double aIn      =  1.0;       // amplitude of input
  double drive    =  2.0;       // feedback drive
  double y1       =  1.0;       // feedback saturator gain at x=1 (parameter for sigmoid)
  double lo       = -1.0;       // feedback saturation low-limit
  double hi       = +1.0;       // feedback saturation high-limit

  // create and set up the filter:
  rsLadderResoShaped2DD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(fc);
  flt.setResonanceDecay(1.0);
  flt.setFeedbackSaturationGainAt1(y1);
  flt.setFeedbackDrive(drive);
  flt.setFeedbackLowerLimit(lo); 
  flt.setFeedbackUpperLimit(hi);

  // set the position for the saturation function(s):
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::PRE_FB_GAIN);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_1ST_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_2ND_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_3RD_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_4TH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_EACH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_1ST_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_2ND_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_3RD_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_4TH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_EACH_STAGE);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create output signals (filtered, resonance, complete)
  vector<double> y(N), yr(N), yf(N);
  int n;
  for(n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays(N, &yf[0]);
  plt.addDataArrays(N, &y[0]);
  plt.addDataArrays(N, &yr[0]);
  plt.plot();


  int dummy = 0;
}





void resoReplace()
{
  // Tests the replacement of the resonance's sine waveform by an arbitrary waveform implemented
  // in rsLadderResoShaped3

  // user parameters:
  int    N     = 3000;       // number of samples
  double fs    = 44100;      // samplerate
  double fIn   = 50;         // frequency of input sawtooth

  // create and set up the filter:
  rsResoReplacerDD flt;
  flt.setSampleRate(fs);
  flt.setCutoff(1000.0);
  flt.setResonanceDecay(0.01);
  flt.setResonanceAttack(0.0);
  flt.setResonancePhase(PI); 
  flt.setResonanceGain(1.0);
  flt.setResonanceWaveform(4);
  flt.setResoCutoffMultiplier(0.5);
  flt.setResoCutoffModulation(20.0);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);

  // get filter output:
  vector<double> y(N), yr(N), yf(N);
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
  }

  //// normalize filter output and write to wavefile:
  //rsScale(&y[0], N, 0.3);
  //writeToMonoWaveFile("ResoReplaceSaw2.wav", &y[0], N, (int)fs, 16);

  // plot outputs:
  GNUPlotter plt;
  plt.setPixelSize(1000, 300);
  //plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  //plt.addDataArrays(N, &yf[0]);
  //plt.addDataArrays(N, &yr[0]);
  plt.plot();

  // todo: we may try to model the resonance signal as an enveloped sine plus a transient and
  // replace only the enveloped sine. this should get rid or reduce the artifacts at the edge

  // apply gating/enveloping based on the input signal to simulate the feedback saturation
  // effect...
}

void resoReplacePhaseBumping()
{
  // We bump the phase of the reso-replace-filter using a signal derived from
  // hp   = highpass(input);      // turn edges into spikes
  // tmp  = fabs(hp) * (hp-res);  // maintain spikes, when spike and reso are antiphase
  //   or: tmp  = fabs(hp) * (in-res);
  // 
  // tmp  = lowpass1(tmp);        // spike -> exponential decay
  // tmp  = lowpass2(tmp);        // decay -> attack/decay
  // phs += tmp;                  // bump phase


  int    N   = 10000;        // number of samples
  double fs  = 44100;       // samplerate
  double fc  = 1037;        // cutoff
  double fIn = 100;         // input saw frequency
  double aIn = 1.0;         // input saw amplitude


  // create and set up filter:
  rsResoReplacerPhaseBumpedDD flt;
  flt.setCutoff(fc);
  flt.setResonanceDecay(0.02);
  flt.setBumpFactor(5.0);
  flt.setBumpCutoffs(20, 100);
  flt.setChaosParameter(1.0);

  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // get filter output:
  vector<double> y(N), yr(N), yf(N);
  vector<double> b(N);  // phase bump signal
  for(int n = 0; n < N; n++)
  {
    flt.getSignalParts(x[n], &yf[n], &yr[n]);
    y[n] = yf[n] + yr[n];
    b[n] = flt.getPhaseBump();
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.setPixelSize(1000, 300);
  plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays(N, &y[0]);
  //plt.addDataArrays(N, &yf[0]);
  plt.addDataArrays(N, &yr[0]);
  plt.addDataArrays(N, &b[0]);
  plt.plot();

  // Observations:
  // -fc=1kHz, decay=0.02
  //  -fIn=109Hz: saw edge and resonance maximum coincide -> strong resonance
  //  -fIn=102Hz: saw edge and resonance minimum coincide -> weak resonance
}

void resoReplaceScream()
{
  // We try to recreate the scream sound from Elan's youtube demonstrations of an analog filter
  // with the resoreplace filter. We use a sawtooth input and produce a filter-sweep.

  // user parameters:
  double length = 5.0;         // length of the output in seconds
  double fc1    = 2000;        // cutoff at start
  double fc2    = 1000;        // cutoff at end
  double fIn    = 50;          // input sawtooth frequency
  double aIn    = 0.5;         // input sawtooth amplitude
  double fs     = 44100;       // sample rate


  // create and set up the filter:
  // create and set up the filter:
  rsResoReplacerPhaseBumpedDD flt;
  flt.setSampleRate(fs);
  flt.setResonanceDecay(1.0);
  flt.setFeedbackDrive(4.0);            // seems to control the growl
  flt.setFeedbackLowerLimit(-0.2);
  flt.setFeedbackUpperLimit(+0.2);
  flt.setFeedbackSaturationGainAt1(1.0);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_FB_GAIN);
  flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_INPUT_ADD);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::POST_EACH_STAGE);
  //flt.setFeedbackSaturationPlace(rsLadderFeedbackSaturatedDD::IN_EACH_STAGE);
  //flt.setResonanceAttack(0.0);
  //flt.setResonancePhase(PI); 
  //flt.setResonanceGain(1.0);
  flt.setResonanceWaveform(2);        // square 
  flt.setResoCutoffMultiplier(2.0);
  flt.setResoCutoffModulation(0.0);  // maybe this cutoff modulation is useless?
  flt.setInputRange(0.2);             // rename this function (name is nondescriptive)
  flt.setAmplitudeLimit(1.0);



  int N = (int) ceil(fs*length);  // length in samples
  
  // create input sawtooth:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fIn, fs, 0., false);
  RAPT::rsArrayTools::scale(&x[0], N, aIn);

  // create frequency sweep:
  vector<double> f(N);
  RAPT::rsArrayTools::fillWithRangeExponential(&f[0], N, fc1, fc2);

  // compute output:
  vector<double> y(N);
  for(int n = 0; n < N; n++)
  {
    flt.setCutoff(f[n]);
    y[n] = flt.getSample(x[n]);
  }


  //// plot outputs:
  //GNUPlotter plt;
  //plt.setPixelSize(1000, 300);
  ////plt.addDataArrays(N, &x[0]);
  //plt.addDataArrays(N, &y[0]);
  //plt.plot();

  // write wavefile
  writeToMonoWaveFile("ResoReplaceScream.wav", &y[0], N, (int)fs, 16);
}

void fakeResonance()
{
  // We test the fake-resonance filter by feeding a sawtooth into it and plotting the result:

  int    N     = 2000;        // number of samples
  double fs    = 44100;       // samplerate
  double fi    = 50;          // input signal frequency
  double ai    = 1.0;         // input amplitude
  double fc    = 500;         // lowpass cutoff
  double fh    = 20000;       // highpass cutoff (maybe use a differencer as highpass later)
  double shift = 0.0;         // resonance pitch shift in semitones
  double dec   = 0.003;       // resonance decay
  double att   = 0.01;        // resonance attack, scalefactor k apllied to dec, 0 < k < 1
  double ar    = 1.0;         // resonance amplitude
  double pr    = 0.0*PI;      // resonance start phase in radians
  double dl    = 1.0;         // resonance delay (as factor for optimal dealy)

  // create and set up the filter object:
  rsFakeResonanceFilterDD flt;
  flt.setSampleRate(fs);
  flt.setLowpassCutoff(fc);
  flt.setHighpassCutoff(fh);
  flt.setResonanceShift(shift);
  flt.setResonanceAmplitude(ar);
  flt.setResonanceAttack(att);
  flt.setResonanceDecay(dec);
  flt.setResonancePhase(pr);
  flt.setResonanceDelay(dl);
  
  // create sawtooth wave as test input:
  vector<double> x(N);
  createWaveform(&x[0], N, 1, fi, fs, 0.0, false);
  //createWaveform(&x[0], N, 1, fi, fs, 0.0, true);

  // create lowpass, resonance and summed output:
  vector<double> yl(N), yr(N), y(N);
  for(int n = 0; n < N; n++)
  {
    yl[n] = flt.getLowpassOutput(x[n]);
    yr[n] = flt.getResonanceOutput(x[n]);
    y[n]  = yl[n] + yr[n];
  }

  // plot input and output:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &yl[0]);
  plt.addDataArrays(N, &yr[0]);
  plt.addDataArrays(N, &y[0]);

  // Observations:
  // The resonance phase parameter doesn't seem to have any effect - why?

  // todo: make a test how the resonance signal looks like with different frcational delay
  // settings 0.0, 0.25, 0.5, 0.75, 1.0 (for the frcational part) - compare linear with allpass
  // interpolation. Take special care to see, if at high resonance frequency there's any kind
  // of damping with linear interpolation when f=0.5 - if so, check if allpass interpolation
  // solves this (it should)
  // maybe look at the resonance impulse responses

  plt.plot();
}

void fakeResoLowpassResponse()
{
  // We plot the impulse and step response of the lowpass path of our fake resonance filter

  int    N  = 1000;        // number of samples
  double fs = 44100;       // samplerate
  double fc = 500;         // lowpass cutoff

  rsFakeResonanceFilterDD flt;
  flt.setSampleRate(fs);
  flt.setLowpassCutoff(fc);

  // get impulse- and step response:
  vector<double> h(N), s(N);
  h[0] = flt.getLowpassOutput(1.0);
  for(int n = 1; n < N; n++)
    h[n] = flt.getLowpassOutput(0.0);
  RAPT::rsArrayTools::normalize(&h[0], N, 1.0);

  flt.reset();
  for(int n = 0; n < N; n++)
    s[n] = flt.getLowpassOutput(1.0);

  GNUPlotter plt;
  plt.addDataArrays(N, &h[0]);
  plt.addDataArrays(N, &s[0]);
  plt.plot();
}

void fakeResoDifferentDelays()
{
  // We test the fake-resonance filter with different settings for the delay, such that we can 
  // check for changes of the resonance depending on the fractional part of the delaytime.

  int    N     = 2000;        // number of samples
  double fs    = 40000;       // samplerate

  // create and set up the filter object:
  rsFakeResonanceFilterDD flt;
  flt.setSampleRate(fs);
  flt.setLowpassCutoff(5000);
  flt.setResonanceAmplitude(1.0);
  flt.setResonanceDecay(0.003);
  flt.setResonanceAttack(0.01);

  // create output arrays:
  vector<double> y1(N), y2(N), y3(N), y4(N), y5(N);
  int k = 400; // basic offset/delay
  flt.setResonanceDelay(k/fs);
  getImpulseResponse(flt, &y1[0], N);
  flt.setResonanceDelay((k+0.25)/fs);
  getImpulseResponse(flt, &y2[0], N);
  flt.setResonanceDelay((k+0.5)/fs);
  getImpulseResponse(flt, &y3[0], N);
  flt.setResonanceDelay((k+0.75)/fs);
  getImpulseResponse(flt, &y4[0], N);
  flt.setResonanceDelay((k+1)/fs);
  getImpulseResponse(flt, &y5[0], N);

  // extract the envelopes:
  rsOnePoleFilterDD smoother;
  smoother.setCutoff(50.0);
  smoother.setSampleRate(fs);
  smoother.setMode(rsOnePoleFilterDD::LOWPASS_IIT);
  int n;
  smoother.reset();
  for(n = 0; n < N; n++)
    y1[n] = smoother.getSample(y1[n]*y1[n]);
  smoother.reset();
  for(n = 0; n < N; n++)
    y2[n] = smoother.getSample(y2[n]*y2[n]);
  smoother.reset();
  for(n = 0; n < N; n++)
    y3[n] = smoother.getSample(y3[n]*y3[n]);
  smoother.reset();
  for(n = 0; n < N; n++)
    y4[n] = smoother.getSample(y4[n]*y4[n]);
  smoother.reset();
  for(n = 0; n < N; n++)
    y5[n] = smoother.getSample(y5[n]*y5[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, &y1[0]);
  plt.addDataArrays(N, &y2[0]);
  plt.addDataArrays(N, &y3[0]);
  plt.addDataArrays(N, &y4[0]);
  plt.addDataArrays(N, &y5[0]);
  plt.plot();

  // With linear interpolation, the signal for the fractional offset of 0.5 is indeed much quieter
  // than with 0.0. Warped allpass interpolation solves this and should indeed be used here.
}
