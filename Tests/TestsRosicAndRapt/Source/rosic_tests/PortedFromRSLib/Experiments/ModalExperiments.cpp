#include "ModalExperiments.h"

void twoPoleFilter()
{
  // rsTwoPoleFilter is a sort of simplified version of the modal filters...

  double fs = 44100;
  double f  = 1000;
  double bw = 10;

  rsTwoPoleFilter<double, double> flt;
  double s = 2*PI/fs;  // conversion factor from frequency in Hz to omega
  double w = s*f;
  flt.setFrequencyAndAbsoluteBandwidth(w, s*bw);

  double g  = flt.getMagnitudeAt(w);
  double dB = rsAmp2dB(g);

  plotImpulseResponse(flt, 500, 1.0);
  plotFrequencyResponse(flt, 1001, 50.0, 20000.0, fs, true);

  // Observations:
  // -the maximum amplitude of the decaying sine seems to be (roughly) independent from the 
  //  bandwidth
}

void modalFilter()
{
  //static const int N = 5000;  // number of samples to plot

  double fs  = 44100;  // samplerate in Hz
  double frq = 100;    // frequency in Hz
  double amp = 1.0;    // amplitude as raw factor
  double dec = 0.1;    // decay time constant in seconds
  double phs = 0;      // phase in degrees

  // create and set up the modal filter:
  rsModalFilterDD mf;
  mf.setModalParameters(frq, amp, dec, phs, fs);

  plotImpulseResponse(mf, 10000, 1.0);
  plotFrequencyResponse(mf, 2000, 20.0, fs/2, fs, true);
}

void modalFilterFreqResp()
{
  // filter parameters:
  double fs  = 44100;  // samplerate in Hz
  double frq = 100;    // frequency in Hz
  double amp = 1.0;    // amplitude as raw factor
  double att = 0.05;   // attack time in seconds
  double dec = 0.1;    // decay time constant in seconds
  double phs = 45;     // phase in degrees


  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(frq, amp, att, dec, phs, fs, 1.0);

  plotImpulseResponse(  mf, 10000, 1.0);
  plotFrequencyResponse(mf, 5000, 20.0, fs/2, fs, true);
}

// hmm...this is now a bit redundant:
void attackDecayFilter()
{
  static const int N = 20000;  // number of samples to plot

  double fs  = 44100;  // samplerate in Hz
  double ta  = 0.05;   // attack time in seconds
  double td  = 0.2;    // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = 45;     // phase in degrees
  double a   = 1.0;    // amplitude as raw factor

  //ta = 0.00001; // for comparison - todo: allow exact zero attack value

  // create and set up the filter:
  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(f, a, ta, td, phs, fs);

  // generate and plot impulse-response:
  double h[N];
  getImpulseResponse(mf, h, N);
  plotData(N, 0, 1/fs, h);
  int dummy = 0;
}


void modalTwoModes()
{
  double fs  = 44100;
  double at  = 0.05;  // attack - needs to be nonzero for numerical reasons (->try to fix)

  // mode parameter arrays:
  typedef std::vector<double> Vec;
  Vec frq = { 1.0, 4.0 };
  Vec amp = { 1.0, 1.0 };
  Vec phs = { 0.0, 0.0 };
  Vec dec = { 0.2, 0.2 };
  Vec att = { at,  at  };

  // create and set up object:
  rosic::rsModalFilterBankDD mfb;
  mfb.setSampleRate(fs);
  mfb.setReferenceFrequency(100);
  mfb.setReferenceDecay(1.0);
  mfb.setReferenceAttack(1.0);
  mfb.setModalParameters(frq, amp, att, dec, phs);

  // plots:
  plotImpulseResponse(  mfb, 5000, 1.0);
  plotFrequencyResponse(mfb, 5000, 0.0, 500.0, fs, false);

  // Observations:
  // frq: 1,4; amp: 1,1; dec: 0.2,0.2, att: 0.0001:
  // phs: 0,0: 
  //  -two sharp resonances at 100 and 400 Hz, 1 sharp antiresonance at 200 Hz 
  //   (geometric mean of 100 and 400)
  //  -phase response goes from 0° down to -180° and back up to 0°, passing through 90° at the
  //   resonances and antiresonance, transition is quite steep/squarish
  // phs: 0,pi or pi,0
  //  -two sharp resonances at 100,400
  //  -phase monotonically decresing from 0 to -90 then to -180
  //   ...wait - it starts at +180 - something must be wrong about the plot/unwrapping

  // interesting: when increasing the attack time, the freq-resp first deviates away from how it 
  // initially (at=0.0001) looks but seems to come back to a similar look at at=0.1 with some weird
  // wiggling in between
}




/** Like rsDampedSineFilter, but with the global gain factor factored out in the transfer function 
numerator, such that:

         1 + b1*z^-1
H(z) = g ---------------------
         1 + a1*z^-1 + a2*z^-2

which can be implemeneted as the difference equation:

w[n] = x[n] + b1*x[n-1] - a1*w[n-1] - a2*w[n-2], y[n] = g*w[n] */
void rsDampedSineFilterNormalizedB0(double w, double A, double d, double p, double *g, double *b1, 
  double *a1, double *a2)
{
  // calculate intermediate variables:
  double P, pp, ri, R;
  P  = exp(-1.0/d);
  p  = rsWrapToInterval(p, 0, 2*PI);
  pp = p-PI/2;
  ri = 0.5*tan(pp);
  R  = rsSqrt(0.25+ri*ri);

  // todo: use rsSinCos, double-angle formula:
  double c1, c2, s1, s2;
  c1 = cos(w);
  c2 = cos(2*w);
  s1 = sin(w);
  s2 = sin(2*w);

  // calculate coefficients:
  *a2 = P*P;
  *a1 = -2*P*c1;
  *b1 = -(P/2)*(2*(1-c2)*ri+s2)/s1; 
  *g  = A/(2*R);
  if( p > PI )
    *g = -*g;
}
// this implementation is obsolete - the new one has simpler formulas


/** Retrieves damped sine filter design parameters from its coefficients. See 
@see rsDampedSineFilter for meaning of parameters. The phase p is returned in the interval 
0...2pi. */
/*
// moved to Prototypes:
void rsDampedSineFilterAnalysis(double b0, double b1, double a1, double a2, double *w, double *A, 
  double *d, double *p)
{
  rsAssert(0.25*a1*a1-a2 < 0.0, "no damped sine filter, poles not complex conjugate");
  double P, cw;
  P  = sqrt(a2);
  cw = -0.5*a1/P;
  *d = -1.0/log(P);
  *w = acos(cw);
  *p = atan2(sin(*w),b1/(P*b0)+cw);
  if( rsAbs(b0) > rsAbs(b1) )
    *A = b0/sin(*p);
  else
    *A = b1/(P*sin(*w-*p));
  if( *A < 0.0 )
  {
    *A  = -*A;
    *p += PI;
  }
}

void rsDampedSineFilterAnalysis2(double b0, double b1, double a1, double a2, double *w, double *A, 
  double *d, double *p)
{
  rsAssert(0.25*a1*a1-a2 < 0.0, "no damped sine filter, poles not complex conjugate");
  rsComplexDbl j(0.0, 1.0);                // imaginary unit
  double P = sqrt(a2);                     // pole radius
  *w = acos(-0.5*a1/P);                    // pole angle
  rsComplexDbl q = P * exp(j * *w);        // pole location
  rsComplexDbl r = (b1+b0*q)/(2*q.imag()); // residue location
  *d = -1.0/log(P);                        // normalized decay time constant
  *A = 2*abs(r);                           // amplitude
  *p = arg(r);                             // start phase...
  if( *p < 0.0 )                           // ...in interval 0...2pi instead of -pi...pi
    *p += 2*PI;
  // Remark: There are actually two mathematical errors in this sequence of assignments which 
  // conveniently cancel each other and streamline the implementation, that's why I left them in.
  // Actually, it should be r = (b1+b0*q)/(2.0*j*q.im) and *p = arg(r) + 0.5*PI, so we have 
  // first missed a division by j (corresponding to a rotation by -pi/2) in the computation of the
  // residue r and that's why we later don't need to add pi/2 to the startphase value ;-)
}
*/

void rsDampedSineFilterOld(double w, double A, double d, double p, double *b0, double *b1, 
  double *a1, double *a2)
{
  double g;
  rsDampedSineFilterNormalizedB0(w, A, d, p, &g, b1, a1, a2);
  *b0  = g;
  *b1 *= g;
}

void dampedSineFilterDesign()
{
  static const int N = 5000;  // number of samples to plot

  // user parameters:
  double fs  = 44100;  // samplerate in Hz
  double td  = 0.02;   // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  //double p   = PI/3;   // start-phase in radians
  double p   = 5.0;   // start-phase in radians
  //double p   = 2*PI+0.001;   // start-phase in radians
  double A   = 1.0;    // amplitude as raw factor

  // compute normalized variables:
  double w   = 2*PI*f/fs;  // normalized radian frequency
  double d   = td*fs;      // normalized decay time

  // compute coefficients:
  double b0, b1, a1, a2;
  double B0, B1, A1, A2;
  rsDampedSineFilterOld(   w, A, d, p, &B0, &B1, &A1, &A2);
  rsDampedSineFilterCoeffs(w, A, d, p, &b0, &b1, &a1, &a2);


  // recover design parameters:
  double ww, AA, dd, pp;
  rsDampedSineFilterAnalysis( b0, b1, a1, a2, &ww, &AA, &dd, &pp);
  rsDampedSineFilterAnalysis2(b0, b1, a1, a2, &ww, &AA, &dd, &pp);


  int dummy = 0;
}

void dampedSineFilterImpResp()
{
  // We plot the desired analytic impulse response of a damped sine filter as a pseudo-continuous
  // graph and on top of it the samples of the actual damped sine filter. ...and then also an
  // oversampled version - this was to figure out, why the oversampled signal seemed time-squished
  // with respect to the analytic and non-oversampled signal - it turned out, that the time axis
  // was wrong - stupid mistake!

  // decaying sine parameters:
  double amplitude = 1.0;   // overall amplitude
  double phase     = 45;    // start-phase in degrees
  double decay     = 50;    // number of samples to decay to A/e
  double freq      = 0.1;   // normalized frequency (freq/sampleRate)

  int N = 100;

  // create pseudo-continuous impulse response:
  typedef RAPT::rsArrayTools AR;
  int Nc = N * 20;  // 20 times the sample-rate
  std::vector<double> tc(Nc), yc(Nc);
  AR::fillWithRangeLinear(&tc[0], Nc, 0.0, N-1.0);
  for(int n = 0; n < Nc; n++)
    yc[n] = amplitude * exp(-tc[n]/decay) * sin(2*PI*freq*tc[n] + rsDegreeToRadiant(phase));

  // create sampled impulse-response:
  double a[3], b[2];
  a[0] = 1.0;
  rsDampedSineFilterCoeffs(2*PI*freq, amplitude, decay, rsDegreeToRadiant(phase),
    &b[0], &b[1], &a[1], &a[2]);
  std::vector<double> y(N);
  AR::fillWithZeros(&y[0], N); y[0] = 1;
  AR::filter(&y[0], N, &y[0], N, b, 1, a, 2);

  // create oversampled impulse-response:
  int oversampling = 4;
  int No = N * oversampling;
  rsDampedSineFilterCoeffs(2*PI*freq/oversampling, amplitude, decay*oversampling,
    rsDegreeToRadiant(phase), &b[0], &b[1], &a[1], &a[2]);
  std::vector<double> to(No), yo(No);
  //AR::fillWithRangeLinear(&to[0], No, 0.0, N-1.0);               // wrong!
  AR::fillWithRangeLinear(&to[0], No, 0.0, (No-1.0)/oversampling); // correct!
  AR::fillWithZeros(&yo[0], No); yo[0] = 1;
  AR::filter(&yo[0], No, &yo[0], No, b, 1, a, 2);


  GNUPlotter plt;

  // add pseudo-continuous data (from analytic expression)
  plt.addDataArrays(Nc, &tc[0], &yc[0]);  
  plt.addGraph("index 0 using 1:2 with lines lw 2 lc rgb \"#808080\" notitle");

  // add sampled data:
  plt.addDataArrays(N, &y[0]);          
  plt.addGraph("index 1 with points pt 7 ps 0.8 lc rgb \"#000080\" notitle");

  // add oversampled data:
  plt.addDataArrays(No, &to[0], &yo[0]);
  plt.addGraph("index 2 with points pt 7 ps 0.8 lc rgb \"#008000\" notitle");

  plt.plot();
}



/** A genaral time-domain design of biquad filters can be thought of as a damped-sine filter 
(two-pole/one-zero) plus a scaled delta-impulse at the beginning of the impulse-response
...check, if this is actually true... */
void rsDelayedDampedSinePlusScaledDeltaFilter(double w, double A, double d, double p, double c,
  double *b0, double *b1, double *b2, double *a1, double *a2)
{
  double g, c1;
  rsDampedSineFilterNormalizedB0(w, A, d, p, &g, &c1, a1, a2);
  *b0 = c;
  *b1 = c * (*a1) + g;
  *b2 = c * (*a2) + g*c1;
}

void biquadImpulseResponseDesign()
{
  static const int N = 5000;  // number of samples to plot

  // design damped-sinusoid biquad filter:
  double fs  = 44100;  // samplerate in Hz
  double td  = 0.02;   // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = PI/4;   // start-phase in radians
  double A   = 1.0;    // amplitude as raw factor
  double c   = 0.2;    // scaler for the unit impulse

  // design the filter:
  double g, a1, a2, b0, b1, b2, c1;
  rsDampedSineFilterNormalizedB0(2*PI*f/fs, A, td*fs, phs, &g, &c1, &a1, &a2);
  b0 = c;
  b1 = b0*a1+g;
  b2 = b0*a2+g*c1;


  //rsDelayedDampedSinePlusScaledDeltaFilter(

  // copy filter coefficients into arrays, suitable for rsFilter:
  double a[3], b[3];
  a[0] = 0;  a[1] = a1; a[2] = a2;
  b[0] = b0; b[1] = b1; b[2] = b2;

  // generate time-axis and impulse-response:
  double t[N], x[N];
  createTimeAxis(N, t, fs);
  RAPT::rsArrayTools::fillWithZeros(x, N);
  x[0] = 1;
  RAPT::rsArrayTools::filter(x, N, x, N, b, 2, a, 2);

  // plot the impulse-response (versus the time-axis):
  plotData(N, t, x);
}

void modalBankTransient()
{
  double fs  = 44100;  // samplerate in Hz
  //double fs  = 2000;      // samplerate in Hz
  double f   = 100;       // fundamental frequency
  double length = 4.0;    // total length in seconds

  double attack = 0.1;   // master attack
  double decay  = 0.5;   // master decay time in seconds
  double nonlin = +0.12;  // nonlinear feedback strength

  // transient parameters:
  double tr = 0.5 * (1+sqrt(5.0)); // golden ratio - relative frequency of "transient mode"
  double ta = 1.0;                 // transient strength/amplitude
  double tl = 0.2;                 // transient length
  double tp = 0;                   // transient phase - try pi too -> phase rresponse more monotonic

  // use 3 modes at relative frequencies 1,2,3 where 2 decays faster than 3 (generally, we may have
  // and even/odd decay balance.

  // use an additional fast-decaying mode at an irrational relative frequency, like the golden 
  // ratio (=1.618... which is the number which is hardest to approximate by a rational number)
  // ..then use nonlinear feedback to produce chaos in the transient

  // try "lowpass" attack/decay filters to model a transient - or do my regular modal filters
  // admit zero frequency? ...maybe if not, i could just make a dispatcher to switch to another 
  // design formula, if w==0

  typedef std::vector<double> Vec;


  // relative mode frequencies, amplitudes, phases, decay-times, attack times
  Vec frq = { 1,    2,    3, tr }; // use last mode as transient (preliminary - later have a
  Vec amp = { 1, 1./2, 1./3, ta }; // dedicated interface for setting up transients
  Vec phs = { 0,    0,    0, tp };
  //Vec phs = { 0,    0,    PI, PI }; // alternating (the last is actualy at index 1, if sorted
  Vec dec = { 1, 1./2,    1, tl };
  //Vec att = 

  // set up modal filter bank:
  int N = ceilInt(length*fs);
  Vec x(N);
  rosic::rsModalFilterBankDD mfb;
  mfb.setSampleRate(fs);
  mfb.setNonLinearFeedback(nonlin);
  mfb.setReferenceFrequency(f);
  mfb.setReferenceDecay(decay);
  mfb.setReferenceAttack(attack);
  mfb.setModalParameters(frq, amp, 0.1*dec, dec, phs);

  // synthesize signal as impulse-response of modal filter bank:
  x[0] = mfb.getSample(1);
  for(int n = 1; n < N; n++)
    x[n] = mfb.getSample(0);
  RAPT::rsArrayTools::normalize(&x[0], N, 1.0, true);

  // plot and/or write to wavefile:
  //plotImpulseResponse(mfb, 5000, 1.0);
  plotFrequencyResponse(mfb, 5000, 0.0, 400.0, fs, false);
  rosic::writeToMonoWaveFile("ModalWithTransient.wav", &x[0], N, (int)fs);
    // todo: make convenience function that takes a std::vector as signal, double for the 
    // sample-rate

  // maybe have a feedback-drive and feedback-amount factor

  // using x / (1 + x^(2*N)) with large N as nonlinearity gives a sort of sawtooth shape 
  // (one cycle only)

  // idea for more controlled feedback:
  // -evaluate the complex frequency response (magnitude and phase) H(z) of the modal filter bank
  //  without feedback (just add up all individual complex(!) contributions)
  // -figure out at which w the phase response exp(-j*w) * H(exp(j*w)) goes through 180° (the first
  //  exp(-j*w) is the phase response coming from the unit delay in the feedback path)
  //  ...hmm...i'm not sure, if a closed form solution for w is possible
  // -this is our would-be resonance frequency (...or maybe there are several?)
  // -now put in an allpass into the feedback path that let's us shift this resonance frequency to
  //  any desired value (which, as fraction/multiple of the fundamental, is a user parameter)
  // -evaluate the magnitude of H(z) at that resonance frequency and use a normalized feedback 
  //  parameter fb, such that the net loop gain becomes 1, when fb=1
  // -all this is very similar to the resonance tuning of the ladder
  // -when this is done, introduce a nonlinearity (1/(1+x^2) or whatever) into the feedback path
  // -maybe the Q of the feedback allpass can be adjusted too (maybe it makes the resonance more
  //  narrow or broad?)
  // -maybe make the nonlinearity adjustable (by pre/post multiplication by a and 1/a for some a)
  // -maybe it's necessary to factor in the allpass response right from the beginning and solve the
  //  full response (including allpasses) for w
  // -or maybe w has to be found numerically (bisection or whatever)

  // let's define:
  // -i-th modal transfer function:             Gi(z) 
  // -modal bank transfer function:             G(z) = sum_i Gi(z) 
  // -feedback allpass transfer function:       A(z) = ?...to be found
  // -total feedback transfer function:         F(z) = G(z) * A(z) * (z^-1)
  // -complete transfer function with feedback: H(z) = G(z) / (1 + k * F(z)) 
  //  where k is the feedback gain
  // -todo: evaluate that and plot, first evaluate G(z) - that should somehow tell us, how to set 
  //  up A(z) in order to achieve a 180° phase shift at some particular frequency and |G(z)| at 
  //  that frequency should tell us how to set up k to achieve a certain feedback amount
  // -maybe experiment with other types of filters in the feedback path (lowpass, highpass, etc.)
  //  -especially a highpass could be interesting because transients are characterized by high
  //   frequencies...but maby also low/high shelvers...whatever - just make sure, the phase 
  //   response comes out right and compensate for amplitude losses/gains in feedback factor
  //  -maybe just use a general cookbook biquad with all available filter-types there
  //  -maybe later extend it to a biquad chain a la EasyQ
  // -maybe the feedback gain should be defined in terms of a decay-time - the transient decay 
  //  time
  // -maybe use a modulatable version and use the (summed, delayed) output signal for 
  //  phase-modulation
  // -modulatable version would lend itself well for a monophonic synthesizer with glide

  // todo: 
  //  -set up test project where we put the feedback path externally around the object - easier
  //   for experimentation (shorter build times)
  //  -maybe make a wrapper class rsModalBankWithFeedback
  //  -

  // Observations:
  // -when all modes are in phase, the phase response is alternating between 0 and -180 and at the 
  //  modes it is close to 90° (but not exactly..or maybe it is exact? - hard to say due to messed 
  //  up axis values in plot -> figure out)
  // -when mode phases are alternating, the phase response is monotonic. the 1st mode is at -90°
  //  and the other ones are at distances -n*180 form that (quite precisely)

  // ...maybe make an experiment with just 2 modes to figure out these things ...then maybe 10 for
  // nice plots - figure also out, how the responses depend on attack
  //

  int dummy = 0;
}

/** Energy of the envelope given by a weighted sum of 4 exponential decays:
f(t) = A*exp(a*t) + B*exp(b*t) - C*exp(c*t) - D*exp(d*t)
all uppercase coeffs are assumed to be positive and lowercase coeffs negative */
template<class T>
T fourExpEnergy(T A, T a, T B, T b, T C, T c, T D, T d)
{
  //// Sage code that was used to find this monster formula:
  //var("t A B C D a b c d")
  //f(t) = A*exp(a*t) + B*exp(b*t) - C*exp(c*t) - D*exp(d*t)
  //e(t) = (f(t))^2
  //assume(A > 0)
  //assume(B > 0)
  //assume(C > 0)
  //assume(D > 0)
  //assume(a < 0)
  //assume(b < 0)
  //assume(c < 0)
  //assume(d < 0)
  //assume(exp(b+a)-1 < 0)
  //energy = integral(e(t), t, 0, oo)

  T A2 = A*A,  B2 = B*B,  C2 = C*C,  D2 = D*D;
  T A3 = A*A2, B3 = B*B2, C3 = C*C2, D3 = D*D2;
  T A4 = A*A3, B4 = B*B3, C4 = C*C3, D4 = D*D3;

  T a2 = a*a,  b2 = b*b,  c2 = c*c,  d2 = d*d;
  T a3 = a*a2, b3 = b*b2, c3 = c*c2, d3 = d*d2;
  T a4 = a*a3, b4 = b*b3, c4 = c*c3, d4 = d*d3;

  T E = -1.0/2.0*((D2*a3*b2 + D2*a2*b3)*c4 + (C2*a3*b2 + C2*a2*b3 + (B2*a2 + A2*b2 
    + (A2 + 4*A*B + B2)*a*b)*c3 + (B2*a3 + A2*b3 + (A2 + 4*A*B + 2*B2 - 4*(A + B)*C + C2)*a2*b 
    + (2*A2 + 4*A*B + B2 - 4*(A + B)*C + C2)*a*b2)*c2 + ((B2 - 4*B*C + C2)*a3*b 
    + (A2 + 4*A*B + B2 - 4*(A + B)*C + 2*C2)*a2*b2 + (A2 - 4*A*C + C2)*a*b3)*c)*d4 
    + (D2*a4*b2 + 2*D2*a3*b3 + D2*a2*b4)*c3 + (C2*a4*b2 + 2*C2*a3*b3 + C2*a2*b4 
    + (B2*a2 + A2*b2 + (A2 + 4*A*B + B2)*a*b)*c4 + (2*B2*a3 + 2*A2*b3 
    + (2*A2 + 8*A*B + 4*B2 - 4*(A + B)*C + C2 - 4*(A + B - C)*D + D2)*a2*b 
    + (4*A2 + 8*A*B + 2*B2 - 4*(A + B)*C + C2 - 4*(A + B - C)*D + D2)*a*b2)*c3 
    + (B2*a4 + A2*b4 + (A2 + 4*A*B + 4*B2 - 4*(A + 2*B)*C + 2*C2 - 4*(A + B - C)*D + D2)*a3*b 
    + 2*(2*A2 + 6*A*B + 2*B2 - 6*(A + B)*C + 2*C2 - 4*(A + B - C)*D + D2)*a2*b2 
    + (4*A2 + 4*A*B + B2 - 4*(2*A + B)*C   + 2*C2 - 4*(A + B - C)*D + D2)*a*b3)*c2 
    + ((B2 - 4*B*C + C2)*a4*b + (A2 + 4*A*B + 2*B2 - 4*(A + 2*B)*C 
    + 4*C2 - 4*(A + B - C)*D + D2)*a3*b2 + (2*A2 + 4*A*B + B2 - 4*(2*A + B)*C 
    + 4*C2 - 4*(A + B - C)*D + D2)*a2*b3 + (A2 - 4*A*C + C2)*a*b4)*c)*d3 
    + (D2*a4*b3 + D2*a3*b4)*c2 + (C2*a4*b3 + C2*a3*b4 + (B2*a3 + A2*b3 
    + (A2 + 4*A*B + 2*B2 - 4*(A + B)*D + D2)*a2*b + (2*A2 + 4*A*B + B2 - 4*(A + B)*D + D2)*a*b2)*c4 
    + (B2*a4 + A2*b4 + (A2 + 4*A*B + 4*B2 - 4*(A + B)*C + C2 - 4*(A + 2*B - C)*D + 2*D2)*a3*b 
    + 2*(2*A2 + 6*A*B + 2*B2 - 4*(A + B)*C + C2 - 2*(3*A + 3*B - 2*C)*D + 2*D2)*a2*b2 
    + (4*A2 + 4*A*B + B2     - 4*(A + B)*C + C2 - 4*(2*A + B - C)*D 
    + 2*D2)*a*b3)*c3 + ((2*B2 - 4*B*C + C2 - 4*(B - C)*D + D2)*a4*b 
    + 2*(A2 + 4*A*B + 2*B2 - 2*(2*A + 3*B)*C + 2*C2 - 2*(2*A + 3*B - 3*C)*D + 2*D2)*a3*b2 
    + 2*(2*A2 + 4*A*B + B2 - 2*(3*A + 2*B)*C + 2*C2 - 2*(3*A + 2*B - 3*C)*D + 2*D2)*a2*b3 
    +   (2*A2 - 4*A*C + C2 - 4*(A - C)*D + D2)*a*b4)*c2 
    + ((B2 - 4*B*C + 2*C2 - 4*(B - C)*D + D2)*a4*b2 
    + (A2 + 4*A*B + B2 - 4*(A + B)*C + 4*C2 - 4*(A + B - 2*C)*D + 2*D2)*a3*b3 
    + (A2 - 4*A*C + 2*C2 - 4*(A - C)*D + D2)*a2*b4)*c)*d2 + (((B2 - 4*B*D + D2)*a3*b 
    + (A2 + 4*A*B + B2 - 4*(A + B)*D + 2*D2)*a2*b2 + (A2 - 4*A*D + D2)*a*b3)*c4 
    + ((B2 - 4*B*D + D2)*a4*b + (A2 + 4*A*B + 2*B2 - 4*(A + B)*C + C2 - 4*(A + 2*B - C)*D 
    + 4*D2)*a3*b2 + (2*A2 + 4*A*B + B2 - 4*(A + B)*C + C2 - 4*(2*A + B - C)*D 
    + 4*D2)*a2*b3 + (A2 - 4*A*D + D2)*a*b4)*c3 + ((B2 - 4*B*C + C2 - 4*(B - C)*D + 2*D2)*a4*b2 
    + (A2 + 4*A*B + B2 - 4*(A + B)*C + 2*C2 - 4*(A + B - 2*C)*D + 4*D2)*a3*b3 
    + (A2 - 4*A*C + C2 - 4*(A - C)*D + 2*D2)*a2*b4)*c2 + ((C2 + 4*C*D + D2)*a4*b3 
    + (C2 + 4*C*D + D2)*a3*b4)*c)*d)
    /
    (((a2*b + a*b2)*c3 + (a3*b + 2*a2*b2 + a*b3)*c2 
    + (a3*b2 + a2*b3)*c)*d4 + ((a2*b + a*b2)*c4 + 2*(a3*b + 2*a2*b2 + a*b3)*c3 
    + (a4*b + 4*a3*b2 + 4*a2*b3 + a*b4)*c2 + (a4*b2 + 2*a3*b3 + a2*b4)*c)*d3 
    + ((a3*b + 2*a2*b2 + a*b3)*c4 + (a4*b + 4*a3*b2 + 4*a2*b3 + a*b4)*c3 
    + 2*(a4*b2 + 2*a3*b3 + a2*b4)*c2 + (a4*b3 + a3*b4)*c)*d2 + ((a3*b2 + a2*b3)*c4 
    + (a4*b2 + 2*a3*b3 + a2*b4)*c3 + (a4*b3 + a3*b4)*c2)*d);

  // frequent patterns: 4*(A + B)*C, 4*(A - C)*D, 4*(A + B - C)*D

  return E;
}

void fourExponentials()
{
  // f(t) = A*exp(a*t) + B*exp(b*t) - C*exp(c*t) - D*exp(d*t)
  // A,a,B,b are reposnibel for the decay, C,c,D,d for the attack

  double A, B, C, D, a, b, c, d;

  A = 0.2;  a = -1/10.0;  // late decay, a = -1/tau_a, 
  B = 1-A;  b = -1/2.0;   // early decay
  C = 0.8;  c = -1/3.0;
  D = 1-C;  d = -1/0.2;

  double E = fourExpEnergy(A, a, B, b, C, c, D, d); 
  double normalizer = 1/sqrt(E);

  // create time axis and evaluate function, evaluate also it's square and plot it
  double tMin = 0, tMax = 30;
  static const int N = 500;
  double t[N], env[N], env2[N], energy[N]; // energy is accumulated energy up to t
  RAPT::rsArrayTools::fillWithRangeLinear(t, N, tMin, tMax);
  for(int n = 0; n < N; n++) {
    env[n]  = A*exp(a*t[n]) + B*exp(b*t[n]) - C*exp(c*t[n]) - D*exp(d*t[n]);
    env[n] *= normalizer; // energy normalization
    env2[n] = env[n]*env[n];
  }

  RAPT::rsNumericIntegral(t, env2, energy, N);
  // if the envelope is decayed to zero at the end, the final value in the energy array should be
  // close to 1 - the accumulated energy from 0 to infinity of the whole (normalized) envelope 
  // should be 1

  // plot envelope and its square:
  GNUPlotter plt;
  plt.addDataArrays(N, t, env, env2, energy);
  plt.plot();
 
  // the additional parameters are actually quite useful to shape the envelope - but i still need
  // to figure out how to parameterize it in the most user-friendly way. maybe EarlyDecay should
  // be a fraction of LateDecay and the EarlyAttack a fraction of LateAttack
  // EarlyAttack, LateAttack, AttackBlend - same for Decay - maybe have an overall time-scale
}

void modalWithFancyEnv()
{
  // Here, we create a mode with an envelope consisting of 4 exponentials

  // user parameters:
  double length      = 10.0;
  double sampleRate  = 48000;

  //double frequency   = 220;
  double frequency   = 100*GOLDEN_RATIO;
  double amplitude   = 2.5;
  double phase       = 0.0;

  double attackEarly = 0.01;  // maybe rename to attackFast/Slow
  double attackBlend = 0.7;
  double attackLate  = 0.3;

  double decayEarly  = 0.2;
  double decayBlend  = 0.1;
  double decayLate   = 2.0;

  // derived parameters and abbreviations:
  double f = frequency;
  double w = 2*PI*frequency/sampleRate;  // omega
  double A = amplitude;                  // include energy normalizer later
  double p = phase, fs = sampleRate;     // shortcuts

  rsModalFilterDD f0, f1, f2, f3;
  f0.setModalParameters(f, -(1-attackBlend)*A, attackEarly, p, fs);
  f1.setModalParameters(f, -   attackBlend *A, attackLate,  p, fs);
  f2.setModalParameters(f,  (1-decayBlend) *A, decayEarly,  p, fs);
  f3.setModalParameters(f,     decayBlend  *A, decayLate,   p, fs);

  // synthesize the sound:
  int numSamples = (int) ceil(length*sampleRate);
  std::vector<double> x(numSamples);
  x[0] = f0.getSample(1) + f1.getSample(1) + f2.getSample(1) + f3.getSample(1);
  for(int n = 1; n < numSamples; n++)
    x[n] = f0.getSample(0) + f1.getSample(0) + f2.getSample(0) + f3.getSample(0);

  // now synthesize the sound again using the rsModalFilterFloatSSE2 class - the result may be
  // slightyl different due to single precision processing:

  rsModalFilterFloatSSE2 mf;
  mf.setParametersTwoEnvs(w, amplitude, p, 
    fs*attackEarly, fs*attackLate, attackBlend,
    fs*decayEarly,  fs*decayLate,  decayBlend);
  std::vector<double> y(numSamples); // we convert the floats back to double on the fly
  y[0] = mf.getSample(1.f);
  for(int n = 1; n < numSamples; n++)
    y[n] = mf.getSample(0.f);

  // experiment: add a signal with a frequency slightly offset but otherwise the same parameters
  // to get mode-beating:
  mf.reset();
  double df = 7.0;
  double dw = 2*PI*df/fs;
  mf.setParametersTwoEnvs(w+dw, amplitude, p, 
    fs*attackEarly, fs*attackLate, attackBlend,
    fs*decayEarly,  fs*decayLate,  decayBlend);
  std::vector<double> z(numSamples); // we convert the floats back to double on the fly
  z[0] = mf.getSample(1.f) + y[0];
  for(int n = 1; n < numSamples; n++)
    z[n] = 0.125*mf.getSample(0.f) + y[n];
  // question what is the perceived frequency as function of the two beating mode freqs and 
  // amplitudes - when the amplitudes are equal, it's the arithemtic mean but if they are 
  // unequal, it's supposed to be skewed toward the louder mode - but by what function? maybe
  // just a weighted arithmetic mean? that would be simple enough and make sense - if we want to 
  // provide a mode-beating parameter at some stage, we need a formula to adjust the master
  // frequency according to the delta-f. 
  // maybe have a beating parameter b in 0..1 and do
  // freq1 = masterFreq -    b  * beatFreq
  // freq2 = masterFreq + (1-b) * beatFreq
  // out = (1-b)*sineWithFreq1 + b*sineWithFreq2
  // ...or something
  // we can also adjust the phase of the second wave...but then it more and more turns into a 
  // full blown additional mode with its full parameter set...hmmm
  // ...maybe we should have a bunch of simple modes and another bunch of more complex modes and 
  // the user can choose to model a given mode with whatever version of the modal filter he wants
  // or better: just auto-detect, if a given feature is used or not and allocate the appropriate
  // mode filter for it ...maybe especially the faster decaying modes do not need such a fancy
  // envelope?

  // ToDo next: perfomrance tests and tweaking to get the best possible performance, try:
  // -addition instead of subtraction in difference equation
  // -(transposed) direct form 2
  // -make the y a member variable...hmm - we actually already have y1 - oh - no - that
  //  can't be used because it's nnede later to update y2



  // compute error due to single precision floating point precision in optimized filter:
  std::vector<double> err(numSamples);
  for(int n = 0; n < numSamples; n++)
    err[n] = x[n] - y[n];

  //GNUPlotter plt;
  //plt.addDataArrays(5000, &err[10000]);
  //plt.plot();
  // the error is also in the shape of a sinewave with amplitude of order 0.01

  // maybe normalize the error for writing to wavefile:
  RAPT::rsArrayTools::normalize(&err[0], numSamples, 1.0);
  // absolute error is greatest when signals are loudest ...roughly - maybe consider relative
  // error...at least, the error doesn't seem to get much worse over time maybe try frequency
  // that is "more irrational" ...however, all in all, it looks good - single precision seems
  // to be good enough for rock'n'roll

  rosic::writeToMonoWaveFile("ModalWithFancyEnvDbl.wav", &x[0],   numSamples, (int)fs);
  rosic::writeToMonoWaveFile("ModalWithFancyEnvFlt.wav", &y[0],   numSamples, (int)fs);
  rosic::writeToMonoWaveFile("ModalWithFancyEnvErr.wav", &err[0], numSamples, (int)fs);
  rosic::writeToMonoWaveFile("ModalWithFancyEnvBeating.wav", &z[0],   numSamples, (int)fs);
}
/*
// move to plotting tools:
void stemPlot(int N, double *x, double *y)
{
  GNUPlotter plt;
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(N, x, y); // can probably be done without adding the data twice
  plt.setGraphStyles("impulses", "points pt 7 ps 1.2");
  plt.plot();
}
*/
void plotModalLevelSpectrum(const rosic::rsModalSynth& ms, int key, int vel)
{
  int mL = ms.getLowestMode();
  int mH = ms.getHighestMode();
  int M = mH - mL + 1;
  std::vector<double> f(M), a(M);
  for(int m = mL-1; m < mH; m++) {
    int i = m-mL+1;
    f[i]  = ms.getModeFreqRatio(m, key, vel);  // relative frequency
    a[i]  = ms.getModeLevel(    m, key, vel);  // decay time in milliseconds
  }
  stemPlot(M, &f[0], &a[0]);
}

// todo: merge bothe functions - too much repetition!
void plotModalDecaySpectrum(const rosic::rsModalSynth& ms, int key, int vel)
{
  int mL = ms.getLowestMode();
  int mH = ms.getHighestMode();
  int M = mH - mL + 1;
  std::vector<double> f(M), d(M);
  for(int m = mL-1; m < mH; m++) {
    int i = m-mL+1;
    f[i]  = ms.getModeFreqRatio(m, key, vel);  // relative frequency
    d[i]  = ms.getModeDecay(    m, key, vel);  // decay time in milliseconds
  }
  stemPlot(M, &f[0], &d[0]);
}

void modalSynthSpectra()
{
  // A testbed to set up various parameters in rosic::rsModalSynth and plot resulting spectra in 
  // order to check that the implementations produces the expected results.

  // set up the modal synth:
  rosic::rsModalSynth ms;

  // freq ratio parameters:
  ms.setFreqRatioProfileTopLeft(    ms.ALL_HARMONICS);
  ms.setFreqRatioProfileTopRight(   ms.ALL_HARMONICS);
  ms.setFreqRatioProfileBottomLeft( ms.ALL_HARMONICS);
  ms.setFreqRatioProfileBottomRight(ms.ALL_HARMONICS);
  ms.setFreqRatioMixX(0.0);            // -1..+1, 0 means eqaul mix
  ms.setFreqRatioMixY(0.0);   
  ms.setInharmonicity(0.0);            // relevant only for stiff string tuning

  // global parameters:
  ms.setSampleRate(44100);
  ms.setLowestMode(1);
  ms.setHighestMode(32);

  ms.setLevel(0.0);                    // in dB
  ms.setLevelByKey(0.0);               // in +- dB at extreme keys
  ms.setLevelByVel(0.0);               // in +- dB at extreme velocities

  // amplitude spectrum parameters:
  ms.setAmpSlope(-2.0);                // in dB/oct
  ms.setAmpSlopeByKey(-3.0);           // in +- dB/oct at extreme keys
  ms.setAmpSlopeByVel(+4.0);           // in +- dB/oct at extreme velocities

  // envelope parameters:
  ms.setAttack(10.0);                  // in ms
  ms.setAttackByRatio(-50);            // in %
  ms.setAttackByKey(-50);              // in %
  ms.setAttackByVel(-50);              // in %
  ms.setDecay(1000);                   // in ms
  ms.setDecayByRatio(-50);             // in %
  ms.setDecayByKey(-50);               // in %
  ms.setDecayByVel(-50);               // in %

  // plot:
  plotModalLevelSpectrum(ms, 64, 64);
  //plotModalDecaySpectrum(ms, 64, 64);

  // todo: figure out why there are comb filtering artifacts on retrigger
}

void setupHarmonicAnalyzerForModal(RAPT::rsHarmonicAnalyzer<double>& analyzer, double fs)
{
  typedef rsWindowFunction::WindowType WT;
  analyzer.setSampleRate(fs);
  analyzer.setSincInterpolationLength(64);
  analyzer.setNumCyclesPerBlock(4);
  //analyzer.setWindowType(WT::hamming);
  analyzer.setWindowType(WT::blackman);
  analyzer.setSpectralOversampling(8);  // zero padding
  analyzer.setAllowInharmonics(true);
  analyzer.setSpectralPeakSearchWidth(0.5);       // default: 1 - blackman needs a value less than 1
  analyzer.setMinPeakToMainlobeWidthRatio(0.75);  // default: 0.75
}

void modalDecayFit1()
{
  // Computes the A,tau values for an exponential decay function f(t) = A * exp(-t/tau) given two
  // time instants t1,t2 and associated amplitudes a1,a2. Then it creates the actual decay function
  // and plots it.

  int    N  = 2000;    // number of samples
  double fs = 10000;   // sample rate
  double t1 = 0.04;    // first time instant
  double a1 = 0.7;     // first amplitude value
  double t2 = 0.14;    // second time instant
  double a2 = 0.3;     // second amplitude value

  // compute parameters for exponential:
  double A, tau;
  rsExpDecayParameters(t1, a1, t2, a2, &A, &tau);

  // create and plot exponential decay:
  typedef std::vector<double> Vec;
  Vec t(N), x(N);
  RAPT::rsArrayTools::fillWithIndex(&t[0], N);
  t = (1.0/fs) * t;
  for(int n = 0; n < N; n++)
    x[n] = A * exp(-t[n]/tau);
  rsPlotVectorsXY(t, x);  // OK - looks good - goes through (0.04, 0.7) and (0.14, 0.3)
}

void modalDecayFit2()
{
  // We create attack/decay envelope with non-equidistant sample instants and from that data, we 
  // try to estimate/recover the parameters (attack, decay, amplitude). Attack and amplitude are
  // obtained from the location and height of the peak, decay is obtained from the total energy 
  // (or maybe total area - we'll see which gives better fits on real-world data - with this 
  // artificial data, both variants should give the same result). The plan is to use this 
  // algorithm for automatic determination of optimal settings for the peak-shadowing algorithm...

  // user parameters:
  int    N      = 1000; // number of samples
  double length = 2.5;  // length of signal in seconds
  double att    = 0.1;  // attack in seconds
  double dec    = 0.5;  // decay in seconds
  double amp    = 1.5;  // peak amplitude

  // compute algo parameters
  double tau1 = dec, tau2, scl;
  expDiffScalerAndTau2(tau1, att, &tau2, &scl);
  // amp scales the env such that the peak height is 1.0 - todo: use another scaler
  // also let the suer choose the total length (in seconds)

  // compute envelope:
  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;
  Vec t = randomSampleInstants(N, 0.2, 1.8, 0);  // maybe have a randomness parameter
  //Vec t = randomSampleInstants(N, 1.0, 1.0, 0);
  AT::scale(&t[0], N, length/t[N-1]);
  Vec x(N);
  int n;
  for(n = 0; n < N; n++)
    x[n] = amp * scl * (exp(-t[n]/tau1) - exp(-t[n]/tau2));

  // todo: 
  // -recover att, dec, amp from data in t,x
  //  -find peak location - that fixes att and amp
  //  -find area or energy - that fixes dec - maybe energy is better because it uses the squares 
  //   which makes the result less susceptible to cutting off the tail

  int i = AT::maxIndex(&x[0], N);
  double att2 = t[i];  // todo: refine by quadratic interpolation (maybe)
  double amp2 = x[i];

  // for estimating dec, there are several possibilities:
  // -use the total area under the curve (from 0 to the end)
  // -use the area only from the peak to the end
  // -we could use the area-formula for the exponential decay or for an actual attack/decay curve
  // -we could use the energy instead of the area under the curve (i.e. square the curve before
  //  taking the integral)
  // -maybe try them all and see, which approach gives the best results - maybe try them with real
  //  world signals, too
  // -test the decay formula also with the exact values of att and amp (not the recovered ones) to 
  //  decouple the decay-error from the att/amp error in this experiment

  // create the squared signal (for energy computation):
  Vec x2(N);
  for(n = 0; n < N; n++)
    x2[n] = x[n] * x[n];

  // Compute total and partial areas and energies by numerical integration:
  double E_t, E_p, A_t, A_p;
  Vec tmp(N);    
  rsNumericIntegral(&t[0], &x[0],  &tmp[0], N);   A_t = tmp[N-1];
  rsNumericIntegral(&t[i], &x[i],  &tmp[i], N-i); A_p = tmp[N-1];
  rsNumericIntegral(&t[0], &x2[0], &tmp[0], N);   E_t = tmp[N-1]; // total energy (from 0 to the end)
  rsNumericIntegral(&t[i], &x2[i], &tmp[i], N-i); E_p = tmp[N-1]; // partial energy (from i to the end)
  // todo: get rid of the tmp array - write a numerical integration routine that just returns a 
  // single value instead of filling an array - we only need the final value here

  double dec2_Ap = A_p / amp;
  // not that bad! 0.5025... correct value is 0.5 - only 5% error

  double dec2_At = A_t / amp;
  // this is actually the wrong area formula - if we use the total energy, we should use a formula 
  // that takes into account the nonzero attack - the result is far too big

  double dec2_Ep = 2*E_p / (amp*amp);
  // this is very wrong - check derivation - it should be better, i think

  double dec2_Et = 2*E_t / (amp*amp);
  // again, this is actually the wrong formula - use one that takes into account the attack - this
  // result is even more off

  // todo: figure out correct formulas that take attack into account - we may also have to estimate
  // the "scl" - but that is coupled with the decay
  // -figure out, how the error behaves when increasing the attack - it goes up, as expected

  // Observations:
  // -it's a bit surprising that dec2_Ep is less accurate than dec2_Ap - i would have expected that
  //  the squaring in the energy computation makes the calculation more robust with respect to 
  //  cutting off a portion from the tail 
  // -when increasing the attack, the error in the decay estimate dec2_Ap goes up - which is 
  //  totally what is expected because this formula is based on assuming a zero attack, decaying
  //  exponential (starting at the peak-instant with peak amplitude). When attack is increased,
  //  the shape becomes more rounded at the start, increasing the overall area and thereby the
  //  decay estimate
  //  -maybe it's possible to figure out a function, how to re-adjust the decay estimate in terms
  //   of the attack estimate (plot relative error of decay estimate as function of attack, fit 
  //   some function to it and use that for compensation)



  // plot:
  //rsPlotVectorsXY(t, x);
  rsPlotVectorsXY(t, x, x2);
}

void modalDecayFit()
{
  //modalDecayFit1();
  modalDecayFit2();
}



void modalAnalysis1()
{
  //double f   = 1000;
  int N  = 5000;
  int fs = 22050;
 
  rsModalFilterParameters<double> p;
  p.freq  = 300;
  p.amp   = 0.5;
  p.phase = 45.0;
  p.att   = 0.02;
  p.dec   = 0.2;

  double peak = p.att * fs; // time-instant of the peak

  typedef std::vector<double> Vec;

  Vec x = synthesizeModal(p, fs, N);

  // create a sinusoidal model and resynthesize sinusoidally:
  RAPT::rsHarmonicAnalyzer<double> sineAnalyzer;
  setupHarmonicAnalyzerForModal(sineAnalyzer, fs);
  RAPT::rsSinusoidalModel<double> sineModel = sineAnalyzer.analyze(&x[0], N);
  sineModel.keepOnly({1});
  //plotSineModel(sineModel, fs);
  Vec ys = synthesizeSinusoidal(sineModel, fs);


  rsModalAnalyzer<double> modeAnalyzer;
  //std::vector<rsModalFilterParameters<double>> modeModel
  //  = modeAnalyzer.getModalModel(sineModel);
  rsModalFilterParameters<double> p2 = modeAnalyzer.getModalModel(sineModel.getPartial(0));

  Vec ym = synthesizeModal(p2, fs, N);



  //rsPlotVector(x);
  rsPlotVectors(x, ys, ym);
  //rosic::writeToMonoWaveFile("ModalOriginal.wav", &x[0], N, fs);
}

void modalAnalysisPluck()
{
  // We create a plucked type of sound using modal synthesis, then analyze that sound via the 
  // sinusoidal analysis framework and then further analyze the output of the sinusoidal analysis
  // in terms of modal parameters. The goal is to recover the original modal parameters as
  // accurately as possible. If we are able to correctly recover modal parameters, we may also
  // apply the same analysis algorithm to natural sounds.

  // user parameters:
  double key = 65;
  double sampleRate = 44100;
  int length = 44100;

  // create the input signal:
  typedef std::vector<double> Vec;
  Vec x = createModalPluck(key, sampleRate, length);

  // create a sinusoidal model and resynthesize sinusoidally:
  RAPT::rsHarmonicAnalyzer<double> sineAnalyzer;
  setupHarmonicAnalyzerForModal(sineAnalyzer, sampleRate);
  RAPT::rsSinusoidalModel<double> sineModel = sineAnalyzer.analyze(&x[0], length);
  sineModel.removePartialsAbove(66);  // model for key=64 shows partials up to 40 kHz - see below
  sineModel.removePartial(0);         // DC confuses modal model
  //plotSineModel(sineModel, sampleRate);
  Vec ys = synthesizeSinusoidal(sineModel, sampleRate);

  // create a modal model and resynthesize modally:
  rsModalAnalyzer<double> modeAnalyzer;
  std::vector<rsModalFilterParameters<double>> modeModel
    = modeAnalyzer.getModalModel(sineModel);
  Vec ym = synthesizeModal(modeModel, sampleRate, length);
  // -decay is still only coarsly estimated
  // -phases seem to be wrong - use phase from peak-amp (or somewhere else after the transient)
  //  and obtain start-phase from that (phase data at the start may be inaccurate)
  // -todo: try it first with a simpler sound with only one mode - tweak the phase estimation
  //  algo until it works with a single-mode sound, then try again with the pluck

  rosic::writeToMonoWaveFile("ModalPluckOriginal.wav", &x[0],  length, (int)sampleRate);
  rosic::writeToMonoWaveFile("ModalPluckSinusoidal.wav", &ys[0], (int)ys.size(), (int)sampleRate);
  rosic::writeToMonoWaveFile("ModalPluckModal.wav", &ym[0], (int)ym.size(), (int)sampleRate);

  // Observations:
  // -with key=64, the model contains partials with frequencies up to 40kHz - this is because the 
  //  maximum measured cycle-length is around 135 and the pitch-flattening step stretches this out 
  //  to the next power of two which is 256 - so the cycle-lenth that goes into the analysis is 
  //  almost twice the original cycle length, leading to a nomial doubling of the bandwidth - the 
  //  amplitudes of the partials above fs/2 are numerically close to zero, though
  //  -maybe the analyzer should restrict itself to analyze only frequencies that fall into the
  //   range up to the original fs/2, not to the upsampled fs/2
  // -with key=65, the maximum measured cycle-length is 127.73 and the pitch-flattener selects 128
  //  as its target length
}


// maybe move this before modalAnalysisPluck - it's sort of a preliminary - it analyzes only a 
// single partial
void modalPartialResynthesis() // maybe rename to exponentialTailModeling
{
  // Tries to resynthesize a partial via an exponentially decaying sinusoid. The difference to the
  // function above is that we use the original phase of the analyzed partials. This is mostly for 
  // Elan to splice an exponential decay to some sample where the  initial section is taken from
  // the original - maybe that means, we should match amplitude and phase exactly at the 
  // splice-point, which the user can select

  // user parameters:
  double key = 65;
  double sampleRate = 44100;
  int N = 44100;              // length in samples
  int partialIndex = 2;       // partial index on which the analysis/resynthesis is tested
  int spliceIndex  = 200;     // frame index, at which we want to match amplitude and phase:

  // create the input signal:
  typedef std::vector<double> Vec;
  Vec x = createModalPluck(key, sampleRate, N);

  // create the sinusoidal model:
  RAPT::rsHarmonicAnalyzer<double> sineAnalyzer;
  setupHarmonicAnalyzerForModal(sineAnalyzer, sampleRate);
  RAPT::rsSinusoidalModel<double> sineModel = sineAnalyzer.analyze(&x[0], N);

  // pick the partial with which we do our experiments and extract envelopes:
  RAPT::rsSinusoidalPartial<double> partial = sineModel.getPartial(partialIndex);
  Vec timeArray  = partial.getTimeArray();
  Vec ampArray   = partial.getAmplitudeArray();
  Vec freqArray  = partial.getFrequencyArray();
  Vec phaseArray = partial.getPhaseArray();
  int numFrames  = (int) timeArray.size();
  //rsPlotVectorsXY(timeArray, ampArray);    // plot amplitude envelope


  int maxIndex  = RAPT::rsArrayTools::maxIndex(&ampArray[0], numFrames);

  // estimate decay time and amplitude of an exponential decay that fits the amp-env:
  double A, tau;
  //rsExpDecayParameters(numFrames, &timeArray[0], &ampArray[0], maxIndex, spliceIndex, &A, &tau);
  rsExpDecayParameters(timeArray[maxIndex], ampArray[maxIndex], 
    timeArray[spliceIndex], ampArray[spliceIndex], &A, &tau); // verify, if called correctly


  // create an exponential decay with given A,tau - the time axis is obtained the analysis frames:

  Vec a2(numFrames);
  for(int k = 0; k < numFrames; k++)
    a2[k] = A * exp(-timeArray[k] / tau);

  // plot actual amplitude envelope and the exponential decay that is supposed to fit it:
  //rsPlotVectorsXY(timeArray, ampArray, a2);
  // ok - this looks ok - maybe it could be tweaked to match the original envelope even better by
  // choosing the first match point not exactly at the peak but some time later (because later, the
  // influence of the attack exponential has decayed away and doesn't disturb the exponential decay
  // anymore)

  Vec tail = rsExpDecayTail(partial, spliceIndex, sampleRate);
  //rsPlotVector(tail);

  // generate the partial from the sinusoidal model and write both into a stereo wavefile for 
  // comparison:
  sineModel.keepOnly({ (size_t)partialIndex });
  Vec ys = synthesizeSinusoidal(sineModel, sampleRate);


  double spliceTime = timeArray[spliceIndex];
    // this is the time-instant (in seconds) at which instantaneous phase and amplitude of the 
    // exponentially decaying sinusoid will match the values from the sine-model - for splicing
    // together original and tail, it's best to use a crossfade centered at that instant

  rosic::writeToMonoWaveFile("PartialResynth.wav", &ys[0],   (int)ys.size(),   (int)sampleRate);
  rosic::writeToMonoWaveFile("PartialExpTail.wav", &tail[0], (int)tail.size(), (int)sampleRate);
    // in the audiofiles, zoom in to the spliceTime to verify that both signals indeed do match
    // phase- and amplitude-wise at that instant


  // todo: maybe instead of reading out the amplitude at two selected points, estimate A,tau from
  // the energy-decay-relief - see section 2.2.5 here:
  // https://www2.ak.tu-berlin.de/~akgroup/ak_pub/abschlussarbeiten/2014/KaapPascal_MasA.pdf

  // todo: in addition to estimate decay and amplitude, also estimate attack....actually, attack is 
  // the only thing that is missing to make a full modal model of the partial



  // todo: estimate attack, decay, amplitude from (t,a)

  // create a modal model of the partial:
  //rsModalAnalyzer<double> modeAnalyzer;
  //rsModalFilterParameters<double> modalParams = modeAnalyzer.getModalModel(partial);
  // params are frequency, phase, amplitude, attack, decay - in this context, only 
  // amplitude, attack, decay are relevant - actually only amplitude and decay because the portion,
  // before but it's nicer to include the attack - hmm - maybe for the splicing of an exponential 
  // decay, this is not what we need....maybe just work with the t,a,f,p arrays here directly
 

  /*
  std::vector<rsModalFilterParameters<double>> modeModel
    = modeAnalyzer.getModalModel(sineModel);
  Vec ym = synthesizeModal(modeModel, sampleRate, length);
  */

  //rosic::writeToMonoWaveFile("ModalPluckOriginal.wav", &x[0],  length, (int)sampleRate);
}


/*
Ideas:
-user defines modal parameters for various keys and velocities (at least two keys at two 
 velocities)
-for key/vel in between the defined ones, an appropriate interpolation is used (maybe
 log -> natural cubic spline -> exp), outside the range, extrapolation is used - for the spline, 
 it may make sense to just use the a0, a1 polynomial coeffs of the endpoints - a2 is zero anyway 
 due to the "natural" end conditions and in extrapolation we artificially set a3 also to zero in 
 order to prevent the polynomial from going crazy
-maybe we can emulate mode-beating by using two modes of nearby frequencies (not necessarily with 
 the same envelope/amplitude
-have an SFZ export function - user selects keyrange, numKeysPerSample, NumVelsPerSample
*/