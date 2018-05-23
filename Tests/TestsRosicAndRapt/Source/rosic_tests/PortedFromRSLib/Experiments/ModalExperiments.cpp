#include "ModalExperiments.h"

void modalFilter()
{
  static const int N = 5000;  // number of samples to plot

  double fs  = 44100;  // samplerate in Hz
  double td  = 0.1;    // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = 0;     // phase in degrees
  double A   = 1.0;    // amplitude as raw factor

  // create and set up the modal filter:
  rsModalFilterDD mf;
  mf.setModalParameters(f, A, td, phs, fs);

  //plotImpulseResponse(mf, 5000, 1.0);
  plotFrequencyResponse(mf, 2000, 20.0, fs/2, fs, true);
}

void modalFilterFreqResp()
{
  // filter parameters:
  double fs  = 44100;  // samplerate in Hz
  double frq = 1000;   // frequency in Hz
  double amp = 1.0;    // amplitude as raw factor
  double phs = 45;     // phase in degrees
  double att = 0.01;   // attack time in seconds
  double dec = 0.1;    // decay time constant in seconds


  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(frq, amp, att, dec, phs, fs, 1.0);

  //plotImpulseResponse(mf, 5000, 1.0);
  plotFrequencyResponse(mf, 5000, 20.0, fs/2, fs, true);
}

void modalTwoModes()
{
  double fs  = 44100;
  double at  = 0.01;  // attack - needs to be nonzero for numerical reasons (->try to fix)

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
  //plotImpulseResponse(  mfb, 5000, 1.0);
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

  // create and set up the filter:
  rsModalFilterWithAttackDD mf;
  mf.setModalParameters(f, a, ta, td, phs, fs);

  // generate and plot impulse-response:
  double h[N];
  getImpulseResponse(mf, h, N);
  plotData(N, 0, 1/fs, h);


  int dummy = 0;
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
// this implementation is obsolete


/** Retrieves damped sine filter design parameters from its coefficients. See 
@see rsDampedSineFilter for meaning of parameters. The phase p is returned in the interval 
0...2pi. */
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
  rsDampedSineFilterOld(w, A, d, p, &B0, &B1, &A1, &A2);
  rsDampedSineFilter(   w, A, d, p, &b0, &b1, &a1, &a2);


  // recover design parameters:
  double ww, AA, dd, pp;
  rsDampedSineFilterAnalysis( b0, b1, a1, a2, &ww, &AA, &dd, &pp);
  rsDampedSineFilterAnalysis2(b0, b1, a1, a2, &ww, &AA, &dd, &pp);


  int dummy = 0;
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
  RAPT::rsArray::fillWithZeros(x, N);
  x[0] = 1;
  RAPT::rsArray::filter(x, N, x, N, b, 2, a, 2);

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
  RAPT::rsArray::normalize(&x[0], N, 1.0, true);

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