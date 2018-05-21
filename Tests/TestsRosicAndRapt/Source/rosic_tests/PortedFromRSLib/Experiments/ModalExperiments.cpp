#include "ModalExperiments.h"

void modalFilter()
{
  static const int N = 5000;  // number of samples to plot

  double fs  = 44100;  // samplerate in Hz
  double td  = 0.1;    // decay time constant in seconds
  double f   = 100;    // frequency in Hz
  double phs = 45;     // phase in degrees
  double A   = 1.0;    // amplitude as raw factor

  // create and set up the modal filter:
  rsModalFilterDD mf;
  mf.setModalParameters(f, A, td, phs, fs);

  // generate time-axis and impulse-response:
  double t[N], x[N];
  createTimeAxis(N, t, fs);
  getImpulseResponse(mf, x, N);

  // plot the impulse-response (versus the time-axis):
  plotData(N, t, x);
}

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
  double length = 2.0;    // total length in seconds

  double attack = 0.1;  // master attack
  double decay  = 0.3;  // master decay time in seconds
  double nonlin = 0.0;  // nonlinear feedback strength

  // transient parameters:
  double tr = 0.5 * (1+sqrt(5.0)); // golden ratio - relative frequency of "transient mode"
  double ta = 1.0;                 // transient strength/amplitude
  double tl = 0.2;                 // transient length
  double tp = 0.0;                 // transient phase

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
  Vec dec = { 1, 1./2,    1, tl };
  //Vec att = 

  // set up modal filter bank:
  int N = ceilInt(length*fs);
  Vec x(N);
  rosic::rsModalFilterBankDD mfb;
  mfb.setSampleRate(fs);
  mfb.setReferenceFrequency(f);
  mfb.setReferenceDecay(decay);
  mfb.setReferenceAttack(attack);
  mfb.setModalParameters(frq, amp, 0.1*dec, dec, phs);

  // synthesize signal as impulse-response of modal filter bank:
  x[0] = mfb.getSample(1);
  for(int n = 1; n < N; n++)
    x[n] = mfb.getSample(0);
  RAPT::rsArray::normalize(&x[0], N, 1.0, true);

  // plot or write to wavefile:
  rosic::writeToMonoWaveFile("ModalWithTransient.wav", &x[0], N, (int)fs);
    // todo: make convenience function that takes a std::vector as signal, double for the 
    // sample-rate


  int dummy = 0;
}