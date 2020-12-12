#ifndef RAPT_MODALFILTERBANK_H
#define RAPT_MODALFILTERBANK_H

//=================================================================================================
// utility functions for damped-sine filter design:
// maybe move them as static functions into class rsModalFilter or rsModalFilterWithAttack

/** Designs (i.e. computes coefficients from user parameters for) a two-pole-one-zero filter in 
terms of its impulse response which is given as the damped sinusoid:

h[n] = A * exp(-n/d) * sin(w*n + p) * u[n]

where u[n] is the unit step function. The filter can be implemented as:

y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] - a2*y[n-2]

The parameters are:
w: normalized radian frequency (= 2*pi*f/fs, f: frequency, fs: samplerate)
A: amplitude
d: normalized decay time constant (= tau*fs, tau: time (in s) to decay to A/e = A*0.3678...)
p: start phase in radians.

You can use different datatypes for the parameters and coefficients (for example, double and float)

move to FilterDesignFormulas, change sign convention for a-coeffs */
template<class TPar, class TCof>
void rsDampedSineFilterCoeffs(
  TPar w, TPar A, TPar d, TPar p, TCof* b0, TCof* b1, TCof* a1, TCof* a2);
// have inverse function rsDampedSineFilterParams




/** Finds the value k such that k*exp(-c*k) - exp(-c) = 0  for the given parameter c in the range
0 < c < 1. k = 1 is always a solution to this root finding problem - the function will try to
find the other solution, if any - otherwise it will return 1. */
template<class T>
T findDecayScalerLess1(T c);
// \todo make a function findDecayScalerGreater1 which assumes c > 1 (or c >= 1) and a function
// findDecayScaler which dispatches between the two versions - check the old matlab files - iirc,
// we need a different initial guess in this case and/or a different iteration algorithm 
// (fixed-point- instead of newton-raphson or something...) - but maybe move as static function
// into the class that needs this function

/** Given a desired asymptotic decay time constant tau1 and the desired time instant of the peak
tPeak, this function computes the required other time constant tau2 and a scale factor to be used
for an attack-decay envelope that can be obtained by the scaled difference of two exponential
decays, like so:

  x(t) = scaler * ( exp(-t/tau1) - exp(-t/tau2) ) 

For this to work, we must have tPeak < tau1 because otherwise the decay will also be governed by
the second term (and the computational methods, we use here may not converge because they also
rely on this assumption). */
template<class T>
void expDiffScalerAndTau2(T tau1, T tPeak, T* tau2, T* scaler);

//=================================================================================================

/** Implements a two-pole filter (without any zeros) decribed by the difference equation:
y[n] = x[n] - a1*y[n-1] - a2*y[n-2]

\todo maybe make a subclass rsNormalizedTwoPoleFilter that uses a b0 (or g) coefficient which
      implements the gain
\todo maybe implement a rotating-phasor version
*/

template<class TSig, class TPar>
class rsTwoPoleFilter
{

public:

  /** Constructor. Creates a neutral filter. */
  rsTwoPoleFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the filter coefficients from a desired normalized radian frequency w = 2*PI*f/fs
  (f: resonance frequency, fs: samplerate) and a normalized decay time-constant d = tau*fs
  (tau: time constant in seconds, time to deacy to 1/e = 0.36..). */
  void setFrequencyAndDecay(TPar w, TPar d);

  /** Sets a gain-factor for the output. */
  void setOutputGain(TPar newGain);


  /** UNDER CONSTRUCTION */
  void setFrequencyAndAbsoluteBandwidth(TPar w, TPar b)
  {
    setFrequencyAndDecay(w, TPar(1)/b);
  }
  // todo: check, if we need a proportionality factor other than 1, i.e. c/b instead of 1/b for 
  // some c, bandwidth definition should be based on the two -3.01 dB points left and right from 
  // the peak...i think it could be 2*Pi or something

  /** UNDER CONSTRUCTION */
  void setFrequencyAndRelativeBandwidth(TPar w, TPar b)
  {
    setFrequencyAndAbsoluteBandwidth(w, w/b);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the magnitude response of the filter at the given normalized radian frequency 
  excluding the output gain as set by setOutputGain. */
  TPar getMagnitudeAt(TPar w);


  std::complex<TPar> getTransferFunctionAt(std::complex<TPar> z)
  {
    return biquadTransferFunctionAt(TPar(1), TPar(0), TPar(0), a1, a2, z);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Computes an output sample from an input sample. */
  RS_INLINE TSig getSample(TSig in);

  /** Resets the internal state variables to zero. */
  void reset();



protected:

  TPar a1, a2;   // feedback coeffs
  TPar g;        // output gain
  TSig y1, y2;   // past ouputs

};

template<class TSig, class TPar>
RS_INLINE TSig rsTwoPoleFilter<TSig, TPar>::getSample(TSig in)
{
  TSig y = in - a1*y1 - a2*y2;
  y2 = y1;
  y1 = y;
  return g*y;
}

//=================================================================================================

/** This class implements a filter that realizes an exponentially decaying sinusoid that may
represent a single mode of vibration in some sound.

\todo templatize this class for different types of input signals (for mono, stereo, multichannel)

\todo write another implementation based on a complex-valued one-pole:
z[n] = a*z[n-1] + x[n];
y[n] = cr*z.re  + ci*z.im;
or, in code:
z = a*z + x;
return cr*z.re + ci*z.im;

where "a" is a complex coefficient that determines frequency and damping, "z" is the filter's
state as complex number (representing a spiraling phasor) and cr, ci determine start-phase and
amplitude of the oscillation. such an implemention will lead itself much better to frequency- and
phase modulation because whenever "a" changes, the phasor will just continue to rotate from where
it was, just changing frequency and possibly damping. also the computation of the coeffs might be
simpler. on the other hand, running the filter might be more expensive (but we need to measure,
if this is actually the case - if not, the current implementation would actually be obsolete)

maybe include a conversion of the decay-time into a decay-rate - if tau is the decay-time, the
decay-rate would be 20/(tau*ln(10)) dB/s (decibels per second) - verify this formula

*/

template<class TSig, class TPar>
class rsModalFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsModalFilter();


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets all the mode parameters - this triggers a calculation of the filter coefficients.
  Parameters: frequency is in Hz, amplitude is a raw factor, decayTime in seconds, 
  startPhase in degrees and sampleRate also in Hz. */
  void setModalParameters(TPar frequency, TPar amplitude, TPar decayTime, TPar startPhase, 
    TPar sampleRate);

  /** Copies the filter coefficients g, b1, a1, a2 from another instance into this one. */
  void copyCoefficientsFrom(const rsModalFilter &other);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the decay time constant (tau) in seconds. Because the filter object is agnostic of
  the sample-rate at which it runs (for economy reasons), you must pass it along to the
  function. */
  TPar getDecayTime(TPar sampleRate);

  /** Returns the filter's z-domain transfer function value at the given value of z. */
  std::complex<TPar> getTransferFunctionAt(std::complex<TPar> z);

  /** Returns the magnitude response of the filter at the normalized radian frequency
  w = 2*pi*f/fs. */
  TPar getMagnitudeAt(TPar w);


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one sample at a time. */
  RS_INLINE TSig getSample(TSig in);


  void processBlock(TSig in[], TSig out[], int blockSize);

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the internal state. */
  void reset();


protected:

  /** \name Data */

  TSig x1, y1, y2;      /** filter state variables */
  TPar b0, b1, a1, a2;  /** filter parameters */

  //double b1, a1, a2, g; /** filter parameters */ // old

};


//=================================================================================================

/**

Alternative implementation based on a rotating phasor in the complex plane (incomplete)

This implementation allows for some interesting extensions to the original modal filter, like
amplitude- and phase-modulation by various linear and nonlinear combinations of the real and
imaginary part. Phase-modulation by the imaginary part introduces overtones, by selecting the
imaginary part alone, we leave the start phase unmodified. Phase modulation by the instantaneous
envelope (== absoulte value == sqrt(re*re+im*im), or some power thereof, leads to a sweep of the
filters frequency - the frequency is coupled to the amplitude. The chosen power determines the
speed of the sweep [check this].

ideas for user parameters (aside from frequency, amplitude, phase, decay):
FfPmA: feedforward phase-modulation amount (modulation index)
FfPmR: feedforward phase-modulation ratio (carrier/modulator)
FbPmA: feedback phase-modulation amount
FbPmR: feedback phase-modulation ratio
sweep: causes a frequency sweep by feedforward phase-modulation by (a power of) the instantaneous
       envelope

call this approach "Nonlinear Modal Synthesis" (NLM Synthesis)

\todo: implement a version of this with attack parameter

todo:
-test freq-mod, phase-mod, amp-mod and decay-mod by external signal (use a long decay time, 
 feed impulse, switch values during impulse-response)
-check, what happens to the modulation response, if zeros and/or gain are processed before the 
 poles - it will probably make the impulse-response unmodulatable (for phase and amp) during it 
 runs - but what happens, if there's a continuous input signal like a saw wave?

note: the nice thing about such a system as compared to a more general additive synthesis
approach is that it is an autonomous system, that is, it does not need to know what time it is
and thus depends only on the input signal and not explictly on time.

decayRate in dB/s: decayRate = 20/(tau*ln(10))  ...check this

// factor out ModalFilterViaPhasor as alternative implementation of the linear variant, move the
nonlinear extenstions into a subclass, maybe create a simple Modal-PM module consisting of two
phasor filters that phase-modulate each other

*/

template<class TSig, class TPar>
class rsNonlinearModalFilter
{

public:


  rsNonlinearModalFilter();

  void setModalParameters(TPar frequency, TPar amplitude, TPar decayTime,
    TPar startPhase, TPar sampleRate);

  void setAmplitude(TPar newAmplitude);

  void setPhaseModulation(TPar newPhaseModulation);

  void copyCoefficientsFrom(const rsNonlinearModalFilter &other);

  //std::complex<TPar> getTransferFunctionAt(std::complex<TPar> z);




  /** Produces the complex output sample for a given complex input sample. A linear combination
  of the real and imaginary part can be used to calculate a real output sample from that, but you
  can also calculate the instantaneous envelope and/or phase. */
  RS_INLINE std::complex<TSig> getComplexSample(std::complex<TSig> in);

  RS_INLINE TSig getSample(TSig in);


  /** Calculates the instantaneous envelope of the mode (== absoulte value == sqrt(re*re+im*im).
  This is useful for the analysis of the modal envelope in input sounds that should be modeled.
  \todo maybe factored out into the linear superclass.*/
  RS_INLINE TSig getInstantaneousEnvelope(TSig in);

  void reset();

protected:

  std::complex<TPar> a;   // recursion-coeff
  std::complex<TSig> z;   // state
  TPar cr, ci;            // output weights for z.re, z.im
  TPar amplitude;         // overall amplitude
  TPar startPhase;        // start-phase in degrees
  TPar phaseModByIm;      // phase-modulation by imaginary part
  TPar phaseModByAbs;     // phase-modulation by absolute value (== instantaneous envelope)
                          // ...maybe use re*re+im*im instead of the sqrt thereof (faster decay)
                          // ...or some general power p (re*re+im*im)^p

  // maybe introduce amplitude- (and/or ring-) modulation by imaginary and/or real part,
  // phase/amp modulation by angle, maybe by other quantities derived by nonlinear combinations
  // of  the real and imaginary part (for example re*im)

};

//===============================================================================================

/**

This class implements a filter that realizes an impulse response that is a sinusoid that is
enveloped by an attack/decay envelope.

\todo: move the ModalFilter classes into the directory with the filters

// is this a Gammatone filter? http://en.wikipedia.org/wiki/Gammatone_filter

*/

template<class TSig, class TPar>
class rsModalFilterWithAttack
{

public:

  /** \name Setup */

  /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
  void setModalParameters(TPar frequency, TPar amplitude, TPar attackTime,
    TPar decayTime, TPar startPhase, TPar sampleRate, TPar detuneFactor = 1.0);


  /** \name Inquiry */

  /** Returns the filter's z-domain transfer function value at the given value of z. */
  std::complex<TPar> getTransferFunctionAt(std::complex<TPar> z);

  /** Returns the time (in seconds) it takes for the sound to decay down to the given
  decayLevel. */
  TPar getLength(TPar decayLevel, TPar sampleRate);


  /** \name Audio Processing */

  /** Calculates one sample at a time. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal state. */
  void reset();


protected:

  /** \name Data */

  rsModalFilter<TSig, TPar> modalFilter1, modalFilter2;

  // \todo maybe introduce a delayline here - this allows the attack-phase to start at some
  // specified delayed time instant - only one delayline is needed to feed the two modal filters
  // ...we may need a templatized delayline class as well
  // ...or make a delegating wrapper-class for stereo signals - allows to detune left and right
  // and some other stuff L/R stuff(phase, amplitude)

};

//===============================================================================================

/**

This class implements a filter that realizes an impulse response that is a sinusoid that is
enveloped by an attack/decay envelope. In contrast to the other implementation, this filter
combines the two parallel decaying sinusoids into a single 4th order filter (implemented in
direct from 2) - hopefully this leads to a performance gain. If not, it may nevertheless lend
itself better to manual optimization via SSE2.

*/

template<class TSig, class TPar>
class rsModalFilterWithAttack2
{

public:

  /** \name Construction  */

  /** Default constructor. */
  rsModalFilterWithAttack2();


  /** \name Setup */

  /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
  void setModalParameters(TPar frequency, TPar amplitude, TPar attackTime,
    TPar decayTime, TPar startPhase, TPar sampleRate);


  /** \name Inquiry */

  /** Returns the time (in seconds) it takes for the sound to decay down to the given
  decayLevel. */
  //double getLength(double decayLevel, double sampleRate);


  /** \name Audio Processing */

  /** Calculates one sample at a time. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the internal state. */
  void reset();

protected:

  /** \name Data */

  TSig w1, w2, w3, w4;   /** filter state variables */
  TPar a1, a2, a3, a4;   /** feedback coefficients */
  TPar b1, b2, b3;       /** feedforward coefficients (b0 is zero) */

};

//===============================================================================================

/**

This is a bank (i.e. parallel connection) of modal filters, each with its own set of parameters.

*/

template<class TSig, class TPar>
class rsModalFilterBank
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsModalFilterBank();

  /** Destructor. */
  ~rsModalFilterBank();


  /** \name Setup */

  /** \todo: use const references for the parameter vectors */

  /** Sets up the sample-rate fopr this filter. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the reference frequency with respect to which all the absolute mode frequencies are
  computed like absoluteModeFrequency = referenceFrequency * relativeModeFrequency. For harmonic
  sounds, this will be the fundamental frequency. */
  void setReferenceFrequency(TPar newFrequency);

  void setReferenceAttack(TPar newAttack);
  void setReferenceDecay(TPar newDecay);

  /** Sets the strength of the nonlinear feedback. This parameter is important to shape the 
  transient. */
  void setNonLinearFeedback(TPar newFeedback) { nonLinFeedback = newFeedback; }



  /** Sets the frequencies of the modes. */
  //void setModalFrequencies(std::vector<TPar> newFrequencies);
    // use a reference for the parameter

  /** Sets the amplitudes of the modes. */
  //void setModalAmplitudes(std::vector<TPar> newAmplitudes);

  /** Sets the decay time-constants of the modes. */
  //void setModalDecayTimes(std::vector<TPar> newDecayTimes);

  /** Sets the start-phases of the modes (in degrees). */
  //void setModalStartPhases(std::vector<TPar> newStartPhases);

  /** Sets all mode parameters at once. */
  void setModalParameters(
    std::vector<TPar> newFrequencies, 
    std::vector<TPar> newAmplitudes,
    std::vector<TPar> newAttackTimes, 
    std::vector<TPar> newDecayTimes, 
    std::vector<TPar> newStartPhases);


  /** \name Inquiry */

  /** Returns the filter's z-domain transfer function value at the given value of z. */
  std::complex<TPar> getTransferFunctionAt(std::complex<TPar> z);

  /** Returns the length (in seconds) until the produced sound will have decayed to the specified
  level. */
  TPar getLength(TPar decayLevel);

  /** Returns the currently active number of modal filters. */
  RS_INLINE int getNumModes()
  {
    int M = (int)rsMin((size_t)numModes, frequencies.size(), amplitudes.size(), decayTimes.size());
    return (int)rsMin((size_t)M, startPhases.size());
    // M: number of modes - optimize this, use a member variable
  }



  /** \name Audio Processing */

  /** Calculates one sample at a time. You should pass an exitatition signal in the input
  argument - if this is a unit impulse, the response will be a superposition of decaying
  sinusoids with parameters (frequencies, amplitudes, decay-times and start-phases) as specified
  via the respective set-functions. */
  RS_INLINE TSig getSample(TSig in);

  /** Processes a block of samples. The inputs are passed as pointer via the in[] parameter,
  likewise for the output block. The pointers must be distinct. */
  void processBlock(TSig in[], TSig out[], int blockSize);


  /** \name Misc */

  /** Resets the internal states of the filters. */
  void reset();

  void calculateModalFilterCoefficients();


  /** \name Static member functions  */
  // useful for setting up vectors of modal parameters - these should go into a class
  // rsModalParameterGenerator


  /** Gives the relative mode decay time for mode with relative frequency f given a (relative)
  cutoff frequency fc and an ultimate slope of the decay-time function (with respect to f)
  determined by the exponent p
  d(f) = a / (b + (f/fc)^p) where a and b are adjusted such that d(f=1)=1 and d(f=fc)=1/2
  fc must be > 1 */
  static TPar modeDecayTime(TPar f, TPar fc, TPar p);

  /** Creates a vector with pseudo-random start-phases subject to the constraint that the first
  sample in the generated sound should be at zero. You must pass a vector of modal amplitudes
  that correspond to the phases to be created. */
  static std::vector<TPar> randomModePhases(const std::vector<TPar> &amplitudes, 
    TPar randomness = 1.0, int seed = 0);

  /** Creates a vector of decay-times for the modes that have relative frequecies given in the
  vector f. The rule for the decay time for a mode with relative frequency f is:
  d(f) = a / (b + (f/fc)^p) where the constants a and b are adjusted such that
  d(f=1)=1 and d(f=fc)=1/2, -> fc > 1 must be satisfied (hmm - maybe fc != 1 suffices) */
  static std::vector<TPar> modeDecayTimes(std::vector<TPar> f, TPar fc, TPar p);

  /** Scales the values at a given interval by the given scaler starting at the given startIndex.
  This is useful for scaling the decay times of even or odd harmonics, for example. */
  static std::vector<TPar> scaleAtIntervals(std::vector<TPar> v, int startIndex, int interval,
    TPar scaler);

protected:

  /** Feedback saturation function. */
  RS_INLINE TSig saturate(TSig x);

  /** \name Data */

  static const int maxNumModes = 1000;    // get rid of this - allow an arbitrary number

  std::vector<rsModalFilterWithAttack<TSig, TPar>> modalFilters;    // maybe use c-arrays instead


  //  consolidate into class rsModalBankParameters:
  TPar referenceFrequency;
  TPar referenceDecay;
  TPar referenceAttack;
  std::vector<TPar> frequencies;
  std::vector<TPar> amplitudes;
  std::vector<TPar> attackTimes;
  std::vector<TPar> decayTimes;
  std::vector<TPar> startPhases;
    // maybe use f, g, a, d, p
  int    numModes;     // restricts the number of modes to be generated - if -1, there's no
                       // restriction other than the minimum of the dimensionalities of the
                       // parameter vectors
  TPar sampleRate;

  TPar nonLinFeedback = TPar(0);
  TSig out = 0;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsModalFilter<TSig, TPar>::getSample(TSig in)
{
  TSig y = b0*in + b1*x1 - a1*y1 - a2*y2; // todo: use all plusses (more efficient)
  x1 = in;  // maybe multiply by b0 at the output instead of input for better reponse to amplitude
  y2 = y1;  // modulation
  y1 = y;
  return y;
}

template<class TSig, class TPar>
RS_INLINE std::complex<TSig> rsNonlinearModalFilter<TSig, TPar>::getComplexSample(std::complex<TSig> in)
{
  z = a*z + in;
  return amplitude * z;
}

template<class TSig, class TPar>
RS_INLINE TSig rsNonlinearModalFilter<TSig, TPar>::getSample(TSig in)
{
  z = a*z + in;
  return amplitude * (cr*z.real() + ci*z.imag());  // test (all is linear)

  // maybe have a 3rd coeff cIn and use
  // amplitude * (cr*z.real() + ci*z.imag() + cIn*in);
  // maybe this filter could replace a general biquad?
  // and/or maybe we should have a parametrized injection injection
  // z.re += a * in; z.im += b * in instead of injecting into the real part only
  // maybe, if a,b = 1/sqrt(2) (equal injection into both parts), we could also model two real
  // poles in parallel (maybe we would then have to compensate in the phase parameter by 
  // subtracting 45° when computing cr, ci



  /*
  bool linear = false;
  if( linear == true )
    return amplitude * (cr*z.re + ci*z.im); // cr, ci are fixed
  else
  {
    // obtain a modulator signal that rotates with a different speed (multiplied by 1/ratio where
    // ratio = c/m - the usual carrier/modulator ratio in phase-modulation synthesis
    //rsComplexDbl zm = z;
    rsComplexDbl zm = z*z;  // later: rsExp(z, 1/ratio); // ratio = carrierFreq/modFreq
    //rsComplexDbl zm = z*z*z;

    double phi = (PI/180.0)*startPhase + phaseModByIm*zm.im;
    cr = sin(phi);
    ci = cos(phi);
    // later: rsSinCos((PI/180.0)*startPhase + phaseMod*z.im, &cr, &ci);

    //...but this is only feedforward phase-modulation of the output - to obtain feedback
    // phase-modulation, we need to multiply z by exp(i*pm*z.im)
    // z *= expC(rsComplexDbl(0.0, 1.0) * fbPhaseModByIm * z.im) ...or something

    return amplitude * (cr*z.re + ci*z.im);
  }
  */

}

template<class TSig, class TPar>
RS_INLINE TSig rsNonlinearModalFilter<TSig, TPar>::getInstantaneousEnvelope(TSig in)
{
  z = a*z + in;
  return amplitude * abs(z);
}

template<class TSig, class TPar>
RS_INLINE TSig rsModalFilterWithAttack<TSig, TPar>::getSample(TSig in)
{
  return modalFilter1.getSample(in) - modalFilter2.getSample(in);
}

template<class TSig, class TPar>
RS_INLINE TSig rsModalFilterWithAttack2<TSig, TPar>::getSample(TSig in)
{
  TSig w0 = in - a1*w1 - a2*w2 - a3*w3 - a4*w4; // todo: use positive sign convention
  TSig y  =      b1*w1 + b2*w2 + b3*w3;
  w4 = w3;
  w3 = w2;
  w2 = w1;
  w1 = w0;
  return y;
  //return b1*w2 + b2*w3 + b3*w4; // can also be used
}

template<class TSig, class TPar>
RS_INLINE TSig rsModalFilterBank<TSig, TPar>::saturate(TSig x)
{
  // later: x = (1-k)*x + k*x^3; k between 0..1, high k make small feedback values even smaller, so
  // it sort of controls the decay-shape of the nonlinear feedback

  //return rsClip(x, TSig(-1), TSig(+1));  // use a parameteric softClip later
  //return rsClip(x, TSig(-0.1), TSig(+0.1));

  //return x / (1 + abs(x)); // cheap, saturating, not very smooth at origin

  x *= x;  // test
  x *= x;

  return x = x / (1 + x*x); // non-monotonic, smooth (continuous derivatives of all order (verify))
}

template<class TSig, class TPar>
RS_INLINE TSig rsModalFilterBank<TSig, TPar>::getSample(TSig in)
{
  //int M = (int)rsMin((size_t)numModes, frequencies.size(), amplitudes.size(), decayTimes.size());
  //M     = (int)rsMin((size_t)M, startPhases.size());
  //// M: number of modes - optimize this, use a member variable
  //// factor out inot getNumModes

  int M = getNumModes();

  in += saturate(nonLinFeedback*out);
  out = 0.0;
  for(int m = 0; m < M; m++)
    out += modalFilters[m].getSample(in); // factor out into getBody, make a similar loop for getTransient
  return out;

  // old:
  //double out = 0.0;
  //for(int m = 0; m < M; m++)
  //  out += modalFilters[m].getSample(in);
  //return out;

  // idea for nonlinear feedback: apply a nonlinearity (saturation or similar) to the final sum
  // and feed it back to the input with some feedback factor. OR: obtain the signal to be fed back
  // as a differently weighted sum of the modal filter outputs, in which the weights depend 
  // (inversely) on the decay time of the respective mode - faster decaying modes produce higher
  // nonlinear feedback -> should emphasize/chaosify the transient more than the body

  // maybe introduce special "transient" modes - with inharmonic frequencies and fast decay

  // maybe the nolinearity could also have some sigmoid(x^3) chracteristic -> lessen the feedback for 
  // low-level signals

  // maybe do all of these things in ModalFilterNonLinear...or actually, we need a 
  // ModalFilterBankNonlinear ......for this, we should also have a nonlinear/modulatable version
  // with attack...check performance of plain modal filter vs the modulatable version - if the
  // modulatable has similar performance, there's actually no need to use the plain version

  // maybe to further shape the transient, we could use attack/decay lowpass filters - maybe let's
  // have a second array of filters for this purpose - this has parameters attack/decay/amplitude
  // per "mode" (amplitude can be negative) - so transient building blocks are simpler than body
  // building blocks (which additionally have frequency and phase)
  
  // ...OR use a completely general pole/zero 
  // filter to model the transient

  // check, if the transient behaves similar at different sample-rates (due to chaos, it may not)

  // give the user a transient/body balance slider

  // make a stereo version - channels may differ in start-phase and/or frequency (maybe provide
  // a mid/side adjustment afterwards)

  // allow the user to use a combination of (filtered) impulse and (filtered) noise as input
  // signal

  // the velocity should control the transient strength (among other things)

  // provide a dedicated API for setting up the transient
}

#endif