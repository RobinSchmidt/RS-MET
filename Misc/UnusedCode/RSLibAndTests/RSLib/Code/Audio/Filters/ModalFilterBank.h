#ifndef RS_MODALFILTERBANK_H
#define RS_MODALFILTERBANK_H

namespace RSLib
{


/** Implements a two-pole filter (without any zeros) decribed by the difference equation: 
y[n] = x[n] - a1*y[n-1] - a2*y[n-2]

\todo maybe make a subclass rsNormalizedTwoPoleFilter that uses a b0 (or g) coefficient which 
      implements the gain
\todo maybe implement a rotating-phasor version
*/

class RSLib_API rsTwoPoleFilter
{

public:

  /** Constructor. Creates a neutral filter. */
  rsTwoPoleFilter();

  /** Sets up the filter coefficients from a desired normalized radian frequency w = 2*PI*f/fs 
  (f: resonance frequency, fs: samplerate) and a normalized decay time-constant d = tau*fs 
  (tau: time constant in seconds, time to deacy to 1/e = 0.36..). */
  void setFrequencyAndDecay(double w, double d);

  /** Sets a gain-factor for the output. */
  void setOutputGain(double newGain);

  /** Returns the magnitude response of the filter at the given normalized radian frequency. */
  double getMagnitudeAt(double w);

  /** Resets the internal state variables to zero. */
  void reset();

  /** Computes an output sample from an input sample. */
  RS_INLINE double getSample(double in);

protected:

  double a1, a2;   // feedback coeffs
  double g;        // output gain
  double y1, y2;   // past ouputs

};

RS_INLINE double rsTwoPoleFilter::getSample(double in)
{
  double y = in - a1*y1 - a2*y2;
  y2 = y1;
  y1 = y;
  return g*y;
}

//=================================================================================================


  /** Designs a two-pole one-zero filter in terms of its impulse response which is given as the
  damped sinusoid:

  h[n] = A * exp(-n/d) * sin(w*n + p) * u[n]

  where u[n] is the unit step function. The filter can be implemented as:

  y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] - a2*y[n-2]

  The parameters are:
  w: normalized radian frequency (= 2*pi*f/fs, f: frequency, fs: samplerate)
  A: amplitude
  d: normalized decay time constant (= tau*fs, tau: time (in s) to decay to A/e = A*0.3678...)
  p: start phase in radians. */
  void RSLib_API rsDampedSineFilter(double w, double A, double d, double p, double *b0, 
    double *b1, double *a1, double *a2);

  /** This is actually almost the same as rsDampedSineFilter, just that here the phase is given in
  degrees. 
  \todo get rid of that function
  */
  //void RSLib_API modalParamsToFilterCoeffs(double omega, double amplitude, double tau,
  //  double startPhase, double *g, double *b1, double *a1, double *a2);

  /** Finds the value k such that k*exp(-c*k) - exp(-c) = 0  for the given parameter c in the range
  0 < c < 1. k = 1 is always a solution to this root finding problem - the function will try to
  find the other solution, if any - otherwise it will return 1. */
  double RSLib_API findDecayScalerLess1(double c);

  // \todo make a function findDecayScalerGreater1 which assumes c > 1 (or c >= 1) and a function
  // findDecayScaler which switches between the two versions

  /** Given a desired asymptotic decay time constant tau1 and the desired time instant of the peak
  tPeak, this function computes the required other time constant tau2 and a scale factor to be used
  for an attack-decay envelope that can be obtained by the scaled difference of two exponential
  decays, like so:
  x(t) = scaler * ( exp(-t/tau1) - exp(t/tau2) )  \todo use latex-markup
  For this to work, we must have tPeak < tau1 because otherwise the decay will also be governed by
  the second term (and the computational methods, we use here may not converge because they also
  rely on this assumption). */
  void RSLib_API expDiffScalerAndTau2(double tau1, double tPeak, double *tau2, double *scaler);


  //===============================================================================================

  /**

  This class implements a filter that realizes an exponentially decaying sinusoid that may
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

  class RSLib_API rsModalFilter
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsModalFilter();


    /** \name Setup */

    /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
    void setModalParameters(double frequency, double amplitude, double decayTime,
      double startPhase, double sampleRate);

    /** Copies the filter coefficients g, b1, a1, a2 from another instance into this one. */
    void copyCoefficientsFrom(const rsModalFilter &other);


    /** \name Inquiry */

    /** Returns the decay time constant (tau) in seconds. Because the filter object is agnostic of
    the sample-rate at which it runs (for economy reasons), you must pass it along to the
    function. */
    double getDecayTime(double sampleRate);

    /** Returns the magnitude response of the filter at the normalized radian frequency 
    w = 2*pi*f/fs. */
    double getMagnitudeAt(double w);



    /** \name Audio Processing */

    /** Calculates one sample at a time. */
    RS_INLINE double getSample(double in);


    void processBlock(double in[], double out[], int blockSize);


    /** \name Misc */

    /** Resets the internal state. */
    void reset();


  protected:

    /** \name Data */

    double x1, y1, y2;      /** filter state variables */
    double b0, b1, a1, a2;  /** filter parameters */

    //double b1, a1, a2, g; /** filter parameters */ // old

  };


  //===============================================================================================

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

  note: the nice thing about such a system as compared to a more general additive synthesis
  approach is that it is an autonomous system, that is, it does not need to know what time it is
  and thus depends only on the input signal and not explictly on time.

  decayRate in dB/s: decayRate = 20/(tau*ln(10))  ...check this

  // factor out ModalFilterViaPhasor as alternative implementation of the linear variant, move the
  nonlinear extenstions into a subclass, maybe create a simple Modal-PM module consisting of two
  phasor filters that phase-modulate each other

  */

  class RSLib_API rsNonlinearModalFilter
  {

  public:


    rsNonlinearModalFilter();

    void setModalParameters(double frequency, double amplitude, double decayTime,
      double startPhase, double sampleRate);

    void setAmplitude(double newAmplitude);

    void setPhaseModulation(double newPhaseModulation);

    void copyCoefficientsFrom(const rsNonlinearModalFilter &other);

    /** Produces the complex output sample for a given complex input sample. A linear combination 
    of the real and imaginary part can be used to calculate a real output sample from that, but you 
    can also calculate the instantaneous envelope and/or phase. */
    RS_INLINE rsComplexDbl getComplexSample(rsComplexDbl in);

    RS_INLINE double getSample(double in);






    /** Calculates the instantaneous envelope of the mode (== absoulte value == sqrt(re*re+im*im). 
    This is useful for the analysis of the modal envelope in input sounds that should be modeled. 
    \todo maybe factored out into the linear superclass.*/
    RS_INLINE double getInstantaneousEnvelope(double in);

    void reset();

  protected:

    rsComplex<double> a, z;    // recursion-coeff and state
    double            cr, ci;  // output weights for z.re, z.im


    double amplitude;      // overall amplitude
    double startPhase;     // start-phase in degrees
    double phaseModByIm;   // phase-modulation by imaginary part
    double phaseModByAbs;  // phase-modulation by absolute value (== instantaneous envelope)
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

  class RSLib_API rsModalFilterWithAttack
  {

  public:

    /** \name Setup */

    /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
    void setModalParameters(double frequency, double amplitude, double attackTime,
      double decayTime, double startPhase, double sampleRate, double detuneFactor = 1.0);


    /** \name Inquiry */

    /** Returns the time (in seconds) it takes for the sound to decay down to the given
    decayLevel. */
    double getLength(double decayLevel, double sampleRate);


    /** \name Audio Processing */

    /** Calculates one sample at a time. */
    RS_INLINE double getSample(double in);


    /** \name Misc */

    /** Resets the internal state. */
    void reset();


  protected:

    /** \name Data */

    rsModalFilter modalFilter1, modalFilter2;

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

  class RSLib_API rsModalFilterWithAttack2
  {

  public:

    /** \name Construction  */

    /** Default constructor. */
    rsModalFilterWithAttack2();


    /** \name Setup */

    /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
    void setModalParameters(double frequency, double amplitude, double attackTime,
      double decayTime, double startPhase, double sampleRate);


    /** \name Inquiry */

    /** Returns the time (in seconds) it takes for the sound to decay down to the given
    decayLevel. */
    //double getLength(double decayLevel, double sampleRate);


    /** \name Audio Processing */

    /** Calculates one sample at a time. */
    RS_INLINE double getSample(double in);


    /** \name Misc */

    /** Resets the internal state. */
    void reset();

  protected:

    /** \name Data */

    double w1, w2, w3, w4;   /** filter state variables */
    double a1, a2, a3, a4;   /** feedback coefficients */
    double b1, b2, b3;       /** feedforward coefficients (b0 is zero) */

  };





  //===============================================================================================

  /**

  This is a bank (i.e. parallel connection) of modal filters, each with its own set of parameters.

  */

  class RSLib_API rsModalFilterBank
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
    void setSampleRate(double newSampleRate);

    /** Sets the reference frequency with respect to which all the absolute mode frequencies are
    computed like absoluteModeFrequency = referenceFrequency * relativeModeFrequency. For harmonic
    sounds, this will be the fundamental frequency. */
    void setReferenceFrequency(double newFrequency);

    void setReferenceAttack(double newAttack);
    void setReferenceDecay(double newDecay);



    /** Sets the frequencies of the modes. */
    void setModalFrequencies(rsVectorDbl newFrequencies);

    /** Sets the amplitudes of the modes. */
    void setModalAmplitudes(rsVectorDbl newAmplitudes);

    /** Sets the decay time-constants of the modes. */
    void setModalDecayTimes(rsVectorDbl newDecayTimes);

    /** Sets the start-phases of the modes (in degrees). */
    void setModalStartPhases(rsVectorDbl newStartPhases);

    /** Sets all mode parameters at once. */
    void setModalParameters(rsVectorDbl newFrequencies, rsVectorDbl newAmplitudes,
      rsVectorDbl newAttackTimes, rsVectorDbl newDecayTimes, rsVectorDbl newStartPhases);


    /** \name Inquiry */

    /** Returns the length (in seconds) until the produced sound will have decayed to the specified
    level. */
    double getLength(double decayLevel);


    /** \name Audio Processing */

    /** Calculates one sample at a time. You should pass an exitatition signal in the input
    argument - if this is a unit impulse, the response will be a superposition of decaying
    sinusoids with parameters (frequencies, amplitudes, decay-times and start-phases) as specified
    via the respective set-functions. */
    RS_INLINE double getSample(double in);

    /** Processes a block of samples. The inputs are passed as pointer via the in[] parameter,
    likewise for the output block. The pointers must be distinct. */
    void processBlock(double in[], double out[], int blockSize);


    /** \name Misc */

    /** Resets the internal states of the filters. */
    void resetModalFilters();

    void calculateModalFilterCoefficients();


    /** \name Static member functions  */
    //        useful for setting up vectors of modal parameters

    /** Gives the relative mode decay time for mode with relative frequency f given a (relative)
    cutoff frequency fc and an ultimate slope of the decay-time function (with respect to f)
    determined by the exponent p
    d(f) = a / (b + (f/fc)^p) where a and b are adjusted such that d(f=1)=1 and d(f=fc)=1/2
    fc must be > 1 */
    static double modeDecayTime(double f, double fc, double p);

    /** Creates a vector with pseudo-random start-phases subject to the constraint that the first
    sample in the generated sound should be at zero. You must pass a vector of modal amplitudes
    that correspond to the phases to be created. */
    static rsVectorDbl randomModePhases(const rsVectorDbl &amplitudes, double randomness = 1.0,
      int seed = 0);

    /** Creates a vector of decay-times for the modes that have relative frequecies given in the
    vector f. The rule for the decay time for a mode with relative frequency f is:
    d(f) = a / (b + (f/fc)^p) where the constants a and b are adjusted such that
    d(f=1)=1 and d(f=fc)=1/2, -> fc > 1 must be satisfied (hmm - maybe fc != 1 suffices) */
    static rsVectorDbl modeDecayTimes(rsVectorDbl f, double fc, double p);

    /** Scales the values at a given interval by the given scaler starting at the given startIndex.
    This is useful for scaling the decay times of even or odd harmonics, for example. */
    static rsVectorDbl scaleAtIntervals(rsVectorDbl v, int startIndex, int interval,
      double scaler);

  protected:

    /** \name Data */

    static const int maxNumModes = 1000;    // get rid of this - allow an arbitrary number

    std::vector<rsModalFilterWithAttack> modalFilters;    // mabe use c-arrays instead

    //  consolidate into class rsModalBankParameters:
    double referenceFrequency;
    double referenceDecay;
    double referenceAttack;
    rsVectorDbl frequencies;
    rsVectorDbl amplitudes;
    rsVectorDbl attackTimes;
    rsVectorDbl decayTimes;
    rsVectorDbl startPhases;
      // maybe use f, g, a, d, p
    int    numModes;     // restricts the number of modes to be generated - if -1, there's no
                         // restriction other than the minimum of the dimensionalities of the
                         // parameter vectors


    double sampleRate;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE double rsModalFilter::getSample(double in)
  {
    double y = b0*in + b1*x1 - a1*y1 - a2*y2;
    x1 = in;
    y2 = y1;
    y1 = y;
    return y;
  }

  RS_INLINE rsComplexDbl rsNonlinearModalFilter::getComplexSample(rsComplexDbl in)
  {
    z = a*z + in;
    return amplitude * z;
  }

  RS_INLINE double rsNonlinearModalFilter::getSample(double in)
  {
    z = a*z + in;
    return amplitude * (cr*z.re + ci*z.im);  // test (all is linear)


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
    
  RS_INLINE double rsNonlinearModalFilter::getInstantaneousEnvelope(double in)
  {
    z = a*z + in;
    return amplitude * z.getRadius();
  }

  RS_INLINE double rsModalFilterWithAttack::getSample(double in)
  {
    return modalFilter1.getSample(in) - modalFilter2.getSample(in);
  }

  RS_INLINE double rsModalFilterWithAttack2::getSample(double in)
  {
    double w0 = in - a1*w1 - a2*w2 - a3*w3 - a4*w4;
    double y  =      b1*w1 + b2*w2 + b3*w3;
    w4 = w3;
    w3 = w2;
    w2 = w1;
    w1 = w0;
    return y;
    //return b1*w2 + b2*w3 + b3*w4; // can also be used
  }

  RS_INLINE double rsModalFilterBank::getSample(double in)
  {
    int M = rsMin(numModes, frequencies.dim, amplitudes.dim, decayTimes.dim);
    M     = rsMin(M, startPhases.dim);
    // M: number of modes - optimize this, use a member variable


    double out = 0.0;
    for(int m = 0; m < M; m++)
      out += modalFilters[m].getSample(in);
    return out;
  }

}

#endif
