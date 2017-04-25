#ifndef rosic_ModalSynthesizer_h
#define rosic_ModalSynthesizer_h

#include<vector>

// rosic-indcludes:
#include "../math/rosic_Vector.h"
#include "../math/rosic_Complex.h"
#include "../math/rosic_SpecialFunctionsReal.h"  // why?

namespace rosic
{

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

  */

  class ModalFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModalFilter(); 

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
    void setModalParameters(double frequency, double amplitude, double decayTime, 
      double startPhase, double sampleRate);

    /** Sets the amplitude (initial value of the exponentially decaying envelope). */
    void setAmplitude(double newAmplitude);

    /** Copies the filter coefficients g, b1, a1, a2 from another instance into this one. */
    void copyCoefficientsFrom(const ModalFilter &other);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    double g, b1, a1, a2; // filter parameters
    double x1, y1, y2;    // state variables  

    // maybe introduce members:
    // omega, amplitude, decayTime, startPhase,
    // w (omega), a (amplitude), d (normalized decay), p (start-phase)/
    // using omega and normalized decay avoids the need to store the samplerate - typically we will 
    // use an array of such filters and it would be wasteful if each array-element would store the 
    // samplerate
    // ...but maybe we should maintain these variables in the AudioModule layer - keeps the DSP 
    // layer lean

  };






  // alternative implementation based on a rotating phasor in the complex plane:

  class ModalFilter2
  {

  public:


    ModalFilter2(); 

    void setModalParameters(double frequency, double amplitude, double decayTime, 
      double startPhase, double sampleRate);

    void setAmplitude(double newAmplitude);

    void setPhaseModulation(double newPhaseModulation);

    void copyCoefficientsFrom(const ModalFilter2 &other);

    INLINE double getSample(double in);

    void reset();

  protected:

    Complex a, z;    // recursion-coeff and state
    double  cr, ci;  // output weights for z.re, z.im

    double startPhase;
    double phaseMod;
    double amplitude;

    //double amplitude, phase

  };







  /**

  This class implements a filter that realizes an impulse response that is a sinusoid that is 
  enveloped by an attack/decay envelope.

  \todo: move the ModalFilter classes into the directory with the filters

  */

  class ModalFilterWithAttack
  {

  public:

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets all the mode parameters - this triggers a calculation of the filter coefficients. */
    void setModalParameters(double frequency, double amplitude, double attackTime, 
      double decayTime, double startPhase, double sampleRate, double detuneFactor = 1.0);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset();

    /** Finds the value k such that k*exp(-c*k) - exp(-c) = 0  for the given parameter c. k = 1 is 
    always a solution to this root finding problem - the function will try to find the other 
    solution, if any - otherwise it will return 1. The function assumes 0 < c < 1 */
    static double findDecayScaler(double c);

    /** Given a desired asymptotic decay time constant tau1 and the desired time instant of the 
    peak, this function computes the required other time constant tau2 and a scale factor to be 
    used for an attack-decay envelope that can be obtained by the scaled difference of two 
    exponential decays, like so:
    x(t) = scaler * ( exp(-t/tau1) - exp(t/tau2) )  \todo use latex-markup
    For this to work, we must have peakTime < decayTime because otherwise the decay will also
    be governed by the second term (and the computational methods, we use here may not converge
    because they also rely on this assumption). */
    static void getScalerAndSecondTimeConstant(double decayTime, double peakTime, double &tau2, 
      double &scaler);

    //=============================================================================================

  protected:

    ModalFilter modalFilter1, modalFilter2;

    // \todo maybe introduce a delayline here - this allows the attack-phase to start at some 
    // specified delayed time instant - only one delayline is needed to feed the two modal filters
    // ...we may need a templatized delayline class as well
    // ...or make a delegating wrapper-class for stereo signals - allows to detune left and right
    // and some other stuff L/R stuff(phase, amplitude)

  };


  /**

  This is a ..

  \todo:
  -include a smooth attack to the modes
  -maybe make modal frequency a function of amplitude - modeling of nonlinear effects
  -rename this class into ModalFilterBank
  -use the version of the ModalFilter class that has the attack-parameter

  */

  class ModalSynthesizer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModalSynthesizer(); 

    /** Destructor. */
    ~ModalSynthesizer();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets up the sample-rate fopr this filter. */
    void setSampleRate(double newSampleRate);

    /** Sets the frequencies of the modes. */
    void setModalFrequencies(Vector newFrequencies);

    /** Sets the amplitudes of the modes. */
    void setModalAmplitudes(Vector newAmplitudes);

    /** Sets the decay time-constants of the modes. */
    void setModalDecayTimes(Vector newDecayTimes);

    /** Sets the start-phases of the modes (in degrees). */
    void setModalStartPhases(Vector newStartPhases);

    /** Sets all mode parameters at once. */
    void setModalParameters(Vector newFrequencies, Vector newAmplitudes, 
      Vector newDecayTimes, Vector newStartPhases);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one sample at a time. You should pass an exitatition signal in the input 
    argument - if this is a unit impulse, the response will be a superposition of decaying 
    sinusoids with parameters (frequencies, amplitudes, decay-times and start-phases) as specified 
    via the respective set-functions. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the filters. */
    void resetModalFilters();

    void calculateModalFilterCoefficients();

    //=============================================================================================

  protected:

    static const int maxNumModes = 1000;

    // mabe use c-arrays instead:
    std::vector<ModalFilter> modalFilters;  // use an array of ModalFilterWithAttack
    Vector frequencies;
    Vector amplitudes;
    Vector decayTimes;
    Vector startPhases;
    double sampleRate;
    int    numModes;     // restricts the number of modes to be generated - if -1, there's no 
                         // restriction other than the minimum of the dimensionalities of the
                         // parameter vectors
  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double ModalFilter::getSample(double in)
  {
    double y = in + b1*x1 - a1*y1 - a2*y2;
    x1 = in;
    y2 = y1;
    y1 = y;
    return g*y;
  }

  INLINE double ModalFilter2::getSample(double in)
  {
    z = a*z + in;

    // phase modulation by signal itself (feedback phase-modulation) and instantaneous envelope:
    double phaseModByEnv = -100.0;
    sinCos((PI/180.0)*startPhase + phaseMod*z.im + phaseModByEnv*(z.getRadius()-1.0), &cr, &ci);

    /*
    // phase modulation by :
    double env = z.getRadius();
    double cr2, ci2;
    sinCos((PI/180.0)*startPhase + phaseModByEnv*z.getRadius(), &cr2, &ci2);
    return env;
    */

    return amplitude * (cr*z.re + ci*z.im);
  }

  INLINE double ModalFilterWithAttack::getSample(double in)
  {
    return modalFilter1.getSample(in) - modalFilter2.getSample(in);
  }

  INLINE double ModalSynthesizer::getSample(double in)
  {
    double accu = 0.0;
    int nm      = rmin(numModes, frequencies.dim, amplitudes.dim, decayTimes.dim);
    nm          = rmin(nm, startPhases.dim);
    for(int m=0; m<nm; m++)
      accu += modalFilters[m].getSample(in);
    return accu;
  }

} // end namespace rosic

#endif // rosic_ModalSynthesizer_h
