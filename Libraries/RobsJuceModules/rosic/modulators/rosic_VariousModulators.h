#ifndef rosic_AttackDecayEnvelope_h
#define rosic_AttackDecayEnvelope_h
// rename to MiscModulators

namespace rosic
{


/** This is a class implements an envelope generator based on feeding the output signal of a
multiplicative accumulator (which is decaying exponential) into a leaky integrator (which smoothes 
out the (otherwise instantaneous) attack). The output of the envelope is normalized to the range 
0...1.  */

class AttackDecayEnvelope
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AttackDecayEnvelope();

  /** Destructor. */
  ~AttackDecayEnvelope();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the time instant at which the peak excursion occurs (in milliseconds). */
  void setPeakTime(double newPeakTime);

  /** Sets the time constant for the multiplicative accumulator (which we consider as primarily
  responsible for the decaying part) in milliseconds. */
  void setDecayTimeConstant(double newTimeConstant);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the time instant at which the peak excursion occurs (in milliseconds). */
  double getPeakTime() const { return 1000.0*tp; }

  /** Returns the length of the decay phase (in milliseconds). */
  double getDecayTimeConstant() const { return 1000.0*td; }

  /** True, if output is below some threshold. */
  bool endIsReached(double threshold);

  //-----------------------------------------------------------------------------------------------
  // \name Processing:

  /** Calculates one output sample at a time. */
  INLINE double getSample();

  void trigger(bool startFromCurrentLevel);

  /** Resets the time variable. */
  void reset();

protected:

  /** Calculates the attack time constant ta from the desired time-instant of the peak excursion
  tp and the decay time constant td. */
  void calculateAttackTimeConstant();

  /** Calculates the filter coefficients and the in initial value for the internal state of the
  decay-filter (this realizes the amplitude normalization). */
  void calculateCoeffsAndInitValue();

  double ca, cd;     // coefficients for attack and decay
  double ya, yd;     // previous outputs of attack and decay filter
  double ydInit;     // initial value for previous output of decay filter
  double ta, td;     // attack and decay time-constants (in seconds)
  double tp;         // time-instant of the peak excursion (in seconds);
  double fs;         // sample-rate
  int    n, np;      // current sample-index and sample-index of the peak

};

INLINE double AttackDecayEnvelope::getSample()
{
  yd *= cd;
  ya  = yd + ca*(ya-yd);
  return ya;
}

//=================================================================================================

/** An envelope generator class based on RAPT::rsTriSawOscillator. Instead of setting it up via
a frequency and an asymmetry parameter, this class is parametrized in terms of attack- and decay 
times. */

class rsTriSawModulator : public RAPT::rsTriSawOscillator<double>
{
  typedef RAPT::rsTriSawOscillator<double> Base; // for conveniently calling baseclass methods

public:

  /** Sets the sample rate in Hz. */
  void setSampleRate(double newSampleRate);

  /** Sets the attack time in milliseconds. */
  void setAttackTime(double newAttack);

  /** Sets the decay time in milliseconds. */
  void setDecayTime(double newDecay);

  /** Scales the overall length of the envelope with given factor. */
  void setTimeScaler(double newScaler);

  /** Sets the floor, i.e. the lowest value. */
  void setFloor(double newFloor) 
  {
    //RAPT::rsAssert(RAPT::rsIsFiniteNumber(newFloor));
    flooor = newFloor; 
  }

  /** Sets the ceiling, i.e. the highest value. */
  void setCeiling(double newCeiling) 
  { 
    //RAPT::rsAssert(RAPT::rsIsFiniteNumber(newCeiling));
    ceiling = newCeiling; 
  }

  /** Produces one output sample at a time. */
  inline double getSample() /* override ...but baseclass method is not virtual */
  {
    return flooor + (ceiling-flooor) * 0.5 * (Base::getSample() + 1);
  }

protected:

  /** Called from setAttackTime, setDecayTime, etc. to update the parameters of the underlying 
  oscillator. */
  void updateOscParameters();

  double sampleRate = 44100, attack = 0.01, decay = 0.99, timeScale = 1;
  double flooor = -1, ceiling = +1;

  // overriden because they should not be used in this subclass:
  inline void setPhaseIncrement(double newInc) { Base::setPhaseIncrement(newInc); }
  inline void setAsymmetry(double newAsym)     { Base::setAsymmetry(newAsym); }
  // use setAttackTime, setDecayTime instead

};


} // end namespace rosic

#endif // rosic_AttackDecayEnvelope_h
