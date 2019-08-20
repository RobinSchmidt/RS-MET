#pragma once
namespace rosic
{

/** A class for generating the frequency ratios of the modes for modal synthesis. In the very
simplest case, the frequency ratios are just all he integers, giving a harmonic series. There are
other formulas, some of which are based on musical considerations, some on physical ones and some
ad-hoc. For the physically inspired mode frequency ratios, such as those of a rod with two free
ends, note that these values are computed from a theoretically idealized situation. In a practical
instrument, the mode frequencies may be shifted away from these idealized ones due to non-ideal or
even purposefully changed shapes - and there seems to be no real musical significance of the
*exact* values computed by some of the physically based formulas which are based on simple,
idealized geometries. In fact, it's the job of a good instrument maker to adjust physical
parameters such as shapes in order tweak the mode frequencies from the ratios obtained from simple
geometries into more musically meaningful ones - so it seems to be a rather pointless excersise to
exactly hit these exact frequency ratios of simply shaped resonators - however, if they are given
in a book, we may want to have code for that. They may serve as a starting point for finding more
musically interesting ratios that preserve some general qualities of the underlying instrument
classes. But the true power of modal synthesis actually lies in the fact that we have *direct*
access to the modal frequencies and may tune them individually to musical taste rather than
indirectly tweaking them in complicated interrelated ways by altering a (possibly complicated)
geometry of a physical resonator. */

class rsModalFrequencyGenerator
{

public:

  /** Fills the array "r" of length "N" with the frequency ratios corresponding to a harmonic
  series. */
  static void allHarmonics(double* r, int N);

  static void oddHarmonics(double* r, int N);

  /** The mode frequencies are powers of the twelfth-root-of-two, but for the lower modes, not all
  powers are used but only those that are close to harmonic ratios. This mode tuning will play
  maximally consonant in 12-TET in the sense explained here:
  http://sethares.engr.wisc.edu/consemi.html */
  static void twelveTone(double* r, int N);

  /** The lower modes are the same as in twelveTone but the upper modes follow a harmonic series
  instead of being yet more powers of the twelfth-root-of-two. */
  static void twelveTonePseudoHarmonic(double* r, int N);

  /** Calculates the relative frequencies of the vibrational modes of a stiff rod with both ends
  free. */
  static void rodFreeFree(double* r, int N);

  /** Calculates the relative frequencies of the vibrational modes of a stiff rod with one free end
  and one clamped end. */
  static void rodFreeClamped(double* r, int N);


  //static void idealBar(double* r, int N);

  //static void circularMembrane(double* r, int N);

  /** Frequency ratios of a stiff string, such as a piano string. the parameter B controls the
  amount of stiffness and hence the amount of inharmonicity. */
  static void stiffString(double* r, int N, double B);

protected:

  /** First 21 modes for the 12-TET based tunings. Used by twelveTone and
  twelveTonePseudoHarmonic */
  static void twelveTone21(double* r, int N);


};

// not yet used, but may (or may not) be useful later
//struct rsModalUserParameters
//{
//  double frequency   = 440;   // in Hz
//  double amplitude   = 1.0;   // as raw factor ...or should we use decibels?
//  double phase       = 0.0;   // in degrees
//  double attack      = 50;    // in milliseconds
//  double decay       = 500;   // in milliseconds
//  double freqSpread  = 0;     // in Hz
//  double phaseDelta  = 0;
//  double blend       = 0.5;
//  double attackScale = 1.0;
//  double decayScale  = 1.0;
//};
//
//class rsModalAlgoParameters
//{
//
//public:
//
//  /** Sets our members by converting a set of user parameters for a modal filter to the
//  corresponding set of algorithm parameters. */
//  void setFromUserParameters(const rsModalUserParameters& userParams, double sampleRate);
//
//protected:
//
//  double w;      // 2*pi*f/fs
//  double A;      // raw amplitude factor
//  double p;      // phase in radians
//  double att;    // attack in samples
//  double dec;    // decay in samples
//  double dw;     // delta omega
//  double dp;     // phase delta in radians
//  double b;      // blend
//  double attScl; // attack scale
//  double decScl; // decay scale
//};

//=================================================================================================

/** A monophonic modal synthesizer. The parameters for the modes (their relative frequencies,
amplitudes, decay times, attack times, etc.) are set up in terms of higher level macro parameters.
When you trigger a note, an impulse is fed inot the bank of modal filters and you hear the impulse
response

ToDo: let it process incoming audio, add noise to the excitation signal while the note is on */

class rsModalSynth
{

public:


  rsModalSynth();


  /** A selection of predefined frequency ratio profiles. */
  enum freqRatioProfiles
  {
    // 0-parametric:
    ALL_HARMONICS = 0,
    ODD_HARMONICS,
    TWELVE_TONE_PSEUDO_HARMONIC,
    TWELVE_TONE_EQUAL,
    ROD_FREE_FREE,
    ROD_FREE_CLAMPED,
    //IDEAL_BAR,
    //CIRCULAR_MEMBRANE,

    // 1-parametric:
    //STIFF_STRING,            // stiffness/inharmonicity
    //POWER_RULE,              // power/exponent
    //RECTANGULAR_MEMBRANE,    // aspect ratio (0..1)
    //ELLIPTIC_MEMBRANE        // dito

    //CUSTOM_FORMULA,
    //CUSTOM_DATA,

    NUM_FREQ_RATIO_PROFILES
  };


  /** \name Setup */

  void setSampleRate(double newRate) { sampleRate = newRate; }

  /** Selects one of the predefined frequency ratio profiles for the top-left slot. For example, a
  "harmonic" profile means that the ratios should be 1,2,3,4,5, etc. i.e. the n-th partial has a
  frequency ratio (with respect to the fundamental) of n. The "stiff-string" setting uses the
  formula  n * sqrt(1+B*n^2) where B is the "inharmonicity" parameter, etc. */
  void setFreqRatioProfileTopLeft(int newProfile);
  // todo: let the user define their own profiles and load them from an xml and/or define a formula

  void setFreqRatioProfileTopRight(int newProfile);
  void setFreqRatioProfileBottomLeft(int newProfile);
  void setFreqRatioProfileBottomRight(int newProfile);

  /** We use a vector mix/morph between 4 frequency ratio profiles. This sets the x-coordinate of the mixing
  vector. */
  void setFreqRatioMixX(double newMix);

  /** @see setFreqRatioMixX - same for the y-coordinate. */
  void setFreqRatioMixY(double newMix);

  /** Sets the inharmonicity parameter for the stiff string frequency ratio profile. */
  void setInharmonicity(double newInharmonicity) { inharmonicity = newInharmonicity; }

  /** Sets a limit for the number of partials - use this to keep CPU load under control. */
  //void setMaxNumPartials(int newMax) { numPartialsLimit = newMax; }
    // rename to setMaxPartialIndex, add setMinPartialIndex -> brickwall highpass


  void setLowestMode(int newLimit) { lowestMode = newLimit; }

  void setHighestMode(int newLimit) { highestMode = newLimit; }



  //void setSpectralSlope(double newSlope) { spectralSlope = newSlope; }
  //void setSpectralSlopeByKey(double newSlopeByKey) { spectralSlopeByKey = newSlopeByKey; }
  //void setSpectralSlopeByVel(double newSlopeByKey) { spectralSlopeByVel = newSlopeByKey; }

  // ByRatio, ByKey, ByVel parameters should be passed as percentage values

  /** Sets the overall level in decibels.*/
  void setLevel(       double newLevel)        { level      = newLevel; }

  /** Sets the key tracking of overall level. The unit is +-dB at extreme keys with respect to the
  neutral reference key 64. When you set LevelByKey to 10 dB, the level will be unchanged, if the
  key is 64, 10 dB louder when the key is 127 and 10 dB quieter when the key is 1. */
  void setLevelByKey(  double newLevelByKey)   { levelByKey = newLevelByKey; }

  /** Velocity tracking of overall level. Unit is +-dB at extreme velocities */
  void setLevelByVel(  double newLevelByVel)   { levelByVel = newLevelByVel; }

  // insert Tune/Key/Vel

  /** The slope of the magnitude spectrum in dB/oct. Actually, due to the different attack and
  decay-rates, this slope changes over time - this here sets the reference slope that applies
  when all attack and decay times are equal for each mode. */
  void setAmpSlope(     double newAmpSlope)      { ampSlope      = newAmpSlope; }

  /** Key tracking of the slope of magnitude spectrum. The unit is +-dB/oct at extreme keys. */
  void setAmpSlopeByKey(double newAmpSlopeByKey) { ampSlopeByKey = newAmpSlopeByKey; }

  /** Velocity tracking of slope of magnitude spectrum. */
  void setAmpSlopeByVel(double newAmpSlopeByVel) { ampSlopeByVel = newAmpSlopeByVel; }

  /** Reference attack time in milliseconds. */
  void setAttack(       double newAttack)        { attack        = newAttack; }

  /** Dependency of the attack time on the modal frequency ratio in percent. At 100%, it means that
  a mode that has twice the frequency of the fundamental will also have twice the attack time, at
  -100% it will have half the attack time (typically you will want to use negative values - higher
  frequency modes develope faster in time). */
  void setAttackByRatio(double newAttackByRatio) { attackByRatio = 0.01*newAttackByRatio; }

  /** Key tracking of the attack time in percent. 100% means, that a note that is an octave above
  the reference key 64 will get twice the attack time. Similarly to setAttackByRatio, you will
  typically want negative values. */
  void setAttackByKey(  double newAttackByKey)   { attackByKey   = 0.01*newAttackByKey; }

  /** Velocity tracking of attack time in percent. 100% means... */
  void setAttackByVel(  double newAttackByVel)   { attackByVel   = 0.01*newAttackByVel; }


  // decay parameters are entirely analogous to attack parameters:
  void setDecay(       double newDecay)        { decay        = newDecay; }
  void setDecayByRatio(double newDecayByRatio) { decayByRatio = 0.01*newDecayByRatio; }
  void setDecayByKey(  double newDecayByKey)   { decayByKey   = 0.01*newDecayByKey; }
  void setDecayByVel(  double newDecayByVel)   { decayByVel   = 0.01*newDecayByVel; }


  void setPhaseRandomness(double newRandomness) { phaseRandomness = newRandomness; }
  void setPhaseRandomSeed(int newSeed) { phaseRandomSeed = newSeed; }



  /** \name Inquiry */


  //int getMaxNumPartials() const { return numPartialsLimit; }

  int getLowestMode() const { return lowestMode; }

  int getHighestMode() const { return highestMode; }

  // these functions are mostly for reference and doing plots to verify the expected behavior - for
  // realtime purposes, optimized computations can be used

  // todo: let the key/vel parameters also be double

  /** Frequency ratio of given mode at given key and velocity */
  double getModeFreqRatio(int modeIndex, double key, double vel) const;

  /** Absolute frequency in Hz of a given mode at given key and velocity. */
  double getModeFrequency(int modeIndex, double key, double vel) const;

  /** Pitch offset of given mode with respect to the fundamental in semitones. */
  double getModeRelativePitch(int modeIndex, double key, double vel) const;

  /** Amplitude of given mode as raw factor. */
  double getModeAmplitude(int modeIndex, double key, double vel) const;

  /** Level of given mode in decibels. */
  double getModeLevel(int modeIndex, double key, double vel) const;

  /** Attack time of given mode in milliseconds. */
  double getModeAttack(int modeIndex, double key, double vel) const;

  /** Decay time of given mode in milliseconds. */
  double getModeDecay(int modeIndex, double key, double vel) const;

  // to compute the total decay time, we use the formula
  // decay = baseDecay * exp(cr*decayByRatio  + ck*decayByKey  + cv*decayByVel);
  // these 3 functions compute the coefficients cr, ck, cv (maybe make protected):
  double getTimeCoeffByIndex(int index) const;
  double getTimeCoeffByKey(double key) const;
  double getTimeCoeffByVel(double vel) const;


  //int getMaxNumModes() const { return maxNumModes; }








  /** \name Processing */

  inline void getSampleFrameStereo(double* outL, double* outR)
  {
    //if(noteAge > decayLength)
    //  return;
    // todo: we should actually somehow reset the note-length whenever we receive a new excitation
    // but then "noteAge" is the wrong name...but maybe we don't need this cutting off of notes
    // anyway - when we use continuous excitation signals, it doesn't make sense anymore anyway


    float x = getExcitation(); // later maybe pass the input signal to the exciter
    rsFloat32x4 y = modalBank.getSample(rsFloat32x4(x));
    *outL = *outR = y.getSum();  // preliminary - later do stereo mixing
    noteAge++;
  }

  inline float getExcitation()
  {
    // preliminary - just a unit impulse at the start of the note
    if(noteAge == 0)
      return 1.f;
    return 0.f;

    // todo: scale impulse by user adjustable factor and add noise while note is on (also with
    // use adjustable factor
  }

  /** A note-on function where key and vel are passed as floating point numbers because in the
  context of experimenting, it may be convenient to allow arbitrary/continuous frequencies. */
  void noteOn(double key, double velocity);

  /** A regular note-on function, accepting integers for key and vel. */
  void noteOn(int key, int velocity) { noteOn(double(key), double(velocity)); }

  void reset()
  {
    modalBank.reset();
    //noteAge = 0;
  }


protected:

  /** Fills the "ratios" array with the frequency ratios of given "profile" and the logRatios array
  with their logarithms. */
  void fillFreqRatios(double* ratios, double *logRatios, int profile);

  void updateFreqRatios();

  static const int maxNumModes = rsModalBankFloatSSE2::maxNumModes; // for convenience

  // global settings:
  double level = 0, levelByKey = 0, levelByVel = 0;

  // modal frequency settings:
  //int numPartialsLimit  = maxNumModes;  // obsolete

  int lowestMode  = 1;
  int highestMode = maxNumModes;

  int freqRatioProfileTopLeft     = ALL_HARMONICS;
  int freqRatioProfileTopRight    = ALL_HARMONICS;
  int freqRatioProfileBottomLeft  = ALL_HARMONICS;
  int freqRatioProfileBottomRight = ALL_HARMONICS;
  double freqRatioMixX  = 0.0;
  double freqRatioMixY  = 0.0;
  double inharmonicity  = 0.0;
  //bool interpolatePitches = true;  // alway do it like this

  // macro parameters for modal filters (maybe wrap into a class):
  //double amplitude = 1.0; // maybe have a level parameter in dB with key and vel scaling as in Straightliner
  //double spectralSlope = 0, spectralSlopeByKey = 0, spectralSlopeByVel = 0;
  //double amp    = 1.0, ampByRatio    = -1.0, ampByKey    = 0, ampByVel    = 0; // get rid of ampByRatio
  double ampSlope = 0, ampSlopeByKey = 0, ampSlopeByVel = 0;
  double attack = 5,   attackByRatio = -1.0, attackByKey = 0, attackByVel = 0;
  double decay  = 500, decayByRatio  = -0.3, decayByKey  = 0, decayByVel  = 0;

  double phaseRandomness = 1.0;
  int phaseRandomSeed = 0;
  int phaseRandomShape = 1; // determines the probability distribution ...later
  // double evenScale, evenScaleByKey, evenScaleByVel
  // double decayCombFreq
  //...more parameters to come....

  rsModalBankFloatSSE2 modalBank;

  RAPT::rsNoiseGenerator<double> phaseGenerator;

  // move to subclass rsModalSynthPoly
  //static const int maxNumVoices = 16;
  //rsModalBankFloatSSE2 modalBanks[maxNumVoices];

  // frequency ratio arrays and their logarithmic versions:
  double freqRatios[maxNumModes];     // the final mix/interpolation of freq ratios
  double freqRatiosLog[maxNumModes];  // logarithms of final freq ratios
  double freqRatiosTopLeft[maxNumModes];
  double freqRatiosTopRight[maxNumModes];
  double freqRatiosBottomLeft[maxNumModes];
  double freqRatiosBottomRight[maxNumModes];
  double freqRatiosLogTopLeft[maxNumModes];
  double freqRatiosLogTopRight[maxNumModes];
  double freqRatiosLogBottomLeft[maxNumModes];
  double freqRatiosLogBottomRight[maxNumModes];
  // maybe keep only the log-ratios...we'll see

  int noteAge = 0;
  int decayLength = INT_MAX;
  double sampleRate = 44100;

  std::atomic_bool freqRatiosAreReady;
  //std::atomic_bool freqRatiosAreReady = false; // doesn't compile on gcc
  //std::atomic_bool freqRatiosAreReady(false);
};

// let the user define the modal parameters at various keys and for each key that is defined, have
// a high and a low velocity setting. on note-on, these datapoints are interpolated for the current
// note and velocity setting. if the incoming note is outside the range of defined key, use linear
// extrapolation
// maybe the user should not specify the low-level per mode parameters but instead macro parameters
// like inharmonicity, attack, decay, attackByFreq, decayByFreq, etc. - perhaps have selectors
// for harmonic-series, stiff-string series, bar-series, 12-TET series, etc and have also an option
// to enter formulas


/** A polyphonic version of the modal synthesizer. */

class rsModalSynthPoly : public rsModalSynth // public rsPolyphonicInstrument
{

public:

protected:

  static const int maxNumVoices = 16; // factor out into rsPolyphonicInstrument
  rsModalBankFloatSSE2 modalBanks[maxNumVoices];

};



}
