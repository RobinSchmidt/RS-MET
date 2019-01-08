#pragma once
namespace rosic
{

/** A class for generating the frequency ratios of the modes for modal synthesis. */

class rsModalFrequencyGenerator
{

public:

  /** Fills the array "r" of length "N" with the frequency ratios corresponding to a harmonic 
  series. */
  static void harmonic(double* r, int N);

  /** The mode frequencies are powers of the twelfth-root-of-two, but for the lower modes, not all
  powers are used but only those that are close to harmonic ratios. This mode tuning will play 
  maximally consonant in 12-TET in the sense explained here:
  http://sethares.engr.wisc.edu/consemi.html */
  static void twelveTone(double* r, int N);

  /** The lower modes are the same as in twelveTone but the upper modes follow a harmonic series 
  instead of being yet more powers of the twelfth-root-of-two. */
  static void twelveTonePseudoHarmonic(double* r, int N);

  /** Calculates the relative frequencies of the vibrational modes of a stiff rod with both ends 
  free. Reference: discussion on the Xylophone in Dave Benson's 'Mathematics and Music' */
  static void rodFreeFree(double* r, int N);

  /** Calculates the relative frequencies of the vibrational modes of a stiff rod with one free end 
  and one clamped end. Reference: discussion on the Mbira in Dave Benson's 'Mathematics and 
  Music' */
  static void rodFreeClamped(double* r, int N);


  //static void idealBar(double* r, int N);

  //static void circularMembrane(double* r, int N);


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

/** A monophonic modal synthesizer.... */

class rsModalSynth
{

public:

  /** A selection of predefined frequency ratio profiles. */
  enum freqRatioProfiles
  {
    // 0-parametric:
    HARMONIC = 0,
    TWELVE_TONE_EQUAL,
    IDEAL_BAR,
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
  formula  n * sqrt(1+B*n^2) where B is the "inharmonicity" parameter, see Eq 10, here:
  http://www.simonhendry.co.uk/wp/wp-content/uploads/2012/08/inharmonicity.pdf
  http://www.jbsand.dk/div/StivStreng.pdf
  ...more profiles are to come, including circular membranes and totally made up fantasy formulas */
  void setFreqRatioProfile1(int newProfile);
  // todo: let the user define their own profiles and load them from an xml and/or define a formula

  void setFreqRatioProfile2(int newProfile);
  void setFreqRatioProfile3(int newProfile);
  void setFreqRatioProfile4(int newProfile);

  /** We use a vector mix/morph between 4 frequency ratio profiles. This sets the x-coordinate of the mixing
  vector. */
  void setFreqRatioMixX(double newMix);

  /** @see setFreqRatioMixX - same for the y-coordinate. */
  void setFreqRatioMixY(double newMix);

  /** Sets the inharmonicity parameter for the stiff string frequency ratio profile. */
  void setInharmonicity(double newInharmonicity) { inharmonicity = newInharmonicity; }

  /** Sets a limit for the number of partials - use this to keep CPU load under control. */
  void setMaxNumPartials(int newMax) { numPartialsLimit = newMax; }


  void setSpectralSlope(double newSlope) { spectralSlope = newSlope; }
  void setAttack(double newAttack) { attack = newAttack; }
  void setDecay(double newDecay) { decay = newDecay; }
  void setPhaseRandomness(double newRandomness) { phaseRandomness = newRandomness; }
  void setPhaseRandomSeed(int newSeed) { phaseRandomSeed = newSeed; }

  void setSpectralSlopeByKey(double newSlopeByKey) { spectralSlopeByKey = newSlopeByKey; }
  void setAttackByKey(double newAttackByKey) { attackByKey = newAttackByKey; }
  void setDecayByKey(double newDecayByKey) { decayByKey = newDecayByKey; }

  void setSpectralSlopeByVel(double newSlopeByKey) { spectralSlopeByVel = newSlopeByKey; }
  void setAttackByVel(double newAttackByVel) { attackByVel = newAttackByVel; }
  void setDecayByVel(double newDecayByVel) { decayByVel = newDecayByVel; }


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
  }

  void noteOn(int key, int velocity);

  void reset()
  {
    modalBank.reset();
    //noteAge = 0;
  }


protected:

  /** Fills the "ratios" array with the frequency ratios of given "profile" and the logRatios array
  with their logarithms. */
  void fillFreqRatios(double* ratios, double *logRatios, int profile);


  //void fillFreqRatiosHarmonic(double* ratios);
  //void fillFreqRatiosStiffString(double* ratios, double B);
  void updateFreqRatios();


  static const int maxNumModes = rsModalBankFloatSSE2::maxNumModes; // for convenience

  // modal frequency settings:
  int numPartialsLimit  = maxNumModes;
  int freqRatioProfile1 = HARMONIC;  // top-left
  int freqRatioProfile2 = HARMONIC;  // top-right
  int freqRatioProfile3 = HARMONIC;  // bottom-left
  int freqRatioProfile4 = HARMONIC;  // bottom-right
  double freqRatioMixX  = 0.5;
  double freqRatioMixY  = 0.5;
  double inharmonicity  = 0.0;


  // macro parameters for modal filters:
  double amplitude = 1.0; // maybe have a level parameter in dB with key and vel scaling as in Straightliner
  double spectralSlope = 0, spectralSlopeByKey = 0, spectralSlopeByVel = 0;
  double attack = 50, attackByKey = 0, attackByVel = 0;
  double decay = 500, decayByKey = 0, decayByVel = 0;
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
  double freqRatios[maxNumModes];   // the final mix/interpolation of freq ratios
  double freqRatios1[maxNumModes];
  double freqRatios2[maxNumModes];
  double freqRatios3[maxNumModes];
  double freqRatios4[maxNumModes];
  double freqRatiosLog1[maxNumModes];
  double freqRatiosLog2[maxNumModes];
  double freqRatiosLog3[maxNumModes];
  double freqRatiosLog4[maxNumModes];
  bool logLinearFreqInterpolation = false;


  int noteAge = 0;
  int decayLength = INT_MAX;
  double sampleRate = 44100;

  std::atomic_bool freqRatiosAreReady = false;
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