#pragma once
namespace rosic
{

/**  */

class rsModalSynth
{

public:

  /** A selection of predefined frequency ratio profiles. */
  enum freqRatioProfiles
  {
    HARMONIC = 0,
    STIFF_STRING,
    //POWER_RULE,
    //RECTANGULAR_MEMBRANE,
    //CIRCULAR_MEMBRANE,

    NUM_FREQ_RATIO_PROFILES
  };

  /** \name Setup */

  /** Selects one of the predefined frequency ratio profiles. For example, a "harmonic" prfile means 
  that the ratios should be 1,2,3,4,5, etc. i.e. the n-th partial has a frequency ratio (with 
  respect to the fundamental) of n. The "stiff-string" setting uses the formula  n * sqrt(1+B*n^2) 
  where B is the "inharmonicity" parameter, see Eq 10, here:
  http://www.simonhendry.co.uk/wp/wp-content/uploads/2012/08/inharmonicity.pdf
  http://www.jbsand.dk/div/StivStreng.pdf
  ...more profiles are to come, including circular membranes and totally made up fantasy formulas */
  void setFreqRatioProfile1(int newProfile);
  // todo: let the user define their own profiles and load them from an xml and/or define a formula

  void setFreqRatioProfile2(int newProfile);

  /** The raw frequency ratios, as determined by the selected frequency ratio profile, may be 
  re-adjusted, for example to coincide with the intervals of 12-tone equal temperament. This 
  parameter sets up, how much they should be re-adjusted (range: 0..1). */
  //void setFreqRatioAdjustment(double newAdjustmentAmount);
  // maybe we should just have two freq-ratio profiles on the same footing and the the user
  // interpolate between the two - or maybe even 4 and use a vector-pad bilinear (or bi-log-linear) 
  // interpolation scheme


  void setFreqRatioMix(double newMix);

  /** Sets the inharmonicity parameter for the stiff string frequency ratio profile. */
  void setInharmonicity(double newInharmonicity);


  /** Sets a limit for the number of partials - use this to keep CPU load under control. */
  //void setMaxNumPartials(int newMax);


  /** \name Processing */

  inline void getSampleFrameStereo(double* outL, double* outR)
  {
    float x = getExcitation(); // later maybe pass the input signal to the exciter
    rsFloat32x4 y = modalBank.getSample( rsFloat32x4(x) );
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


protected:

  /** Fills the "ratios" array with the frequency ratios of given "profile". */
  void fillFreqRatios(double* ratios, int profile);

  void fillFreqRatiosHarmonic(double* ratios);
  void fillFreqRatiosStiffString(double* ratios, double B);

  void updateFreqRatios();

  int freqRatioProfile1 = HARMONIC;
  int freqRatioProfile2 = HARMONIC;
  double freqRatioMix   = 0.0;
  double inharmonicity  = 0.0;

  static const int maxNumModes = rsModalBankFloatSSE2::maxNumModes; // for convenience
  rsModalBankFloatSSE2 modalBank;

  // move to subclass rsModalSynthPoly
  //static const int maxNumVoices = 16;
  //rsModalBankFloatSSE2 modalBanks[maxNumVoices];

  double freqRatios[maxNumModes];   // the final mix/interpolation of freq ratios
  double freqRatios1[maxNumModes];
  double freqRatios2[maxNumModes];
  bool logLinearFreqInterpolation = false;

  int noteAge = 0;

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





}