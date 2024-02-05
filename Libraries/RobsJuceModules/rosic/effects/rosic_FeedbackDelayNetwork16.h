#ifndef rosic_FeedbackDelayNetwork16_h
#define rosic_FeedbackDelayNetwork16_h

namespace rosic
{

/** This class implements a feedback-delay-network with 16 delay-lines which is supposed to be used
mainly as a reverberation effect unit. Each of the 16 delaylines has a damping-filter in series
which serves as a frequency dependent feedback gain. The feedback itself is realized via various
pre-defined feedback-matrices, which are implemented by efficient techniques (for, example a
multiplication by a Hadamard-matrix is done via the Fast Hadamard Transform). The injection- and
output-vectors can also be chosen from a list of pre-defined vectors. The output-signals of the
delay-lines are fed through compensation-filters before actually being routed to the output. These
filters decouple the overall reverberation-time from the output volume - and they even do this in
a frequency dependent way. Roughly speaking, they act opposite to the damping-filters in order to
decoulour (whiten) the wet signal in a steady-state condition.

idea for modulation:
 -introduce a 1st order allpass after each delayline with an allpass coefficient that depends on
  the cross-correlation of the (wet) output signal
  maybe and/or enforce a certain cross-correlation by M/S-processing (apply time varying gain to M
  and S according to cross-correlation of input so as to
 -maybe pass the output signal through narrow bandpass filter, perhaps pass this bandpass signal 
  through a leveller and use it as modulation signal for spread-parameter in the generalized
  hadamard transform - gives a modulation that is derived from the input signal itself but still
  has a user-controllable frequency (via the bandpass frequency) */

class FeedbackDelayNetwork16
{

public:

//#define ONE_OVER_SQRT8 0.35355339059327376220042218105242
#define ONE_OVER_SQRT16 0.25

    /** This is an enumeration of the available injection vectors */
  enum injectionVectors
  {
    IN_ALL_ONES = 0,
    IN_ALTERNATING_PLUSMINUS_01,
    IN_ALTERNATING_PLUSMINUS_02,

    NUM_INJECTION_VECTORS
  };

  /** This is an enumeration of the available feedback matrices */
  enum feedbackMatrices
  {
    IDENTITY = 0,
    MINUS_IDENTITY,
    SERIES_CONNECTION,
    SERIES_WITH_GLOBAL_FEEDBACK,
    HADAMARD,
    MINUS_HADAMARD,

    NUM_FEEDBACK_MATRICES
  };

  enum outputVectors
  {
    OUT_ALL_ONES = 0,
    OUT_01,
    OUT_02,
    OUT_03,
    OUT_04,
    OUT_05,
    //OUT_06,

    NUM_OUTPUT_VECTORS
  };

  enum delayDistributions
  {
    LINEAR = 0,
    DISTANCE_DECAY,
    SIMPLE_RATIO_MODES,
    GEOMETRIC_MEANS,


    EXPONENTIAL,
    DIVISION_ALGO_1,
    PRIME_ALGO_1,
    PRIME_ALGO_2,

    NUM_DELAY_DISTRIBUTIONS
  };

  enum delayOrderings
  {
    ASCENDING = 0,
    DESCENDING,
    ALTERNATING,

    NUM_DELAY_ORDERINGS
  };

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  FeedbackDelayNetwork16();

  /** Destructor. */
  ~FeedbackDelayNetwork16();

  //---------------------------------------------------------------------------------------------
  //parameter settings:

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Scales the reverb time for the low-band respect to the mid-band. */
  void setLowReverbTimeScale(double newLowReverbTimeScale);

  /** Sets the crossover-frequency between low and mid frequencies. */
  void setLowCrossoverFreq(double newLowCrossoverFreq);

  /** sets the time in which the amplitude decays to -60 dB for the mid frequencies. */
  void setMidReverbTime(double newMidReverbTime);

  /** Scales the reverb time for the high-band respect to the mid-band. */
  void setHighReverbTimeScale(double newHighReverbTimeScale);

  /** Sets the crossover-frequency between mid and high frequencies. */
  void setHighCrossoverFreq(double newHighCrossoverFreq);

  /** Selects, if the delay-lines should operate as comb- or as allpass-filter.
  TODO: allow also for splitting the single allpass into a series connection of two (or more
  allpasses - this will dramatically increase the time density) ... maybe pass an integer
  parameter to determine the number of allpasses
  (between 0 (comb-mode) and some maximum number). */
  void setAllpassMode(bool shouldBeAllpass);

  /** Sets the volume of the dry signal in dB - TODO: use setDryWetRation instead */
  void setDryVolume(double newDryVolume);

  /** Sets the volume of the wet signal in dB - TODO: use setDryWetRation instead */
  void setWetVolume(double newWetVolume);

  /** Sets the cutoff frequency of a lowpass filter for the wet signal. */
  void setWetLpfCutoff(double newWetLpfCutoff);

  /** Sets the cutoff frequency of a highpass filter for the wet signal. */
  void setWetHpfCutoff(double newWetHpfCutoff);

  /** Switches swapping of left and right output channel on and off. */
  void setStereoSwapSwitch(bool newStereoSwapSwitch);

  /** Switches a pinking filter for the wet signal on and off. */
  void setWetPinkingSwitch(bool newWetPinkingSwitch);

  /** Sets a reference delay time in milliseconds - the actual delay-times (in milliseconds) of
  the delaylines may then be specified relative to this value via setRelativeDelayTime. */
  void setReferenceDelayTime(double newReferenceDelayTime);

  /** Sets the relative delaytime for one of the delaylines as factor with respect to the
  'referenceDelayTime' member. */
  void setRelativeDelayTime(int delayLineIndex, double newRelativeDelayTime);

  /** Sets the relative delaytimes for all delaylines at once - the passed pointer should point
  to an array containing 'maxNumDelayLines' (== 16) entries). */
  void setAllRelativeDelayTimes(double *newRelativeDelayTimes);

  /** Assigns the relative delaytimes by some algorithm (see delayDistributions), possibly taking
  some control parameters, the meaning of which will vary depending on the chosen
  algorithm/distribution. */
  void assignRelativeDelayTimesAlgorithmically(int distributionIndex, double parameter1,
    double parameter2);

  /** Chooses, how the delay-lines of different length are ordered with respect to their length
  (ascending, descending, alternating, etc. MAY BE DEPRECATED */
  void setDelayOrdering(int newDelayOrdering);

  /** Chooses one of the pre-defined injection-vectors. */
  void setInjectionVector(int newInjectionVectorIndex);

  /** Chooses one of the pre-defined feedback-matrices. */
  void setFeedbackMatrix(int newFeedbackMatrixIndex);

  /** Chooses one of the pre-defined output-vectors. */
  void setOutputVector(int newOutputVectorIndex);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the relative delaytime of one of the delaylines (0 if index is out of range). */
  double getRelativeDelayTime(int index) const
  {
    if(index >= 0 && index < maxNumDelayLines)
      return relativeDelayTimes[index];
    else
      return 0.0;
  }

  /** Writes the impulse response into the passed buffers (for left and right channel), where the
  buffers should have (at least) the passed length. */
  void getImpulseResponse(double *bufferL, double *bufferR, int length);

  //---------------------------------------------------------------------------------------------
  // audio processing:

  /** Feeds in an impulse to the  reverb unit for auditioning the impulse response. */
  void feedInImpulse(float amplitude = 1.f) { impulseAmplitude = amplitude; }

  INLINE void getSampleFrameStereo(double* inOutL, double* inOutR);
  // get rid, use processFrame instead

  INLINE void processFrame(double* inOutL, double* inOutR)
  { getSampleFrameStereo(inOutL, inOutR); }
  // alias to getSampleFrameStereo to match API of FeedbackDelayNetwork

  //---------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers to zero. */
  void reset();

  //=============================================================================================

protected:

  /** Triggers a re-adjustment of one of the read-pointers - this is necesarry when the
  delay-time of the corresponding delay-line changes. It is also necessarry to do this for all
  delay-lines when the sample-rate changes. */
  void adjustReadPointer(int index);

  /** Adjusts the lengths of the delay-lines according to minDelayTime, maxDelayTime,
  delayDistribution and delayDistributionForm.
  TODO: include the new distribution based on exponentially decaying echo-distances
  */
  void adjustDelayTimes();

  /** Causes the damping and correction filters to be updated. */
  void updateDampingAndCorrectionFilters();

  /** This function applies the feedback-matrix by means of a fast matrix multiplication
  algorithms, for example via the fast Hadamard-Transform (FHT). */
  INLINE void applyTheFeedbackMatrix();

  // maximum and actual number of delay-lines:
  static const int maxNumDelayLines = 16;
  int              numDelayLines;

  float impulseAmplitude; // if nonzero we feed in an impulse in the next call to getSample

  // maximum and actual delays in samples 
  // (maximum: 2 seconds @ sampleRate = 192 kHz):
  static const int maxDelayInSamples = 384000;
  //static const int maxDelayInSamples = 8192; // for testing only
  int              delaysInSamples[maxNumDelayLines];

  // the read and write pointers:
  int tapIns[maxNumDelayLines];  // the read-positions
  int tapOuts[maxNumDelayLines]; // the write-positions

  // the actual delaylines:
  double delayLines[maxNumDelayLines][maxDelayInSamples];

  doubleA delayLineOuts[maxNumDelayLines];
  // output signals of the individual delay-lines

  doubleA dryVol, wetVol; // volume of dry and wet signal



  double referenceDelayTime;
  double relativeDelayTimes[maxNumDelayLines];


  // the injection and output-vectors:
  double injectionVectorL[maxNumDelayLines];
  double injectionVectorR[maxNumDelayLines];
  double outputVectorL[maxNumDelayLines];
  double outputVectorR[maxNumDelayLines];



  // switches: 
  bool stereoSwapSwitch;
  bool wetPinking;       // indicates if wet signal should be pinkened
  bool allpassMode;      // when true, the delaylines will operate as 
  // allpass-filters instead of comb-filters

  // the injection- and output-vector and feedback matrix indices:
  int injectionVectorIndex;
  int feedbackMatrixIndex;
  int outputVectorIndex;

  // index the ordering/permutation of the delaylines:
  int delayOrdering;

  // the form parameter for certain distributions:
  //int delayDistributionIndex;
  //double delayDistributionForm;
  //double

  // embedded objects:
  DampingFilter dampingFilters[maxNumDelayLines];
  // these are the filters applied to the outputs of the delay-lines (before
  // the feedback-loop).

  DampingFilter correctionFilterL, correctionFilterR;
  // these filters have the inverse frequency-response of the damping-filters
  // in the feedback loop and compensate the overall loop-gain in order to 
  // decouple the reverb-time from the loudness of the wet signal. They are
  // used to filter the output-signal.

  LowpassHighpass wetFilterL, wetFilterR;
  // the combined lowpass- and highpass-filter for the wet signal fo left and
  // right channel

  WhiteToPinkFilter pinkingFilterL, pinkingFilterR;
  // the filters to pinken the wet-signal, if this option is chosen.

  //PrimeNumbers primeNumbers;
  // this object is used to determine the actual delay-times which should be 
  // prime-numbers

  doubleA sampleRate; // the sample-rate

  //doubleA minDelayTime, maxDelayTime; 
  // desired minimum and maximum delayline-length in milliseconds. the actual 
  // minimum and maximum delayline-lengths will be chosen to be prime-numbers
  // lower than or equal to the number of samples which corresponds to these
  // time-values.

  doubleA midReverbTime;
  // time, it takes for the output to decay to -60 dB for the mid-band

  doubleA lowReverbTimeScale, highReverbTimeScale;
  // scale-factors for the reverb times for low and high-bands with 
  // respect to reverb time for the mid mid-band

  doubleA lowCrossoverFreq, highCrossoverFreq;
  // these are the crossover-points, they are chosen to be the frequencies 
  // at which the reverb-time assumes the geometric mean between low/mid or
  // mid/high respectively

  doubleA feedbackGains[maxNumDelayLines];
  // gain-factors applied after each delay-line, these are needed only in the 
  // allpass mode because in the comb-mode they are absorbed into the 
  // dampingFilters

};

//-----------------------------------------------------------------------------------------------
// from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
// to be called at audio-rate (they can't be put into the .cpp file):

INLINE void FeedbackDelayNetwork16::getSampleFrameStereo(double* inOutL, double* inOutR)
{
  // possibly add the impulse (if nonzero) for auditioning purposes:
  double tmpL = *inOutL + impulseAmplitude;
  double tmpR = *inOutR + impulseAmplitude;

  if(impulseAmplitude != 0.f)
  {
    //DEBUG_BREAK;
  }

  impulseAmplitude = 0.f;


  static doubleA wetL, wetR, tmp;
  static intA    d;                  // index for delayline

  for(d=0; d<numDelayLines; d++)
  {
    // increment the pointers:
    tapIns[d]++;
    tapOuts[d]++;

    // wraparound the pointers:
    tapIns[d]  %= maxDelayInSamples;
    tapOuts[d] %= maxDelayInSamples;

    // pull out the outputs of the delaylines:
    delayLineOuts[d] = delayLines[d][tapOuts[d]];

    // apply damping-filters to the delayline-outputs:
    delayLineOuts[d] = dampingFilters[d].getSample(delayLineOuts[d]); //TEMPORARILY COMMENTED

    // stuff the new input into the delaylines via the injection-vectors:
    delayLines[d][tapIns[d]] = injectionVectorL[d]*tmpL + injectionVectorR[d]*tmpR;
  }

  // stuff also the current outputs into the delaylines via the 
  // feedback matrix:
  applyTheFeedbackMatrix();

  // obtain the left and right wet output-sample via the output-vectors:
  wetL = 0.0;
  wetR = 0.0;
  if(allpassMode == false) // delaylines work as comb-filters
  {
    for(d=0; d<numDelayLines; d++)
    {
      wetL += outputVectorL[d] * delayLineOuts[d];
      wetR += outputVectorR[d] * delayLineOuts[d];
    }
  }
  else                       // delaylines work as allpass-filters
  {
    static double g;
    for(d=0; d<numDelayLines; d++)
    {
      tmp = delayLines[d][tapIns[d]];
      // this is the signal which has gone into delayline d

      g =  feedbackGains[d];
      // this is the feedback-factor applied after delay-line d (not taking 
      // into account frequency dependent damping)

      wetL += outputVectorL[d] * ((1-g*g)*delayLineOuts[d] - g*tmp);
      wetR += outputVectorR[d] * ((1-g*g)*delayLineOuts[d] - g*tmp);
    }
  }

  // apply the compensation gain (or - as a refinement - the correction-filter):
  wetL = correctionFilterL.getSample(wetL);
  wetR = correctionFilterR.getSample(wetR);

  // pinken the wet signal, if pinking is chosen:
  if(wetPinking == true)
  {
    wetL = pinkingFilterL.getSample(wetL);
    wetR = pinkingFilterR.getSample(wetR);
  }

  // apply the lowpass- and highpass-filtering to the wet signal:
  wetL = wetFilterL.getSample(wetL);
  wetR = wetFilterR.getSample(wetR);

  // swap left and right, if selected:
  if(stereoSwapSwitch == true)
  {
    tmp  = wetL;
    wetL = wetR;
    wetR = tmp;
  }

  // mix dry and wet:
  *inOutL = dryVol*(*inOutL) + 0.25*wetVol*wetL;
  *inOutR = dryVol*(*inOutR) + 0.25*wetVol*wetR;
    // 0.25 (=1/sqrt(16)) compensates for the energy-addition of the 16 delayline-outputs
}

INLINE void FeedbackDelayNetwork16::applyTheFeedbackMatrix()
{
  switch(feedbackMatrixIndex)
  {
  case IDENTITY:
  {
    // apply an 16x16 identity-matrix: 
    //  + 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 + 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 + 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 + 0 0 0 0 0 0 0 0 0 0 0 
    //  0 0 0 0 0 + 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 + 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 + 0 0 0 0 0 0 0 0 
    //  0 0 0 0 0 0 0 0 + 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 + 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 + 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 + 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 + 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 + 0 0 
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 + 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 +
    // this corresponds to a parallel connection

    delayLines[0][tapIns[0]] += delayLineOuts[0];
    delayLines[1][tapIns[1]] += delayLineOuts[1];
    delayLines[2][tapIns[2]] += delayLineOuts[2];
    delayLines[3][tapIns[3]] += delayLineOuts[3];
    delayLines[4][tapIns[4]] += delayLineOuts[4];
    delayLines[5][tapIns[5]] += delayLineOuts[5];
    delayLines[6][tapIns[6]] += delayLineOuts[6];
    delayLines[7][tapIns[7]] += delayLineOuts[7];
    delayLines[8][tapIns[8]] += delayLineOuts[8];
    delayLines[9][tapIns[9]] += delayLineOuts[9];
    delayLines[10][tapIns[10]] += delayLineOuts[10];
    delayLines[11][tapIns[11]] += delayLineOuts[11];
    delayLines[12][tapIns[12]] += delayLineOuts[12];
    delayLines[13][tapIns[13]] += delayLineOuts[13];
    delayLines[14][tapIns[14]] += delayLineOuts[14];
    delayLines[15][tapIns[15]] += delayLineOuts[15];
  }
  break;

  case MINUS_IDENTITY:
  {
    // apply an 16x16 negated identity-matrix: 
    //  - 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 - 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 - 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 - 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 - 0 0 0 0 0 0 0 0 0 0 0 
    //  0 0 0 0 0 - 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 - 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 - 0 0 0 0 0 0 0 0 
    //  0 0 0 0 0 0 0 0 - 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 - 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 - 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 - 0 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 - 0 0 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 - 0 0 
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 - 0
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -    
    // this corresponds to a parallel connection

    delayLines[0][tapIns[0]] -= delayLineOuts[0];
    delayLines[1][tapIns[1]] -= delayLineOuts[1];
    delayLines[2][tapIns[2]] -= delayLineOuts[2];
    delayLines[3][tapIns[3]] -= delayLineOuts[3];
    delayLines[4][tapIns[4]] -= delayLineOuts[4];
    delayLines[5][tapIns[5]] -= delayLineOuts[5];
    delayLines[6][tapIns[6]] -= delayLineOuts[6];
    delayLines[7][tapIns[7]] -= delayLineOuts[7];
    delayLines[8][tapIns[8]] -= delayLineOuts[8];
    delayLines[9][tapIns[9]] -= delayLineOuts[9];
    delayLines[10][tapIns[10]] -= delayLineOuts[10];
    delayLines[11][tapIns[11]] -= delayLineOuts[11];
    delayLines[12][tapIns[12]] -= delayLineOuts[12];
    delayLines[13][tapIns[13]] -= delayLineOuts[13];
    delayLines[14][tapIns[14]] -= delayLineOuts[14];
    delayLines[15][tapIns[15]] -= delayLineOuts[15];
  }
  break;

  case SERIES_CONNECTION:
  {
    // apply an 16x16 matrix which corresponds to a series connection:
    //  0 + 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 + 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 + 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 + 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 + 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 + 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 + 0 0 0 0 0 0 0         
    //  0 0 0 0 0 0 0 0 0 + 0 0 0 0 0 0     
    //  0 0 0 0 0 0 0 0 0 0 + 0 0 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 + 0 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 + 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 + 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 + 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 +
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  

    delayLines[1][tapIns[1]] += delayLineOuts[0];
    delayLines[2][tapIns[2]] += delayLineOuts[1];
    delayLines[3][tapIns[3]] += delayLineOuts[2];
    delayLines[4][tapIns[4]] += delayLineOuts[3];
    delayLines[5][tapIns[5]] += delayLineOuts[4];
    delayLines[6][tapIns[6]] += delayLineOuts[5];
    delayLines[7][tapIns[7]] += delayLineOuts[6];
    delayLines[8][tapIns[8]] += delayLineOuts[7];
    delayLines[9][tapIns[9]] += delayLineOuts[8];
    delayLines[10][tapIns[10]] += delayLineOuts[9];
    delayLines[11][tapIns[11]] += delayLineOuts[10];
    delayLines[12][tapIns[12]] += delayLineOuts[11];
    delayLines[13][tapIns[13]] += delayLineOuts[12];
    delayLines[14][tapIns[14]] += delayLineOuts[13];
    delayLines[15][tapIns[15]] += delayLineOuts[14];
  }
  break;

  case SERIES_WITH_GLOBAL_FEEDBACK:
  {
    // apply an 16x16 matrix which corresponds to a series connection with 
    // global feedback:
    //  0 + 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 + 0 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 + 0 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 + 0 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 + 0 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 + 0 0 0 0 0 0 0 0
    //  0 0 0 0 0 0 0 0 + 0 0 0 0 0 0 0         
    //  0 0 0 0 0 0 0 0 0 + 0 0 0 0 0 0     
    //  0 0 0 0 0 0 0 0 0 0 + 0 0 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 + 0 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 + 0 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 + 0 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 + 0  
    //  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 +
    //  + 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0        

    delayLines[0][tapIns[0]] += delayLineOuts[15];
    delayLines[1][tapIns[1]] += delayLineOuts[0];
    delayLines[2][tapIns[2]] += delayLineOuts[1];
    delayLines[3][tapIns[3]] += delayLineOuts[2];
    delayLines[4][tapIns[4]] += delayLineOuts[3];
    delayLines[5][tapIns[5]] += delayLineOuts[4];
    delayLines[6][tapIns[6]] += delayLineOuts[5];
    delayLines[7][tapIns[7]] += delayLineOuts[6];
    delayLines[8][tapIns[8]] += delayLineOuts[7];
    delayLines[9][tapIns[9]] += delayLineOuts[8];
    delayLines[10][tapIns[10]] += delayLineOuts[9];
    delayLines[11][tapIns[11]] += delayLineOuts[10];
    delayLines[12][tapIns[12]] += delayLineOuts[11];
    delayLines[13][tapIns[13]] += delayLineOuts[12];
    delayLines[14][tapIns[14]] += delayLineOuts[13];
    delayLines[15][tapIns[15]] += delayLineOuts[14];
  }
  break;

  case HADAMARD:
  {
    // apply an 16x16 Hadamard-matrix --- obsoloete - edit this: 
    //  + + + + + + + +
    //  + - + - + - + -
    //  + + - - + + - -
    //  + - - + + - - +
    //  + + + + - - - -   * (1/sqrt(8))
    //  + - + - - + - +
    //  + + - - - - + +
    //  + - - + - + + -         
    // by means of the fast Hadamard transform, + is +1, - is -1

    // assign vectors for intermediate results:
    double a[maxNumDelayLines];
    double b[maxNumDelayLines];

    // do the Fast Hadamard Transform:
    a[0] = delayLineOuts[0] + delayLineOuts[1];
    a[1] = delayLineOuts[2] + delayLineOuts[3];
    a[2] = delayLineOuts[4] + delayLineOuts[5];
    a[3] = delayLineOuts[6] + delayLineOuts[7];
    a[4] = delayLineOuts[8] + delayLineOuts[9];
    a[5] = delayLineOuts[10] + delayLineOuts[11];
    a[6] = delayLineOuts[12] + delayLineOuts[13];
    a[7] = delayLineOuts[14] + delayLineOuts[15];
    a[8] = delayLineOuts[0] - delayLineOuts[1];
    a[9] = delayLineOuts[2] - delayLineOuts[3];
    a[10] = delayLineOuts[4] - delayLineOuts[5];
    a[11] = delayLineOuts[6] - delayLineOuts[7];
    a[12] = delayLineOuts[8] - delayLineOuts[9];
    a[13] = delayLineOuts[10] - delayLineOuts[11];
    a[14] = delayLineOuts[12] - delayLineOuts[13];
    a[15] = delayLineOuts[14] - delayLineOuts[15];

    b[0] = a[0] + a[1];
    b[1] = a[2] + a[3];
    b[2] = a[4] + a[5];
    b[3] = a[6] + a[7];
    b[4] = a[8] + a[9];
    b[5] = a[10] + a[11];
    b[6] = a[12] + a[13];
    b[7] = a[14] + a[15];
    b[8] = a[0] - a[1];
    b[9] = a[2] - a[3];
    b[10] = a[4] - a[5];
    b[11] = a[6] - a[7];
    b[12] = a[8] - a[9];
    b[13] = a[10] - a[11];
    b[14] = a[12] - a[13];
    b[15] = a[14] - a[15];

    a[0] = b[0] + b[1];
    a[1] = b[2] + b[3];
    a[2] = b[4] + b[5];
    a[3] = b[6] + b[7];
    a[4] = b[8] + b[9];
    a[5] = b[10] + b[11];
    a[6] = b[12] + b[13];
    a[7] = b[14] + b[15];
    a[8] = b[0] - b[1];
    a[9] = b[2] - b[3];
    a[10] = b[4] - b[5];
    a[11] = b[6] - b[7];
    a[12] = b[8] - b[9];
    a[13] = b[10] - b[11];
    a[14] = b[12] - b[13];
    a[15] = b[14] - b[15];

    delayLines[0][tapIns[0]] += ONE_OVER_SQRT16 * (a[0] + a[1]);
    delayLines[1][tapIns[1]] += ONE_OVER_SQRT16 * (a[2] + a[3]);
    delayLines[2][tapIns[2]] += ONE_OVER_SQRT16 * (a[4] + a[5]);
    delayLines[3][tapIns[3]] += ONE_OVER_SQRT16 * (a[6] + a[7]);
    delayLines[4][tapIns[4]] += ONE_OVER_SQRT16 * (a[8] + a[9]);
    delayLines[5][tapIns[5]] += ONE_OVER_SQRT16 * (a[10] + a[11]);
    delayLines[6][tapIns[6]] += ONE_OVER_SQRT16 * (a[12] + a[13]);
    delayLines[7][tapIns[7]] += ONE_OVER_SQRT16 * (a[14] + a[15]);
    delayLines[8][tapIns[8]] += ONE_OVER_SQRT16 * (a[0] - a[1]);
    delayLines[9][tapIns[9]] += ONE_OVER_SQRT16 * (a[2] - a[3]);
    delayLines[10][tapIns[10]] += ONE_OVER_SQRT16 * (a[4] - a[5]);
    delayLines[11][tapIns[11]] += ONE_OVER_SQRT16 * (a[6] - a[7]);
    delayLines[12][tapIns[12]] += ONE_OVER_SQRT16 * (a[8] - a[9]);
    delayLines[13][tapIns[13]] += ONE_OVER_SQRT16 * (a[10] - a[11]);
    delayLines[14][tapIns[14]] += ONE_OVER_SQRT16 * (a[12] - a[13]);
    delayLines[15][tapIns[15]] += ONE_OVER_SQRT16 * (a[14] - a[15]);
  }
  break;

  case MINUS_HADAMARD:
  {
    // Comment seems obsolete - ...and is it really worth to have an extra branch for that or could we 
    // just use positive Hadamard branch and then negate afterwards?
    // apply an 8x8 negated Hadamard-matrix: 
    //  - - - - - - - -
    //  - + - + - + - +
    //  - - + + - - + +
    //  - + + - - + + -
    //  - - - - + + + +   * (1/sqrt(8))
    //  - + - + + - + -
    //  - - + + + + - - 
    //  - + + - + - - +         
    // by means of the fast Hadamard transform, + is +1, - is -1

    // assign vectors for intermediate results:
    double a[maxNumDelayLines];
    double b[maxNumDelayLines];

    // do the Fast Hadamard Transform:
    a[0] = delayLineOuts[0] + delayLineOuts[1];
    a[1] = delayLineOuts[2] + delayLineOuts[3];
    a[2] = delayLineOuts[4] + delayLineOuts[5];
    a[3] = delayLineOuts[6] + delayLineOuts[7];
    a[4] = delayLineOuts[8] + delayLineOuts[9];
    a[5] = delayLineOuts[10] + delayLineOuts[11];
    a[6] = delayLineOuts[12] + delayLineOuts[13];
    a[7] = delayLineOuts[14] + delayLineOuts[15];
    a[8] = delayLineOuts[0] - delayLineOuts[1];
    a[9] = delayLineOuts[2] - delayLineOuts[3];
    a[10] = delayLineOuts[4] - delayLineOuts[5];
    a[11] = delayLineOuts[6] - delayLineOuts[7];
    a[12] = delayLineOuts[8] - delayLineOuts[9];
    a[13] = delayLineOuts[10] - delayLineOuts[11];
    a[14] = delayLineOuts[12] - delayLineOuts[13];
    a[15] = delayLineOuts[14] - delayLineOuts[15];

    b[0] = a[0] + a[1];
    b[1] = a[2] + a[3];
    b[2] = a[4] + a[5];
    b[3] = a[6] + a[7];
    b[4] = a[8] + a[9];
    b[5] = a[10] + a[11];
    b[6] = a[12] + a[13];
    b[7] = a[14] + a[15];
    b[8] = a[0] - a[1];
    b[9] = a[2] - a[3];
    b[10] = a[4] - a[5];
    b[11] = a[6] - a[7];
    b[12] = a[8] - a[9];
    b[13] = a[10] - a[11];
    b[14] = a[12] - a[13];
    b[15] = a[14] - a[15];

    a[0] = b[0] + b[1];
    a[1] = b[2] + b[3];
    a[2] = b[4] + b[5];
    a[3] = b[6] + b[7];
    a[4] = b[8] + b[9];
    a[5] = b[10] + b[11];
    a[6] = b[12] + b[13];
    a[7] = b[14] + b[15];
    a[8] = b[0] - b[1];
    a[9] = b[2] - b[3];
    a[10] = b[4] - b[5];
    a[11] = b[6] - b[7];
    a[12] = b[8] - b[9];
    a[13] = b[10] - b[11];
    a[14] = b[12] - b[13];
    a[15] = b[14] - b[15];

    delayLines[0][tapIns[0]] -= ONE_OVER_SQRT16 * (a[0] + a[1]);
    delayLines[1][tapIns[1]] -= ONE_OVER_SQRT16 * (a[2] + a[3]);
    delayLines[2][tapIns[2]] -= ONE_OVER_SQRT16 * (a[4] + a[5]);
    delayLines[3][tapIns[3]] -= ONE_OVER_SQRT16 * (a[6] + a[7]);
    delayLines[4][tapIns[4]] -= ONE_OVER_SQRT16 * (a[8] + a[9]);
    delayLines[5][tapIns[5]] -= ONE_OVER_SQRT16 * (a[10] + a[11]);
    delayLines[6][tapIns[6]] -= ONE_OVER_SQRT16 * (a[12] + a[13]);
    delayLines[7][tapIns[7]] -= ONE_OVER_SQRT16 * (a[14] + a[15]);
    delayLines[8][tapIns[8]] -= ONE_OVER_SQRT16 * (a[0] - a[1]);
    delayLines[9][tapIns[9]] -= ONE_OVER_SQRT16 * (a[2] - a[3]);
    delayLines[10][tapIns[10]] -= ONE_OVER_SQRT16 * (a[4] - a[5]);
    delayLines[11][tapIns[11]] -= ONE_OVER_SQRT16 * (a[6] - a[7]);
    delayLines[12][tapIns[12]] -= ONE_OVER_SQRT16 * (a[8] - a[9]);
    delayLines[13][tapIns[13]] -= ONE_OVER_SQRT16 * (a[10] - a[11]);
    delayLines[14][tapIns[14]] -= ONE_OVER_SQRT16 * (a[12] - a[13]);
    delayLines[15][tapIns[15]] -= ONE_OVER_SQRT16 * (a[14] - a[15]);
  }
  break;
  } // end of switch( feedbackMatrixIndex )

}

} // end namespace rosic

#endif // rosic_FeedbackDelayNetwork16_h
