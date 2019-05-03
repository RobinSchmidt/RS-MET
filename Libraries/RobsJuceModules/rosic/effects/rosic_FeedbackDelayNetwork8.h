#ifndef rosic_FeedbackDelayNetwork8_h
#define rosic_FeedbackDelayNetwork8_h

//// rosic-indcludes:
//#include "../math/rosic_PrimeNumbers.h"
//#include "../filters/rosic_DampingFilter.h"
//#include "../filters/rosic_LowpassHighpass.h"
//#include "../filters/rosic_WhiteToPinkFilter.h"

namespace rosic
{

 /**

 This class implements a feedback-delay-network with 8 delay-lines which is 
 supposed to be used mainly as a reverberation effect unit. Each of the 8 
 delaylines has a damping-filter in series which serves as a frequency 
 dependent feedback gain. The feedback itself is realized via various 
 pre-defined feedback-matrices, which are implemented by efficient 
 techniques (for, example a multiplication by a Hadamard-matrix is done via
 the Fast Hadamard Transform). The injection- and output-vectors can also be 
 chosen from a list of pre-defined vectors. The output-signals of the 
 delay-lines are fed through compensation-filters before actually being routed 
 to the output. These filters decouple the overall reverberation-time from the
 output volume - and they even do this in a frequency dependent way. Roughly 
 speaking, they act opposite to the damping-filters in order to decoulour 
 (whiten) the wet signal in a steady-state condition.

 */

 class FeedbackDelayNetwork8
 {

 public:

  #define ONE_OVER_SQRT8 0.35355339059327376220042218105242

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
   OUT_01 = 0,
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

  //---------------------------------------------------------------------------
  // construction/destruction:

  FeedbackDelayNetwork8();   ///< Constructor.
  ~FeedbackDelayNetwork8();  ///< Destructor.

  //---------------------------------------------------------------------------
  //parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Overrides the setSampleRate() method of the AudioModule base class.

  void setLowReverbTimeScale(double newLowReverbTimeScale);
  /**< Scales the reverb time for the low-band respect to the mid-band. */

  void setLowCrossoverFreq(double newLowCrossoverFreq);
  /**< Sets the crossover-frequency between low and mid frequencies. */

  void setMidReverbTime(double newMidReverbTime);
  /**< sets the time in which the amplitude decays to -60 dB for the 
       mid frequencies. */

  void setHighReverbTimeScale(double newHighReverbTimeScale);
  /**< Scales the reverb time for the high-band respect to the mid-band. */

  void setHighCrossoverFreq(double newHighCrossoverFreq);
  /**< Sets the crossover-frequency between mid and high frequencies. */

  void setDensity(int newDensity);
  //*< Selects, how many delay-lines should be used (4 or 8). */

  void setAllpassMode(bool shouldBeAllpass);
  //*< Selects, if the delay-lines should operate as comb- or as allpass-filter. */

  void setDryVolume(double newDryVolume);
  /**< Sets the volume of the dry signal in dB */

  void setWetVolume(double newWetVolume);
  /**< Sets the volume of the wet signal in dB */

  void setWetLpfCutoff(double newWetLpfCutoff);
  /**< Sets the cutoff frequency of a lowpass filter for the wet signal. */

  void setWetHpfCutoff(double newWetHpfCutoff);
  /**< Sets the cutoff frequency of a highpass filter for the wet signal. */

  void setStereoSwapSwitch(bool newStereoSwapSwitch);
  /**< Switches swapping of left and right output channel on and off. */

  void setWetPinkingSwitch(bool newWetPinkingSwitch);
  /**< Switches a pinking filter for the wet signal on and off. */

  void setMinDelayTime(double newMinDelayTime);
  /**< Chooses (approximately) the length of the shortest delay-line in 
       milliseconds - the exact delay-time will be chosen to be the number of
       samples which is the closest prime-number below the number of samples 
       corresponding to the value specified here. */

  void setMaxDelayTime(double newMaxDelayTime);
  /**< Chooses (approximately) the length of the longest delay-line in 
       milliseconds - the exact delay-time will be chosen to be the number of
       samples which is the closest prime-number below the number of samples 
       corresponding to the value specified here. */
         
  void setDelayDistributionIndex(int newDelayDistributionIndex);
  /**< Chooses one of the various distributions for the delay-times of the
       delay-lines - this governs how the delay-times are distributed between
       the minimum and maximum delay-time. */

  void setDelayDistributionForm(double newDelayDistributionForm);
  /**< Some distributions provide a form parameter witch which the distribution
       can be controlled. */

  void setDelayOrdering(int newDelayOrdering);
  /**< Chooses, how the delay-lines of different length are ordered with 
       respect to their length (ascending, descending, alternating, etc. */

  void setInjectionVector(int newInjectionVectorIndex);
  /**< Chooses one of the pre-defined injection-vectors. */

  void setFeedbackMatrix(int newFeedbackMatrixIndex);
  /**< Chooses one of the pre-defined feedback-matrices. */

  void setOutputVector(int newOutputVectorIndex);
  /**< Chooses one of the pre-defined output-vectors. */


  //---------------------------------------------------------------------------
  //audio processing:

  INLINE void getSampleFrameStereo(double* inL, 
                                   double* inR,
                                   double* outL,
                                   double* outR);  

  //---------------------------------------------------------------------------
  //others:

  void resetBuffers();
  ///< Resets the internal buffers to zero.

  //===========================================================================

 protected:

  // maximum and actual number of delay-lines:
  static const int maxNumDelayLines = 8;       
  int              numDelayLines;  

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

  // the actual injection and output-vectors:
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

  // indices for the delay-distribution and the ordering:
  int delayDistributionIndex;
  int delayOrdering;

  // the form parameter for certain distributions:
  double delayDistributionForm;

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

  PrimeNumbers primeNumbers;
   // this object is used to determine the actual delay-times which should be 
   // prime-numbers

  doubleA sampleRate; // the sample-rate

  doubleA minDelayTime, maxDelayTime; 
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

  // internal functions:
  void adjustReadPointer(int index);
   // triggers a re-adjustment of one of the read-pointers - this is necesarry
   // when the delay-time of the corresponding delay-line changes. It is also
   // necessarry to do this for all delay-lines when the sample-rate changes.

  void adjustDelayTimes();
   // adjusts the lengths of the delay-lines according to minDelayTime, 
   // maxDelayTime, delayDistribution and delayDistributionForm.

  void updateDampingAndCorrectionFilters();
   // causes the damping and correction filters to be updated.

  INLINE void applyTheFeedbackMatrix();  
   // this function applies the feedback-matrix by means of a fast 
   // Hadamard-Transform (FHT).
 };

 //-----------------------------------------------------------------------------
 //from here: definitions of the functions to be inlined, i.e. all functions
 //which are supposed to be called at audio-rate (they can't be put into
 //the .cpp file):

 INLINE void FeedbackDelayNetwork8::getSampleFrameStereo(double* inL, 
                                                         double* inR,
                                                         double* outL,
                                                         double* outR)
 {
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
   delayLineOuts[d] = dampingFilters[d].getSample(delayLineOuts[d]);

   // stuff the new input into the delaylines via the injection-vectors:
   delayLines[d][tapIns[d]] =   injectionVectorL[d] * (*inL) 
                              + injectionVectorR[d] * (*inR);
  }

  // stuff also the current outputs into the delaylines via the 
  // feedback matrix:
  applyTheFeedbackMatrix();

  // obtain the left and right wet output-sample via the output-vectors:
  wetL = 0.0;
  wetR = 0.0;
  if( allpassMode == false ) // delaylines work as comb-filters
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

    wetL += outputVectorL[d] * ((1-g*g)*delayLineOuts[d] - g*tmp) ;
    wetR += outputVectorR[d] * ((1-g*g)*delayLineOuts[d] - g*tmp) ;
   }
  }

  // apply the compensation gain (or - as a refinement - the correction-filter):
  wetL = correctionFilterL.getSample(wetL);
  wetR = correctionFilterR.getSample(wetR);

  // pinken the wet signal, if pinking is chosen:
  if( wetPinking == true)
  {
   wetL = pinkingFilterL.getSample(wetL);
   wetR = pinkingFilterR.getSample(wetR);
  }

  // apply the lowpass- and highpass-filtering to the wet signal:
  wetL = wetFilterL.getSample(wetL);
  wetR = wetFilterR.getSample(wetR);

  // swap left and right, if selected:
  if( stereoSwapSwitch == true )
  {
   tmp  = wetL;
   wetL = wetR;
   wetR = tmp;
  }

  // mix dry and wet:
  *outL = dryVol * (*inL)  +  wetVol * wetL;
  *outR = dryVol * (*inR)  +  wetVol * wetR;
 }

 INLINE void FeedbackDelayNetwork8::applyTheFeedbackMatrix()
 {


  switch( feedbackMatrixIndex )
  {
  case IDENTITY:
   {
    // apply an 8x8 identity-matrix: 
    //  + 0 0 0 0 0 0 0
    //  0 + 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0
    //  0 0 0 + 0 0 0 0
    //  0 0 0 0 + 0 0 0  
    //  0 0 0 0 0 + 0 0
    //  0 0 0 0 0 0 + 0
    //  0 0 0 0 0 0 0 +         
    // this corresponds to a parallel connection

    delayLines[0][tapIns[0]] += delayLineOuts[0];
    delayLines[1][tapIns[1]] += delayLineOuts[1];
    delayLines[2][tapIns[2]] += delayLineOuts[2];
    delayLines[3][tapIns[3]] += delayLineOuts[3];
    delayLines[4][tapIns[4]] += delayLineOuts[4];
    delayLines[5][tapIns[5]] += delayLineOuts[5];
    delayLines[6][tapIns[6]] += delayLineOuts[6];
    delayLines[7][tapIns[7]] += delayLineOuts[7];
   }
   break;

  case MINUS_IDENTITY:
   {
    // apply an 8x8 negated identity-matrix: 
    //  - 0 0 0 0 0 0 0
    //  0 - 0 0 0 0 0 0
    //  0 0 - 0 0 0 0 0
    //  0 0 0 - 0 0 0 0
    //  0 0 0 0 - 0 0 0  
    //  0 0 0 0 0 - 0 0
    //  0 0 0 0 0 0 - 0
    //  0 0 0 0 0 0 0 -         
    // this corresponds to a parallel connection

    delayLines[0][tapIns[0]] -= delayLineOuts[0];
    delayLines[1][tapIns[1]] -= delayLineOuts[1];
    delayLines[2][tapIns[2]] -= delayLineOuts[2];
    delayLines[3][tapIns[3]] -= delayLineOuts[3];
    delayLines[4][tapIns[4]] -= delayLineOuts[4];
    delayLines[5][tapIns[5]] -= delayLineOuts[5];
    delayLines[6][tapIns[6]] -= delayLineOuts[6];
    delayLines[7][tapIns[7]] -= delayLineOuts[7];
   }
   break;

  case SERIES_CONNECTION:
   {
    // apply an 8x8 matrix which corresponds to a series connection:
    //  0 + 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0
    //  0 0 0 + 0 0 0 0
    //  0 0 0 0 + 0 0 0
    //  0 0 0 0 0 + 0 0
    //  0 0 0 0 0 0 + 0
    //  0 0 0 0 0 0 0 +
    //  0 0 0 0 0 0 0 0         

    delayLines[1][tapIns[1]] += delayLineOuts[0];
    delayLines[2][tapIns[2]] += delayLineOuts[1];
    delayLines[3][tapIns[3]] += delayLineOuts[2];
    delayLines[4][tapIns[4]] += delayLineOuts[3];
    delayLines[5][tapIns[5]] += delayLineOuts[4];
    delayLines[6][tapIns[6]] += delayLineOuts[5];
    delayLines[7][tapIns[7]] += delayLineOuts[6];
   }
   break;

  case SERIES_WITH_GLOBAL_FEEDBACK:
   {
    // apply an 8x8 matrix which corresponds to a series connection:
    //  0 + 0 0 0 0 0 0
    //  0 0 + 0 0 0 0 0
    //  0 0 0 + 0 0 0 0
    //  0 0 0 0 + 0 0 0
    //  0 0 0 0 0 + 0 0
    //  0 0 0 0 0 0 + 0
    //  0 0 0 0 0 0 0 +
    //  + 0 0 0 0 0 0 0         

    delayLines[0][tapIns[0]] += delayLineOuts[7];
    delayLines[1][tapIns[1]] += delayLineOuts[0];
    delayLines[2][tapIns[2]] += delayLineOuts[1];
    delayLines[3][tapIns[3]] += delayLineOuts[2];
    delayLines[4][tapIns[4]] += delayLineOuts[3];
    delayLines[5][tapIns[5]] += delayLineOuts[4];
    delayLines[6][tapIns[6]] += delayLineOuts[5];
    delayLines[7][tapIns[7]] += delayLineOuts[6];
   }
   break;

  case HADAMARD:
   {
    // apply an 8x8 Hadamard-matrix: 
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
    a[4] = delayLineOuts[0] - delayLineOuts[1];
    a[5] = delayLineOuts[2] - delayLineOuts[3];
    a[6] = delayLineOuts[4] - delayLineOuts[5];
    a[7] = delayLineOuts[6] - delayLineOuts[7];

    b[0] = a[0] + a[1];
    b[1] = a[2] + a[3];
    b[2] = a[4] + a[5];
    b[3] = a[6] + a[7];
    b[4] = a[0] - a[1];
    b[5] = a[2] - a[3];
    b[6] = a[4] - a[5];
    b[7] = a[6] - a[7];

    delayLines[0][tapIns[0]] += ONE_OVER_SQRT8 * ( b[0] + b[1] );
    delayLines[1][tapIns[1]] += ONE_OVER_SQRT8 * ( b[2] + b[3] );
    delayLines[2][tapIns[2]] += ONE_OVER_SQRT8 * ( b[4] + b[5] );
    delayLines[3][tapIns[3]] += ONE_OVER_SQRT8 * ( b[6] + b[7] );
    delayLines[4][tapIns[4]] += ONE_OVER_SQRT8 * ( b[0] - b[1] );
    delayLines[5][tapIns[5]] += ONE_OVER_SQRT8 * ( b[2] - b[3] );
    delayLines[6][tapIns[6]] += ONE_OVER_SQRT8 * ( b[4] - b[5] );
    delayLines[7][tapIns[7]] += ONE_OVER_SQRT8 * ( b[6] - b[7] );
   }
   break;

  case MINUS_HADAMARD:
   {
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
    a[4] = delayLineOuts[0] - delayLineOuts[1];
    a[5] = delayLineOuts[2] - delayLineOuts[3];
    a[6] = delayLineOuts[4] - delayLineOuts[5];
    a[7] = delayLineOuts[6] - delayLineOuts[7];

    b[0] = a[0] + a[1];
    b[1] = a[2] + a[3];
    b[2] = a[4] + a[5];
    b[3] = a[6] + a[7];
    b[4] = a[0] - a[1];
    b[5] = a[2] - a[3];
    b[6] = a[4] - a[5];
    b[7] = a[6] - a[7];

    delayLines[0][tapIns[0]] -= ONE_OVER_SQRT8 * ( b[0] + b[1] );
    delayLines[1][tapIns[1]] -= ONE_OVER_SQRT8 * ( b[2] + b[3] );
    delayLines[2][tapIns[2]] -= ONE_OVER_SQRT8 * ( b[4] + b[5] );
    delayLines[3][tapIns[3]] -= ONE_OVER_SQRT8 * ( b[6] + b[7] );
    delayLines[4][tapIns[4]] -= ONE_OVER_SQRT8 * ( b[0] - b[1] );
    delayLines[5][tapIns[5]] -= ONE_OVER_SQRT8 * ( b[2] - b[3] );
    delayLines[6][tapIns[6]] -= ONE_OVER_SQRT8 * ( b[4] - b[5] );
    delayLines[7][tapIns[7]] -= ONE_OVER_SQRT8 * ( b[6] - b[7] );
   }
   break;
  } // end of switch( feedbackMatrixIndex )

 }

} // end namespace rosic

#endif // rosic_FeedbackDelayNetwork8_h
