#ifndef rosic_Reverb_h
#define rosic_Reverb_h

namespace rosic
{

/** This class implements a feedback-delay-network based reverb with 16 delay-lines. Each of the 16
delaylines has a damping-filter in series which serves as a frequency dependent feedback gain. The 
feedback itself is realized via a Hadamard-matrix which is implemented via the Fast Hadamard 
Transform. The injection- and output-vectors are pre-defined. The output-signals of the delay-lines
are fed through compensation-filters before actually being routed to the output. These filters 
decouple the overall reverberation-time from the output volume - and they even do this in a 
frequency dependent way. Roughly speaking, they act opposite to the damping-filters in order to
decoulour (whiten) the wet signal in a steady-state condition. */

class rsReverb
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor - will allocate the given amount of memory (in samples) in total for all 16 
  delaylines. A value of 1048576 will support maximum delayline lengths of more than 340 ms with 
  samplerates up to 192 kHz. */
  rsReverb(int delayMemoryInSamplesToAllocate = 1048576);

  /** Destructor. */
  ~rsReverb();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

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

  /** Sets the ratio between dry and wet between 0...1. */
  void setDryWetRatio(double newDryWetRatio)
  {
    RAPT::rsEqualPowerGainFactors(newDryWetRatio, &dryVol, &wetVol, 0.0, 1.0); wetVol *= 0.25;
  }
// factor 0.25 (=1/sqrt(16)) compensates for the energy-addition of the 16 delayline-outputs

/** Sets the cutoff frequency of a lowpass filter for the wet signal. */
  void setWetLowpassCutoff(double newCutoff);

  /** Sets the cutoff frequency of a highpass filter for the wet signal. */
  void setWetHighpassCutoff(double newCutoff);

  /** Switches swapping of left and right output channel on and off. */
  void setStereoSwapSwitch(bool newStereoSwapSwitch);

  /** Switches a pinking filter for the wet signal on and off. */
  void setWetPinkingSwitch(bool newWetPinkingSwitch);

  /** Sets a reference delay time in milliseconds - the actual delay-times (in milliseconds) of
  the delaylines may then be specified relative to this value via setRelativeDelayTime. */
  void setReferenceDelayTime(double newReferenceDelayTime);

  /** Sets up a pred-delay time in milliseconds. */
  void setPreDelay(double newPreDelay)
  {
    preDelayLineL.setDelayInMilliseconds(newPreDelay);
    preDelayLineR.setDelayInMilliseconds(newPreDelay);
  }

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers to zero. */
  void reset();


protected:

  /** sets up the output-vectors for left and right channel. */
  void setupOutputVector();

  /** Assigns the relative delaytimes. */
  void assignRelativeDelayTimes();

  /** Triggers a re-adjustment of one of the read-pointers - this is necesarry when the delay-time 
  of the corresponding delay-line  changes. It is also necessarry to do this for all delay-lines 
  when the sample-rate changes. */
  void adjustReadPointer(int index);

  /** Adjusts the lengths of the delay-lines. */
  void adjustDelayTimes();

  /** Causes the damping and correction filters to be updated. */
  void updateDampingAndCorrectionFilters();

  /** Frees the allocated memory, in case it is not shared memory. */
  void freeMemoryIfNotShared();

  /** This function applies the feedback-matrix by means of the fast Hadamard-Transform (FHT). */
  INLINE void applyTheFeedbackMatrix();

  static const int numDelayLines = 16;      // number of delay-lines
  double *delayLines[numDelayLines];        // the delaylines themselves
  double delayLineOuts[numDelayLines];      // outputs of the individual delay-lines
  double dryVol, wetVol;                    // gain factor for dry and wet signal

  int    delaysInSamples[numDelayLines];    // delay per delayline (in samples)
  int    maxDelayInSamples;                 // maximum delay per delayline
  int    tapIn;                             // write-position
  int    tapOuts[numDelayLines];            // read-positions

  double relativeDelayTimes[numDelayLines]; // relative delays with respect to referenceDelayTime
  double outputVectorL[numDelayLines];      // output vector for left channel
  double outputVectorR[numDelayLines];      // output vector for right channel
  double referenceDelayTime;                // some nominal reference delaytime 
  double sampleRate;                        // the sample-rate
  double midReverbTime;                     // time, it takes for the output to decay to -60 dB
  double lowReverbTimeScale;                // scale-factor for the reverb times for low band
  double highReverbTimeScale;               // scale-factor for the reverb times for high band
  double lowCrossoverFreq;                  // crossover frequency between low and mid band  
  double highCrossoverFreq;                 // crossover frequency between mid and high band  

    // the crossover frequencies are chosen to be the frequencies at which the reverb-time 
    // assumes the geometric mean between low/mid or mid/high respectively

  bool   stereoSwapSwitch;                  // indicates if stereo channels should be swapped
  bool   wetPinking;                        // indicates if wet signal should be pinkened

  // delaylines for the pre-delay:
  IntegerDelayLine preDelayLineL, preDelayLineR;

  // filters to be applied to the outputs of the delay-lines (before the feedback-loop):
  DampingFilter dampingFilters[numDelayLines];

  // These filters have the inverse frequency-response of the damping-filters in the feedback 
  // loop and compensate the overall loop-gain  in order to decouple the reverb-time from the 
  // loudness of the wet signal. They are used to filter the output-signal:
  DampingFilter correctionFilterL, correctionFilterR;

  // Filters for the wet signal fo left and right channel:
  LowpassHighpass wetFilterL, wetFilterR;

  // Filters to pinken the wet-signal (if this option is chosen):
  WhiteToPinkFilter pinkingFilterL, pinkingFilterR;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsReverb::getSampleFrameStereo(double* inOutL, double* inOutR)
{
  double tmpL = preDelayLineL.getSample(*inOutL);
  double tmpR = preDelayLineR.getSample(*inOutR);

  // do the I/O for the delaylines:
  tapIn++;
  if(tapIn >= maxDelayInSamples)
    tapIn = 0;
  int d;
  for(d=0; d<numDelayLines; d++)
  {
    tapOuts[d]++;
    if(tapOuts[d] >= maxDelayInSamples)
      tapOuts[d] = 0;
    delayLineOuts[d]     = dampingFilters[d].getSample(delayLines[d][tapOuts[d]]);
    delayLines[d][tapIn] = SQRT2_INV*(tmpL+tmpR);
  }
  applyTheFeedbackMatrix();

  // obtain (via the output vectors) and post-process the wet signal:
  double wetL = 0.0;
  double wetR = 0.0;
  for(d=0; d<numDelayLines; d++)
  {
    wetL += outputVectorL[d] * delayLineOuts[d];
    wetR += outputVectorR[d] * delayLineOuts[d];
  }
  wetL = correctionFilterL.getSample(wetL);
  wetR = correctionFilterR.getSample(wetR);
  if(wetPinking == true)
  {
    wetL = pinkingFilterL.getSample(wetL);
    wetR = pinkingFilterR.getSample(wetR);
  }
  wetL = wetFilterL.getSample(wetL);
  wetR = wetFilterR.getSample(wetR);
  if(stereoSwapSwitch == true)
  {
    double tmp  = wetL;
    wetL        = wetR;
    wetR        = tmp;
  }

  // mix dry/wet and output:
  *inOutL = dryVol*(*inOutL) + wetVol*wetL;
  *inOutR = dryVol*(*inOutR) + wetVol*wetR;
}

INLINE void rsReverb::applyTheFeedbackMatrix()
{
  // apply an 16x16 Hadamard-matrix by means of the fast Hadamard transform:

  // assign vectors for intermediate results:
  double a[numDelayLines];
  double b[numDelayLines];

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

  delayLines[0][tapIn] += 0.25 * (a[0] + a[1]);
  delayLines[1][tapIn] += 0.25 * (a[2] + a[3]);
  delayLines[2][tapIn] += 0.25 * (a[4] + a[5]);
  delayLines[3][tapIn] += 0.25 * (a[6] + a[7]);
  delayLines[4][tapIn] += 0.25 * (a[8] + a[9]);
  delayLines[5][tapIn] += 0.25 * (a[10] + a[11]);
  delayLines[6][tapIn] += 0.25 * (a[12] + a[13]);
  delayLines[7][tapIn] += 0.25 * (a[14] + a[15]);
  delayLines[8][tapIn] += 0.25 * (a[0] - a[1]);
  delayLines[9][tapIn] += 0.25 * (a[2] - a[3]);
  delayLines[10][tapIn] += 0.25 * (a[4] - a[5]);
  delayLines[11][tapIn] += 0.25 * (a[6] - a[7]);
  delayLines[12][tapIn] += 0.25 * (a[8] - a[9]);
  delayLines[13][tapIn] += 0.25 * (a[10] - a[11]);
  delayLines[14][tapIn] += 0.25 * (a[12] - a[13]);
  delayLines[15][tapIn] += 0.25 * (a[14] - a[15]);
}

} // end namespace rosic

#endif // rosic_Reverb_h
