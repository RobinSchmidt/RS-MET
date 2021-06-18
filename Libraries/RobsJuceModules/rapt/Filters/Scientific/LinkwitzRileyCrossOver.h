#ifndef RAPT_LINKWITZRILEYCROSSOVER_H_INCLUDED
#define RAPT_LINKWITZRILEYCROSSOVER_H_INCLUDED

/** This class implements a pair of filters to split an incoming signal into a low an high band. 
The filters used here are of Linkwitz/Riley type which ensures a flat magnitude response of the 
recombined (summed) output. However, the recombined output signal will not reconstruct the input 
signal exactly. Instead, the recombined signal will be an allpass-filtered version of the input. 
Linkwitz/Riley filters are made from a series connection of two identical Butterworth filters which
implies that the magnitude-response is that of a Butterworth filter squared. In particular, the 
gain at the cutoff/crossover frequency is -6.02 dB because the Butterworth filter has a gain of 
-3.01 dB there.

Stability: it has been stability-tested (with highest slope of 96 dB/oct) for crossover-frequencies 
down to 20 Hz with sample-rates up to 800 kHz - beyond that (either higher sample-rate or lower 
crossover-frequency (or both), filters may become numerically unstable). 

maybe rename to rsLinkwitzRileySplitter2

*/

template<class TSig, class TPar>
class rsLinkwitzRileyCrossOver
{

  typedef std::complex<TPar> Complex;       // preliminary
  //friend class rsCrossOver4Way<TSig, TPar>;

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. newMaxButterworthOrder should be even. */
  rsLinkwitzRileyCrossOver(int newMaxButterworthOrder = 8);

  /** Destructor. */
  ~rsLinkwitzRileyCrossOver() = default;

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample rate for this filter. */
  void setSampleRate(TPar newSampleRate);

  /** Switches the crossover on/off. When off, the input signal is passed through to the lowpass 
  output and the highpass output will be empty (zero). */
  void setActive(bool shouldBeActive) { active = shouldBeActive; }

  /** Sets the crossover frequency. */
  void setCrossoverFrequency(TPar newCrossoverFrequency);

  /** Sets the slope of the filters - this must be a multiple of 12 (\todo allow also for 
  6 dB/oct). */
  void setSlope(int newSlope);

  /** Sets the order of the butterworth-filters. */
  void setButterworthOrder(int newOrder);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the samplerate. */
  TPar getSampleRate() const { return sampleRate; }

  /** Informs, whether the crossover is active or not. */
  bool isActive() const { return active; }

  /** Returns the crossover frequency. */
  TPar getCrossoverFrequency() const { return crossoverFrequency; }

  /** Returns the slope of the filter pair. */
  int getSlope() const { return 12*butterworthOrder; }

  /** Fills the 'magnitudes' array with the magnitude response of the lowpass filter evaluated at 
  the frequencies passed in the'frequencies' array. Both arrays are assumed to be numBins long. */
  void getLowpassMagnitudeResponse(TPar* frequencies, TPar* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false) const;

  /** Fills the 'H' array with the complex freqiuency response of the lowpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getLowpassFrequencyResponse(TPar* frequencies, Complex* H, int numBins, 
    bool accumulate = false) const;

  /** Fills the 'magnitudes' array with the magnitude response of the highpass filter evaluated at
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getHighpassMagnitudeResponse(TPar* frequencies, TPar* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false) const;

  /** Fills the 'H' array with the complex freqiuency response of the highpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getHighpassFrequencyResponse(TPar* frequencies, Complex* H, int numBins, 
    bool accumulate = false) const;

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Calculates a pair of lowpass/highpass output samples from a (presumbaly) broadband input 
  signal. */
  inline void getSamplePair(TSig* in, TSig* outLow, TSig* outHigh)
  {
    if(!active)
    {
      *outLow  = *in;
      *outHigh = TSig(0);
      return;
    }

    TSig tmp  = *in;
    *outLow  = lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
    *outHigh = sumAllpass.getSampleDirect2(tmp) - *outLow;
  }

  /** Processes a buffer of samples. */
  inline void processBuffer(TSig* in, TSig* outLow, TSig* outHigh, int length)
  {
    if(!active)
    {
      memcpy(outLow, in, sizeof(TSig));
      memset(outHigh, 0, sizeof(TSig));
      return;
    }

    TSig tmp;
    for(int n = 0; n < length; n++)
    {
      tmp = in[n];
      outLow[n]  = lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
      outHigh[n] = sumAllpass.getSampleDirect2(tmp) - outLow[n];
    }
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the internal buffers of the filters to all zeros. */
  void resetBuffers();  // rename to reset

  //===============================================================================================

  rsBiquadCascade<TSig, TPar> sumAllpass; // temporarily move to public - use friend declaration

protected:

  /** Triggers a re-calculation of the filter coefficients. */
  void updateFilterCoefficients();

  // embedded objects:
  rsBiquadCascade<TSig, TPar> lowpass1, lowpass2; // direct-form performed worse than biquad-cascade


  TPar sampleRate;
  TPar crossoverFrequency;
  int  butterworthOrder;     // order of the Butterworth filters
  int  maxButterworthOrder;  // maximum order of the Butterworth filters
  bool active = true;

};

#endif

/**
\todo: write a class PerfectReconstructionCrossover that uses the following algorithm:

yL1 = LP1(x);   yH1 =       (x  -  yL1);  // x: input, yL1: 1st lowpass output, yH1: 1st highpass output, LP1: 1st lowpass filter
yL2 = LP2(yL1); yH2 = yH1 + (yL1 - yL2);  // yL2: 2nd lowpass output. yH2: 2nd highpass output
yL3 = LP3(yL2); yH3 = yH2 + (yL2 - yL3);
...etc.
the lowpass of the current stage is always applied the lowpass-output of the previous stage and the 
difference between the current lowpass-input and current lowpass-output is added to the 
highpass-output of the previous stage to form the highpass output of the current stage if we assume 
that x = yL1+yH1, then it follows that x = yL2+yH2 = yL3+yL3 = ...
this ensures that after any number of stages, the lowpass and highpass signals sum up to the 
original signal x however, the crossover frequency of the N-th order crossover will be different 
from the cutoff frequency of the individual 1st order stages -> we must adjust this -> find a 
cutoff scaler in terms of the total order N, we may also have to choose different cutoff 
frequencies for each stage -> this requires some research (find the 3-dB point of the N-th lowpass 
stage, scale the cutoff, etc...)

for an arbitrary N-th order crossover, the general algorithm may be stated as:
yL = lowpasses[0].getSample(x);
yH = x - yL;
for(int k = 1; k < N; k++)
{
  tmp = yL;
  yL  = lowpasses[k].getSample(yL);
  yH  = yH + (tmp - yL);
}
*/
