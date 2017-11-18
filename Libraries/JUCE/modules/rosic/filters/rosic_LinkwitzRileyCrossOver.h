#ifndef rosic_LinkwitzRileyCrossOver_h
#define rosic_LinkwitzRileyCrossOver_h

namespace rosic
{

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
crossover-frequency (or both), filters may become numerically unstable). */

class rsLinkwitzRileyCrossOver
{

  friend class rsCrossOver4Way;

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. newMaxButterworthOrder should be even. */
  rsLinkwitzRileyCrossOver(int newMaxButterworthOrder = 8);

  /** Destructor. */
  ~rsLinkwitzRileyCrossOver();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample rate for this filter. */
  void setSampleRate(double newSampleRate);

  /** Sets the crossover frequency. */
  void setCrossoverFrequency(double newCrossoverFrequency);

  /** Sets the slope of the filters - this must be a multiple of 12 (\todo allow also for 
  6 dB/oct). */
  void setSlope(int newSlope);

  /** Sets the order of the butterworth-filters. */
  void setButterworthOrder(int newOrder);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the samplerate. */
  double getSampleRate() const
  {
    return sampleRate;
  }

  /** Returns the crossover frequency. */
  double getCrossoverFrequency() const
  {
    return crossoverFrequency;
  }

  /** Returns the slope of the filter pair. */
  int getSlope() const
  {
    return 12*butterworthOrder;
  }

  /** Fills the 'magnitudes' array with the magnitude response of the lowpass filter evaluated at 
  the frequencies passed in the'frequencies' array. Both arrays are assumed to be numBins long. */
  void getLowpassMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false);

  /** Fills the 'H' array with the complex freqiuency response of the lowpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getLowpassFrequencyResponse(double* frequencies, Complex* H, int numBins, 
    bool accumulate = false);

  /** Fills the 'magnitudes' array with the magnitude response of the highpass filter evaluated at
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getHighpassMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    bool inDecibels = false,
    bool accumulate = false);

  /** Fills the 'H' array with the complex freqiuency response of the highpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getHighpassFrequencyResponse(double* frequencies, Complex* H, int numBins, 
    bool accumulate = false);

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Calculates a pair of lowpass/highpass output samples from a (presumbaly) broadband input 
  signal. */
  INLINE void getSamplePair(double* in, double* outLow, double* outHigh);
  INLINE void getSamplePair(float*  in, float*  outLow, float*  outHigh);

  /** Processes a buffer of samples. */
  INLINE void processBuffer(double* in, double* outLow, double* outHigh, int length);
  INLINE void processBuffer(float*  in, float*  outLow, float*  outHigh, int length);

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the internal buffers of the filters to all zeros. */
  void resetBuffers();

  //===============================================================================================

protected:

  /** Triggers a re-calculation of the filter coefficients. */
  void updateFilterCoefficients();

  // embedded objects:
  rsBiquadCascade lowpass1, lowpass2; // direct-form performed worse than biquad-cascade
  rsBiquadCascade sumAllpass;

  double sampleRate;
  double crossoverFrequency;
  int    butterworthOrder;     // order of the Butterworth filters
  int    maxButterworthOrder;  // maximum order of the Butterworth filters

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsLinkwitzRileyCrossOver::getSamplePair(double* in, double* outLow, double* outHigh)
{
  double tmp  = *in;
  *outLow     = lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
  *outHigh    = sumAllpass.getSampleDirect2(tmp) - *outLow;
}

INLINE void rsLinkwitzRileyCrossOver::getSamplePair(float* in, float* outLow, float* outHigh)
{
  double tmp  = (double)*in;
  *outLow     = (float)lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
  *outHigh    = (float)(sumAllpass.getSampleDirect2(tmp) - *outLow);
}

INLINE void rsLinkwitzRileyCrossOver::processBuffer(double* in, double* outLow, double* outHigh, 
  int length)
{
  double tmp;
  for(int n=0; n<length; n++)
  {
    tmp        = in[n];
    outLow[n]  = lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
    outHigh[n] = sumAllpass.getSampleDirect2(tmp) - outLow[n];
  }
}

INLINE void rsLinkwitzRileyCrossOver::processBuffer(float* in, float* outLow, float* outHigh, 
  int length)
{
  float tmp;
  for(int n=0; n<length; n++)
  {
    tmp        = in[n];
    outLow[n]  = (float)lowpass2.getSampleDirect2(lowpass1.getSampleDirect2(tmp));
    outHigh[n] = (float)(sumAllpass.getSampleDirect2(tmp) - outLow[n]);
  }
}

//=================================================================================================

/** This class implements a stereo-version of the Linkwitz/Riley crossover. */

class rsLinkwitzRileyCrossOverStereo
{

  friend class rsCrossOver4Way;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  rsLinkwitzRileyCrossOverStereo(/*int newMaxButterworthOrder = 8*/)
  {
    active = true;
    mono   = false;
  }

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the sample rate for this filter. */
  void setSampleRate(double newSampleRate)
  {
    crossoverL.setSampleRate(newSampleRate);
    crossoverR.setSampleRate(newSampleRate);
  }

  /** Sets the crossover frequency. */
  void setCrossoverFrequency(double newCrossoverFrequency)
  {
    crossoverL.setCrossoverFrequency(newCrossoverFrequency);
    crossoverR.setCrossoverFrequency(newCrossoverFrequency);
  }

  /** Sets the slope of the filters - this must be a multiple of 12 (\todo allow also for 
  6 dB/oct). */
  void setSlope(int newSlope)
  {
    crossoverL.setSlope(newSlope);
    crossoverR.setSlope(newSlope);
  }

  /** Switches the crossover on/off. When off, the input signal is passed through to the lowpass 
  output and the highpass output will be empty (zero). */
  void setActive(bool shouldBeActive)
  {
    active = shouldBeActive;
  }

  /** Switches into mono-mode where only the left channel is calculated and the copied into the 
  right output channel. */
  void setMono(bool shouldBeMono)
  {
    mono = shouldBeMono;
  }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the samplerate. */
  double getSampleRate() const
  {
    return crossoverL.getSampleRate();
  }

  /** Informs, whether the crossover is active or not. */
  bool isActive() const
  {
    return active; // todo: drag this flag into class LinkwitzReileyCrossover
  }

  /** Returns the crossover frequency. */
  double getCrossoverFrequency() const
  {
    return crossoverL.getCrossoverFrequency();
  }

  /** Returns the slope of the filter pair. */
  int getSlope() const
  {
    return crossoverL.getSlope();
  }

  /** Fills the 'magnitudes' array with the magnitude response of the lowpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getLowpassMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false)
  {
    crossoverL.getLowpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, 
      accumulate);
  }

  /** Fills the 'magnitudes' array with the magnitude response of the highpass filter evaluated at 
  the frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getHighpassMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    bool inDecibels = false, bool accumulate = false)
  {
    crossoverL.getHighpassMagnitudeResponse(frequencies, magnitudes, numBins, inDecibels, 
      accumulate);
  }

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a sample-frame. */
  INLINE void getSampleFrame(double *inL, double *inR, double *outLowL, double *outLowR, 
    double *outHighL, double *outHighR);
  INLINE void getSampleFrame(float  *inL, float  *inR, float  *outLowL, float  *outLowR, 
    float  *outHighL, float  *outHighR);

  /** Processes a buffer of samples. */
  INLINE void processBuffer(double *inL, double *inR, double *outLowL, double *outLowR, 
    double *outHighL, double *outHighR, int length);
  INLINE void processBuffer(float  *inL, float  *inR, float  *outLowL, float  *outLowR, 
    float  *outHighL, float  *outHighR, int length);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers of the filters to all zeros. */
  void resetBuffers()
  {
    crossoverL.resetBuffers();
    crossoverR.resetBuffers();
  }

  //===============================================================================================

protected:

  rsLinkwitzRileyCrossOver crossoverL, crossoverR;

  bool active, mono;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsLinkwitzRileyCrossOverStereo::getSampleFrame(double *inL, double *inR, 
  double *outLowL, double *outLowR, double *outHighL, double *outHighR)
{
  if(!active)
  {
    *outLowL  = *inL;
    *outLowR  = *inR;
    *outHighL = 0.0;
    *outHighR = 0.0;
  }

  crossoverL.getSamplePair(inL, outLowL, outHighL);
  if(mono)
  {
    *outLowR  = *outLowL;
    *outHighR = *outHighL;
  }
  else
    crossoverR.getSamplePair(inR, outLowR, outHighR);
}

INLINE void rsLinkwitzRileyCrossOverStereo::getSampleFrame(float *inL, float *inR, float *outLowL, 
  float *outLowR, float *outHighL, float *outHighR)
{
  if(!active)
  {
    *outLowL  = *inL;
    *outLowR  = *inR;
    *outHighL = 0.0;
    *outHighR = 0.0;
  }

  crossoverL.getSamplePair(inL, outLowL, outHighL);
  if(mono)
  {
    *outLowR  = *outLowL;
    *outHighR = *outHighL;
  }
  else
    crossoverR.getSamplePair(inR, outLowR, outHighR);
}

INLINE void rsLinkwitzRileyCrossOverStereo::processBuffer(double *inL, double *inR, double *outLowL, 
  double *outLowR, double *outHighL, double *outHighR, int length)
{
  if(!active)
  {
    memcpy(outLowL, inL, sizeof(double));
    memcpy(outLowR, inR, sizeof(double));
    memset(outHighL, 0, sizeof(double));
    memset(outHighR, 0, sizeof(double));
  }

  crossoverL.processBuffer(inL, outLowL, outHighL, length);
  if(mono)
  {
    memcpy(outLowR, outLowL, sizeof(double));
    memcpy(outHighR, outHighL, sizeof(double));
  }
  else
    crossoverR.processBuffer(inR, outLowR, outHighR, length);
}

INLINE void rsLinkwitzRileyCrossOverStereo::processBuffer(float *inL, float *inR, float *outLowL, 
  float *outLowR, float *outHighL, float *outHighR, int length)
{
  if(!active)
  {
    memcpy(outLowL, inL, sizeof(float));
    memcpy(outLowR, inR, sizeof(float));
    memset(outHighL, 0, sizeof(float));
    memset(outHighR, 0, sizeof(float));
  }

  crossoverL.processBuffer(inL, outLowL, outHighL, length);
  if(mono)
  {
    memcpy(outLowR, outLowL, sizeof(float));
    memcpy(outHighR, outHighL, sizeof(float));
  }
  else
    crossoverR.processBuffer(inR, outLowR, outHighR, length);
}

}

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
