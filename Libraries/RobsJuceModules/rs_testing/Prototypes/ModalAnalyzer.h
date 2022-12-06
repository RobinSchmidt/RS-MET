#pragma once

/** Data structure to hold parameters of a single modal filter (with attack/decay envelope).
todo: maybe later allow for more general modal models with, for example, with two-stage decay,
etc. */
template<class T>
struct rsModalFilterParameters
{
  T freq = 1000, amp = 1, att = 0.01, dec = 0.1, phase = 0;
};
// ToDo: document the units of the parameters - especially for the phase. I think, it's in degrees.

std::vector<double> synthesizeModal(
  const rsModalFilterParameters<double>& params, double sampleRate, int length);

std::vector<double> synthesizeModal(
  const std::vector<rsModalFilterParameters<double>>& params, double sampleRate, int length);
// move to rapt or rosic



//=================================================================================================

/** A class for creating a modal model of a sound. The input is a sinusoidal model. */

template<class T>
class rsModalAnalyzer  // maybe rename into rsModalAnalyzerFromSineModel
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  // setPhaseEstimationMethod(newMethod);
  // setFrequencyEstimationMethod(newMethod); // mean, median, etc.

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Takes a sinusoidal model as input and creates the set/array of modal filter parameters that 
  can be used to approximate the given model. */
  std::vector<rsModalFilterParameters<T>> getModalModel(const RAPT::rsSinusoidalModel<T>& model);

  /** Takes a sinusoidal partial as input and creates the modal filter parameters that can be used 
  to approximate the given partial. */
  rsModalFilterParameters<T> getModalModel(const RAPT::rsSinusoidalPartial<T>& partial);

  T estimatePhaseAt(const RAPT::rsSinusoidalPartial<T>& partial,
    int dataPointIndex, T frequency, T timeInstant = T(0));

  T estimateFrequency(const RAPT::rsSinusoidalPartial<T>& partial, int startIndex, int endIndex);




protected:

};

//=================================================================================================

/** Another modal analyzer that operates directly on the time domain signal. */

template<class T>
class rsModalAnalyzer2
{


public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample-rate. This determines the values of the frequencies that will be written
  into the model. */
  void setSampleRate(T newRate) { sampleRate = newRate; }

  /** Sets the frequency bandwidth for the peak-masking...tbc... */
  void setMaskWidth(T newWidth) { maskWidth = newWidth; }
  // Intention: If we have a mode at 1000 Hz at amplitude 1 and maskWidth is 100, then modes at 
  // 900 = 1000-100 or 1100 = 1000+100 Hz must be higher than 0.5 in amplitude in order to be 
  // also detected, i.e. to have their peaks taken seriously and not be masked by the higher peak
  // at 1000. -> check if it actually really works that way in a unit-test, then add that to the 
  // documentation.


  /** Determines at which time instant the phase of the resynthesized mode should match the phase
  in the extracted mode. If true is passed, we will match the phase at a zero-crossing that is near
  the peak of the amp-envelope. If false is passed, we'll match it at the very first zero-crossing
  that is encountered. The former option may be better for reducing the overall maximum error in 
  the residual where the latter may be better to more accurately model the attack transients. I'm 
  not yet sure what's the better strategy, so both options are available. Some experimentation is 
  needed....tbc... */
  void setPhaseMatchNearPeak(bool shouldMatchNearPeak) { matchPhaseNearPeak=shouldMatchNearPeak; }



  // ToDo: setLevelThreshold, setMaxNumModes, setDecayMeasurementAmplitudes

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Takes a time domain signal as input and creates the set/array of modal filter parameters that 
  can be used to approximate the given model. */
  std::vector<rsModalFilterParameters<T>> analyze(T* sampleData, int numSamples);
    // has same name, signature and semantics as rsHarmonicAnalyzer::analyze

  /** Given an imput signal x (typically the full, original signal) of length N, this function 
  extract a single mode and writes the result into y (which may point to tze same aray as x). It 
  does this by applying a bidirectional bandpass filter with given center frequency and absolute 
  bandwidth in Hz. The optional wrk array (of length N_wrk) can be used to allow the bidirectional
  filter have a ringout/warmup phase between forward and backward pass. */
  void extractMode(const T* x, T* y, int N, T centerFreqHz, T bandwidthHz, 
    T* wrk = nullptr, int N_wrk = 0);
  // ToDo: Try to get rid of the ringout/warmup buffer by solving the state initialization problem
  // analyitically. I've already done something similar for a general first order filter. This here
  // requires to figure out similar equations for the second order case. Perhaps, it can be done by
  // breaking it down into a parallel connection of two complex 1st order filters via partial 
  // fraction expansion and then using the 1st order solution - we'll see....

  void extractModeEnvelope(const T* x, T* y, int N);


protected:

  // User parameters:
  int maxNumModes = 1024;
  T sampleRate = T(1);
  T threshDb   = T(-60);   // Relative threshold for modes to be taken seriously
  T maskWidth  = T(10);    // Half-bandwidth of the peak-masks in pre-analysis
  T decayAmp1  = T(0.5);
  T decayAmp2  = T(0.25);
  bool matchPhaseNearPeak = true;

  // ToDo:
  // -Have parameters to tune the bandpass filters used in the mode-extraction step. Let the user
  //  select the bandwidth (maybe as a multiplier for the preliminary bandwidth estimate coming from
  //  the FFT pre-analysis) and the order. We may want to use higher order (Gaussian) filters later.
  //  That's a job for EngineersFilter - currently, we use a simple biquad SVF.


  // Temporary data buffers and embedded objects:
  std::vector<T> buf1, buf2;
  std::vector<T> peakPositions, peakHeights;
  std::vector<rsVector2D<T>> peaks;
  rsFourierTransformerRadix2<T> ft;

};