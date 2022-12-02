#pragma once

/** Data structure to hold parameters of a single modal filter (with attack/decay envelope).
todo: maybe later allow for more general modal models with, for example, with two-stage decay,
etc. */
template<class T>
struct rsModalFilterParameters
{
  T freq = 1000, amp = 1, att = 0.01, dec = 0.1, phase = 0;
};

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

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample-rate. This determines the values of the frequencies that will be written
  into the model. */
  void setSampleRate(T newRate) { sampleRate = newRate; }


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


protected:




  // User parameters:
  T sampleRate = T(1);
  int maxNumModes = 1024;


  // Internal data and objects:
  std::vector<T> buf1, buf2;
  std::vector<T> peakPositions, peakHeights;
  std::vector<rsVector2D<T>> peaks;
  rsFourierTransformerRadix2<T> ft;





  // ToDo:
  // -Have parameters that control the sensitivity of peak-detection/picking algorithms

  // T skip

};