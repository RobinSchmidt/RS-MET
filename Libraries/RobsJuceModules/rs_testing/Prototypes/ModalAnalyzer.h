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



protected:

  T sampleRate = T(1);
  int maxNumModes = 1024;

  // ToDo:
  // -Have parameters that control the sensitivity of peak-detection/picking algorithms

  // T skip

};