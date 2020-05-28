#pragma once


/** A class for synthesizing a sound from a sinusoidal model via interpolating the instantaneous 
phases and amplitudes of the partials up to the sample-rate and driving an oscillator bank with 
these interpolated values. The (unwrapped) instantaneous phases are obtained by numerically 
integrating the instantaneous frequency values and then re-adjusting the results according to the 
stored phase values. Then, these instantaneous phase values (which are still at datapoint-rate) are 
interpolated up to the sample-rate. A similar but simpler procedure is done to the instantaneous 
amplitudes (there's no integration or unwrapping involved here). 

\todo:
This synthesis algorithm here is conceptually straightforward (and perhaps the most accurate? 
maybe - verify) but also quite expensive - implement a more efficient synthesis algorithm by 
generating main-lobes, iFFT, overlap/add as in Xavier Serra's SMS framework 

\todo:
Maybe include an option to not synthesize partials above the Nyquist frequency ...but what when a
partial goes above the limit only for a short moment?

*/

template<class T>
class rsSinusoidalSynthesizer
{

public:

  /** Enumeration of the available phase interpolation methods. cubicHermite is the standard 
  method that is described in the literature - it uses the instantaneous phase data to fix the 
  phase-values at the datapoints and the instantaneous frequency data to fix the phase derivative 
  at the datapoints. tweakedFreqIntegral obtains preliminary (unwrapped) instantaneous phases at 
  the datapoints by computing a numerical integral of the instantaneous frequencies (via the 
  trapezoidal method) and then tweaks the results to make them consistent with the stored 
  instantaneous phase values (chooses the closest phase value that is consistent with the stored 
  value). Then, on this unwrapped phase data, it uses natural cubic spline interpolation.
  quinticHermite is supposed to use quintic instead of cubic interpolation - but doesn't work 
  yet */
  enum class PhaseInterpolationMethod
  {
    cubicHermite = 0,    // this is the default
    tweakedFreqIntegral,  
    quinticHermite       // does not seem to work yet 
  };


  /** \name Setup */

  /** Sets the synthesis sample rate that will be used in synthesize(). */
  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; }

  /** Switches the amplitude interpolation method between cubic (true) and linear (false). By 
  default, the amplitude is interpolated linearly between datapoints, mainly to avoid 
  overshooting artifacts, when the model has "fade-in/out" datapoints that are very close to the
  actual data points. When this is not the case, using cubic interpolation may give a "rounder"
  sounding result. However, it's perhaps better to use linear interpolation and a smoothing
  filter (todo: add optional smoothing facilities) */
  void setCubicAmplitudeInterpolation(bool shouldBeCubic) { cubicAmplitude = shouldBeCubic; }
  // todo: maybe optionally interpolate the amplitude linearly on a dB scale - and maybe make it
  // switchable separately for attack and decay sections - for attacks, linear may be more 
  // suitable but for decay log-linear (i.e. exponential) may be better (for attack, actually 
  // inverted exponential may be best - that's what analog envelopes do)

  /** Switches the phase interpolation method between cubic (true) and linear (false). By default,
  the phase is interpolated cubically between datapoints because linear interpolation may lead to
  audible discontinuities in the derivative of the resulting sound (or is it the 2nd derivative? 
  check this!). */
  //void setCubicPhaseInterpolation(bool shouldBeCubic) { cubicPhase = shouldBeCubic; }

  void setPhaseInterpolation(PhaseInterpolationMethod method) 
  { 
    phaseInterpolation = method; 
  }



  /** \name Inquiry */

  /** Returns the sample-rate at which this synthesizer will synthesize a sound when asked to do 
  so by calling synthesize(). */
  T getSampleRate() const { return sampleRate; }


  /** \name Processing */

  /** Synthesizes a sound from a sinusoidal model and returns it as std::vector. The sound will be
  synthesized at a sample-rate that you can set up via setSampleRate. */
  std::vector<T> synthesize(const RAPT::rsSinusoidalModel<T>& model) const;

  /** Writes the given sinusoidal partial into the buffer x which is assumed to be of length 
  xLength. You may also pass an argument that time-shifts the partial's time-stamps with respect
  to their stored values. This is useful for synthesizing models that don't start at time zero. 
  The function is used internally by synthesize() and client code does not need to deal with it
  directly (but can, if needed). */
  void synthesizePartial(const RAPT::rsSinusoidalPartial<T>& partial, T* x, int xLength, 
    T timeShift = T(0)) const;


  std::vector<T> getInterpolatedAmplitudes(const RAPT::rsSinusoidalPartial<T>& partial, 
    const std::vector<T>& timeAxisCoarse, const std::vector<T>& timeAxisFine) const;

  std::vector<T> getInterpolatedPhases(const RAPT::rsSinusoidalPartial<T>& partial, 
    const std::vector<T>& timeAxisCoarse, const std::vector<T>& timeAxisFine) const;

  std::vector<T> phasesViaTweakedIntegral(const RAPT::rsSinusoidalPartial<T>& partial,
    const std::vector<T>& timeAxisCoarse, const std::vector<T>& timeAxisFine) const;

  std::vector<T> phasesHermite(const RAPT::rsSinusoidalPartial<T>& partial,
    const std::vector<T>& timeAxisCoarse, const std::vector<T>& timeAxisFine, 
    bool quintic) const;



  /** Given an array of time-stamps and corresponing frequency and wrapped phase values, this 
  function computes the corresponding array of unwrapped phase values by numerically integrating
  the frequency array and then re-adjusting the resulting (unwrapped) phases to have a value that 
  is consistent with the wrappedPhase values (i.e. a suitable multiple of 2*pi shifted from the 
  stored value). To integrate the frequecy data, we also need the time axis because the data may
  be nonuniformly sampled. */
  //std::vector<T> unwrapPhase(const std::vector<T>& time, 
  //  const std::vector<T>& freq, const std::vector<T>& wrappedPhase) const; 





protected:


  T sampleRate = 44100;

  bool cubicAmplitude = false;

  //PhaseInterpolationMethod phaseInterpolation = PhaseInterpolationMethod::tweakedFreqIntegral;
  PhaseInterpolationMethod phaseInterpolation = PhaseInterpolationMethod::cubicHermite;
  // later use quinticHermite by default ..or maybe make tests, which is better


  // maybe get rid of these:
  //bool cubicPhase = true;
  //bool accumulatePhaseDeltas = true;
  // Notes: accumulation gives rise to an O(M^2) complexity of the phase unwrapping algorithm where 
  // M is the number of data points in the respective partial - but without accumulation there may 
  // be phase errors by 180° at certain datapoints after the unwrapping - figure out what's going 
  // on, maybe change the unwrapping algorithm completely or rather let the user select between 
  // various algorithms

};
