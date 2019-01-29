#pragma once

// todo: move file to library - maybe split into two files analyzer and synthesizer - or maybe 
// merge with the SinusoidalModel file


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
generating main-lobes, iFFT, overlap/add as in Xavier Serra's SMS framework */

template<class T>
class SinusoidalSynthesizer
{

public:

  /** Enumeration of the available phase interpolation methods. */
  enum class PhaseInterpolationMethod
  {
    cubicHermite = 0,
    quinticHermite,
    tweakedFreqIntegral
  };


  /** \name Setup */

  /** Sets the synthesis sample rate that will be used in synthesize(). */
  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; }

  /** Switches the amplitude interpolation method between cubic (true) and linear (false). By 
  default, the amplitude is interpolated linearly between datapoints, mainly to avoid 
  overshooting artifacts, when the model has "fade-in/out" datapoints that are very close to the
  actual data points. When this is not the case, using cubic interpolation may give a "rounder"
  sounding result. */
  void setCubicAmplitudeInterpolation(bool shouldBeCubic) { cubicAmplitude = shouldBeCubic; }

  /** Switches the phase interpolation method between cubic (true) and linear (false). By default,
  the phase is interpolated cubically between datapoints because linear interpolation may lead to
  audible discontinuities in the derivative of the resulting sound (or is it the 2nd derivative? 
  check this!). */
  //void setCubicPhaseInterpolation(bool shouldBeCubic) { cubicPhase = shouldBeCubic; }


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


  /** Given an array of time-stamps and corresponing frequency and wrapped phase values, this 
  function computes the corresponding array of unwrapped phase values by numerically integrating
  the frequency array and then re-adjusting the resulting (unwrapped) phases to have a value that 
  is consistent with the wrappedPhase values (i.e. a suitable multiple of 2*pi shifted from the 
  stored value). To integrate the frequecy data, we also need the time axis because the data may
  be nonuniformly sampled. */
  std::vector<T> unwrapPhase(const std::vector<T>& time, 
    const std::vector<T>& freq, const std::vector<T>& wrappedPhase) const; 





protected:

  T sampleRate = 44100;

  bool cubicAmplitude = false;

  PhaseInterpolationMethod phaseInterpolation = PhaseInterpolationMethod::tweakedFreqIntegral;
  // later use quinticHermite by default


  // maybe get rid of these:
  //bool cubicPhase = true;
  //bool accumulatePhaseDeltas = true;
  // Notes: accumulation gives rise to an O(M^2) complexity of the phase unwrapping algorithm where 
  // M is the number of data points in the respective partial - but without accumulation there may 
  // be phase errors by 180° at certain datapoints after the unwrapping - figure out what's going 
  // on, maybe change the unwrapping algorithm completely or rather let the user select between 
  // various algorithms

};

//=================================================================================================

/** Creates a sinusoidal model of an input sound by means of identifying and tracking stable 
sinusoidal partials in its spectrogram. */

template<class T>
class SinusoidalAnalyzer
{

public:


  /** \name High level setup */


  // Good choices for sinusoidal analysis with a minimum resolvable frequency delta df:
  // Window: Blackman-Harris
  // ...
  // todo: make a function setParameterPreset(int index, double df, double fs)
  // maybe it should slo include setting the maxFreqDeltaBase and maxFreqDeltaSlope variables

  /** Sets the fade-in time in seconds for newborn partial tracks. Whenever a new track is born, we 
  prepend a datapoint with the same frequency, zero amplitude, a time-stamp given by the time of 
  the new track-start minus the fade-in time and a phase set to an appropriate value according to 
  frequency and fade-time. If you set this to zero, the additional fade-in datapoint will be left 
  out. */
  void setFadeInTime(T newTime) { fadeInTime = newTime; }

  /** Sets the fade-out time in seconds for partial tracks that have died. Whenever a track dies, 
  we append an additional datapoint with the same frequency, zero amplitude, time-stamp given by 
  the time of death plus the fade-out time and a phase set to an appropriate value according to 
  frequency and fade-time. If you set this to zero, the additional fade-out datapoint will be left 
  out.  */
  void setFadeOutTime(T newTime) { fadeOutTime = newTime; }

  /** Sets the minimum length (in seconds) that a track must have in order to be considered a 
  stable partial. If you have set up some fade-in and fade-out times, this minimum length should
  be at least the sum of these two fade times (because each track will be guaranteed to be at least
  faedInTime+fadeOutTime seconds long - so a value less than that for minTrackLength will 
  effectively turn the threshold off (well, of course, if you want to turn it off, you can do it 
  like that). */
  void setMinimumTrackLength(T newLength) { minTrackLength = newLength; }



  /** \name Low level setup */

  /** Sets the FFT size for the underlying spectrogram analysis. Should be >= the block size. */
  inline void setTrafoSize(int newSize)      { sp.setTrafoSize(newSize); }

  /** Sets the block size for the underlying spectrogram analysis. Should be <= the FFT size. */
  inline void setBlockSize(int newBlockSize) { sp.setBlockSize(newBlockSize); }

  /** Sets FFT size and block size for the underlying spectrogram at the same time (setting them 
  simultaneously avoids potential temporary violations of trafoSize >= blockSize (which triggers an
  assertion) during setup). */
  inline void setBlockAndTrafoSize(int newBlockSize, int newTrafoSize)
  {
    sp.setBlockAndTrafoSize(newBlockSize, newTrafoSize);
  }
  
  /** Sets the hop size for the underlying spectrogram analysis. Should typically be some fraction 
  of the block size (such as 1/2, 1/4 or something). */
  inline void setHopSize(int newHopSize)     { sp.setHopSize(newHopSize); }

  /** Sets the analysis window type for the underlying spectrogram analysis. Should be one of the 
  type in RAPT::rsWindowFunction::windowTypes. The window type affects the time-frequency tradeoff
  and the precision of the partial frequency estimation.... */
  inline void setWindowType(int newType)     { sp.setAnalysisWindowType(newType); }

  // void setContinuationAlgorithm

  /** Sets the maximum allowed frequency difference that a sinusoidal track may have between 
  successive frames in order to be considered a stable sinusoid. This maximum allowed difference
  may increase towards higher frequencies, so this here is the baseline value that is used at zero
  frequency. At higher frequencies an higher threshold may be used according to 
  setMaxFreqDeltaSlope.  */
  void setMaxFreqDeltaBase(T newDifference) { maxFreqDeltaBase = newDifference; }

  /** Sets the amount by which the maximum allowed frequency difference increases with 
  frequency. */
  void setMaxFreqDeltaSlope(T newSlope) { maxFreqDeltaSlope = newSlope; }

  /** Sets the magnitude threshold, above which a spectral peak must be, to be considered a 
  potential sinusoidal track. The threshold is in decibels and relative to the global maximum of
  the spectrum. A value of -40 means that peaks that are 40 dB below the global maximum are not 
  considered. The value should typically be less than the sidelobe level of the analysis window in
  order to not pick up on sidelobes. */
  void setRelativeLevelThreshold(T newThreshold) { magThreshold = rsDbToAmp(newThreshold); }



  /** \name Inquiry */

  // T getRequiredWindowSize(int windowType, T freqDelta, T sampleRate);

  /** Returns the maximum allowed frequency difference that a sinusoidal track may have between 
  successive frames in order to be considered a stable sinusoid. This maximum allowed difference
  may different at different frequencies, so you should pass a reference frequency. Mainly for 
  internal use. - maybe we should distinguis between a function that returns the raw parameter and 
  one that returns the final computed value  */
  T getMaxFreqDelta(T referenceFreq) const 
  { 
    return maxFreqDeltaBase + maxFreqDeltaSlope*referenceFreq; 
  }
  // todo: check SMS tool which formula is used there (this here is my ad-hoc formula)
  // i think, the value should depend on (be proportional to) the hop-size

  /** Returns the analysis hop size. */
  int getHopSize() const { return sp.getHopSize(); }

  /** Returns the minimum block size that is required to obtain a given frequency resolution at a 
  given sample rate. This resolution gives the minimum frequency difference by which two partials 
  have to be separated in order to be resolved as two distinct partials. If their frequencies are 
  closer together than that, they may appear as a single partial to the analyzer. Because this
  resolvability depends on the mainlobe width of the chosen analysis window, you must also tell 
  the function, which type of window you want to use. */
  static int getRequiredBlockSize(int windowType, T frequencyResolution, T sampleRate, 
    bool oddSize = false);

  /** Returns minimum possible amplitude threshold that can reasonably be used for a given window 
  type without mistakenly picking up the sidelobes of the window as partials. With margin = 0, this 
  is actually just the sidelobe level of the window function itself. A margin of ...dB proved to
  be reasonable in practice */
  static T getRequiredThreshold(int windowType, T dBmargin);
  // should return the sidelobe level of given window in dB 
  // todo: maybe provide a deafult margin -> figure out reasonable value empirically


  /** \name Processing */

  /** Analyzes the given input sound and returns the sinusoidal model object that models the given 
  sound. */
  RAPT::rsSinusoidalModel<T> analyze(T* sampleData, int numSamples, T sampleRate);
  // rename to analyzeSample, maybe don't let the user pass the sampleRate, instead have a 
  // setSampleRate function - this would be consistent with the synthesizer

  /** Analyzes the given spectrogram and return the sinusoidal model for it */
  RAPT::rsSinusoidalModel<T> analyzeSpectrogram(
    const RAPT::rsMatrix<std::complex<T>>& spectrogram, T sampleRate);

  /** Creates and returns a complex spectrogram from the given sample data. Used internally by 
  analyze, so client code needs to call this directly only if it wants to see/plot/investigate the 
  underlying spectrogram analysis result. */
  RAPT::rsMatrix<std::complex<T>> getComplexSpectrogram(T* sampleData, int numSamples);

  /** Returns an array of indices of peaks in the given array x of length N. A peak at index i is 
  defined by the condition x[i-1] < x[i] < x[i+1]. You can also pass a relative threshold (with 
  respect to the highest peak) below which peaks will not be reported. For example, if your 
  maximum peak is 5 and threshToMax = 0.1, then all peaks above 0.5 will be reported. If you pass 
  zero, it effectively turns the thresholding off. */
  static std::vector<int> peakIndices(T* x, int N, T threshToMax = T(0));
  // todo: that function seems to be more generally useful - maybe move out of the class (and maybe 
  // into rsArray)

  /** Given an array of spectral magnitude values x and a peak-index k such that 
  x[k-1] < x[k] < x[k+1], this function computes the exact location of the peak 
  pos = k + d, where d is some number between -1..+1 and the actual y-value of the peak at that 
  exact location by fitting a quadratic parabola to the dB-values around x[k]. */
  static void spectralMaximumPositionAndValue(T *x, int k, T* pos, T* val);

protected:

  /** For a given frequency (belonging to a spectral peak of the current frame), search through the
  tracks array to find the one that is closest in frequency, if it is within a given range around 
  the peak-frequency and for which the trackContinued flag is false (i.e. the track has not been 
  used up by another current peak frequency already - each track can be continued only with one
  partner). Called from continuePartialTracks */
  size_t findBestMatchingTrack(T frequency, std::vector<RAPT::rsSinusoidalPartial<double>>& tracks, 
    const std::vector<bool>& trackContinued) const;

  /** This function implements the peak continuation step - for all current spectral peaks in 
  newPeakData, find a corresponding continuation partner among the activeTracks - 3 situations have 
  to be handled:
  -when a partner is found, continue the track, i.e. append the newPeakData to the corresponding
   track in activeTracks
  -when no partner is found, create a new track ("birth"), i.e. start a new track in activeTracks
  -all active tracks that have not been used in this continuation are killed (i.e. moved to 
   finishedTracks */
  void continuePartialTracks1(
    std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeakData,
    std::vector<RAPT::rsSinusoidalPartial<T>>& activeTracks,
    std::vector<RAPT::rsSinusoidalPartial<T>>& finishedTracks) const;
  // rename to continuePartialTracks1 and let continuePartialTracks be a dispatcher that selects
  // between continuePartialTracks1/continuePartialTracks2
  // or rename to findContinuations

  /** For a given frequency, search through the peaks array to find one that is closest in 
  frequency, subject to the constraint that is within a given maximum deviation and that the 
  peakUnused flag is false. This is the dual function to findBestMatchingTrack for the other way
  of organizing the loops */
  size_t findBestMatchingPeak(T frequency, std::vector<RAPT::rsInstantaneousSineParams<T>>& peaks, 
    const std::vector<bool>& peakUsed) const;

  /** Alternative version of the peak-tracking algoritm. This one loops over all the tracks to find 
  the best match in newPeakData instead of looping over all peaks to find a best match in the 
  activeTracks, i.e. the roles are reversed. This is, how it's described in the literature */
  void continuePartialTracks0(
    std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeakData,
    std::vector<RAPT::rsSinusoidalPartial<T>>& activeTracks,
    std::vector<RAPT::rsSinusoidalPartial<T>>& finishedTracks) const;

  /** Internal function called from continuePartialTracks1/2... */
  void applyContinuations(
    std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeaks,
    std::vector<RAPT::rsSinusoidalPartial<T>>& aliveTracks,
    std::vector<RAPT::rsSinusoidalPartial<T>>& deadTracks,
    std::vector<size_t>& births, std::vector<size_t>& deaths,
    std::vector<std::pair<size_t, size_t>>& continuations) const;

  /** Applies some post-processing to the model to clean it up such as deleting spurious tracks, 
  merging tracks that seem to represent a single partial, etc. */
  void cleanUpModel(RAPT::rsSinusoidalModel<T>& model) const;




  // fade-in and fade-out times for partials (in seconds):
  T fadeInTime  = 0.01;
  T fadeOutTime = 0.01;
  T minTrackLength = 0;

  // frequency difference threshold parameters:
  T maxFreqDeltaBase  = 100;
  T maxFreqDeltaSlope = 0.01;

  T magThreshold = 0.0;

  int contAlgo = 0;  // peak continuation algorithm


  RAPT::rsSpectrogram<T> sp;   // spectrogram processor

};