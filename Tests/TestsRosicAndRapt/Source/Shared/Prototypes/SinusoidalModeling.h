#pragma once

// todo: move file to library - maybe split into two files analyzer and synthesizer - or maybe 
// merge with the SinusoidalModel file


/** A class for synthesizing a sound from a sinusoidal model. */

template<class T>
class SinusoidalSynthesizer
{

public:


  /** \name Setup */

  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; }

  // setAmplitudeInterpolationMethod
  // setPhaseInterpolationMethod
  // etc.




  /** \name Inquiry */

  T getSampleRate() const { return sampleRate; }






  /** \name Processing */

  /** Synthesizes a sound from a sinusoidal model and returns it as std::vector. */
  std::vector<T> synthesize(const RAPT::rsSinusoidalModel<T>& model) const;

  /** Writes the given sinusoidal partial into the buffer x which is assumed to be of length 
  xLength. You may also pass an argument that time-shifts the partial's time-stamps with respect
  to their stored values. This is useful for synthesizing models that don't start at time zero. 
  The function is used internally by synthesize() and client code does not need to deal with it
  directly (but can, if needed). */
  void synthesizePartial(const RAPT::rsSinusoidalPartial<T>& partial, T* x, int xLength, 
    T timeShift = T(0)) const;


protected:


  T sampleRate = 44100;

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
  stable partial. */
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


  /** \name Processing */



  /** Analyzes the given input sound and returns the sinusoidal model object that models the given 
  sound. */
  RAPT::rsSinusoidalModel<T> analyze(T* sampleData, int numSamples, T sampleRate);
  // rename to analyzeSample

  /** Analyzes the given spectrogram and return the sinusoidal model for it */
  RAPT::rsSinusoidalModel<T> analyzeSpectrogram(
    const RAPT::rsMatrix<std::complex<T>>& spectrogram, T sampleRate);

  /** Creates and returns a complex spectrogram from the given sample data. Used internally by 
  analyze, so client code needs to call this directly only if it wants to see/plot/investigate the 
  underlying spectrogram analysis result. */
  RAPT::rsMatrix<std::complex<T>> getComplexSpectrogram(T* sampleData, int numSamples);

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
  of organizing the loops ...  Not yet implemented   */
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