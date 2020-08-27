#pragma once

/** A class to create a sinusoidal model for quasi-harmonic signals. The assumed harmonic 
relationship between the partials allows for a different analysis algorithm that is tailored 
specifically for harmonic signals and is supposed to give better results than a more general
sinusoidal analysis algorithm that doesn't make such an assumption (if the assumption indeed holds
true for the analyzed signal, of course). The algorithm works as follows:

-pre-process audio (pitch flattening):
 -obtain cycle marks
 -obtain a time warping map that makes the cycle-marks equidistant (distance is chosen to be a 
  power of two >= the longest cycle in the input signal)
 -time-warp the signal -> all cycles have the same, FFT-friendly length
-harmonic analysis:
 -analyze each cycle of pitch flattened signal by an FFT
 -the FFT size is equal to the cycle-length -> frequencies don't need to be 
  estimated, they are known in advance
 -only amplitude and phase have to be measured (i.e. simply read off from the FFT data)
 -update: now the FFT size is equal to some power-of-2 multiple of the cycle-length, if it's > 1,
  then again partial frequecies must be estimated (but we know roughly where to look for them)
-post-process model data (account for pitch flattening):
 -move time instants of datapoints according to the inverse time-warping map
 -modify frequencies according to the applied stretch factors

...this is a bit out of date - this was the first version of the algo - it has been refined 
since - but maybe we should actually keep that simpler and easier to understand algo for 
reference...


The so obtained model models any inharmonicity, transients and noise in the input signal as fast 
variations of instantaneous frequencies and amplitudes of the partials - when the envelopes
are lowpass-filtered before resynthesis, we can resynthesize only the quasi-harmonic part of the 
sound. Time domain subtraction from the original should give a residual containing inharmonic, 
noisy and transient parts. The model also includes a DC component which models - apparently - any
DC present in the signal, but if the amplitude of the DC varies, it may also model subharmonic
content. However, subharmonic content may also be partially modeled/encoded by amplitude- and 
frequency modulation of the harmonics.

todo: adjust the phase of the initial cycle - the analysis computes a phase value at the sample 
where the window sample happens to fall on - but if that window center is before time 0 in the 1st
block...wait...the initial datapoint may actually get a time-index < 0 if less than a half-cycle 
is in the initial section...figure out the initial/final

\todo: give the user the option to prepend and append some zeros before the analysis to avoid 
edge-artifacts - this padding, if used, has to be accounted for in the model after analysis is 
finished. the padding should be at least one cycle long, maybe use two to be on the safe side

-maybe have an "oversampling" option: let the hop-size be a half or quarter of the block-size
 -samples amplitudes and phases more densely in time

-make a sample editor module for ToolChain:
 -let the user choose background: blank (some color according to colorscheme), spectrogram, 
  phasogram, etc.. and forground: blank, waveform, partial freqs, etc...
-the user can make edits based on a sinusoidal model of the sound, spectrogram and waveform (and 
 later maybe more representations - what about a wavelet trafo, for example? or a single big 
 spectrum for the whole sound)
-it should be possible to save and load sinusoidal models (devise a file-format)

*/

template<class T>
class rsHarmonicAnalyzer
{

public:

  /*
  enum class Preset
  {
    hamming_4_4,
    blackman_4_4
  };
  */


  rsHarmonicAnalyzer();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample-rate. This determines the values of the frequencies that will be written
  into the model. */
  void setSampleRate(T newRate) 
  { 
    sampleRate = newRate; 
    cycleFinder.setSampleRate(sampleRate);
  }

  //void setRemoveDC(bool shouldBeRemoved);

  /** Sets the length of the sinc interpolator that is used for the pitch flattening. In tests with
  a sum-of-sines input, it turned out that the amplitude of the residual seems to be roughly 
  inverserly proportional to that length. But that's only a rough tendency - for example, there's a 
  huge difference (factor ~3) between 16 and 17 ...more research needed */
  void setSincInterpolationLength(T newLength) { sincLength = newLength; }
  // rename to setSincLength

  /** Sets the amount of oversampling the FFT spectrum, i.e. factor by which the FFT buffer is
  lengthened by zero-padding. Must be a power of two. Zero padding (i.e. spectral oversampling) may 
  be useful when the partials are (slightly) inharmonic to help detect the actual frequency of the 
  partial and/or suppress artifacts from spectral leakage of the inharmonic into other bins (due
  to non-rectangular windowing in case of zero-padding...to be implemented...unpadded FFTs are 
  always be computed using the rectangular window...hmm..but maybe we don't have to?)
  (for totally inharmonic partials, an entirely
  different analysis algorithm may be more appropriate). */
  void setSpectralOversampling(int newFactor) 
  { 
    rsAssert(rsIsPowerOfTwo(newFactor));
    zeroPad = newFactor; 
  }
  // rename to setZeroPadFactor

  /** Sets the number of cycles in each analyzed block. Must be a power of two. */
  void setNumCyclesPerBlock(int newNumCycles) { cyclesPerBlock = newNumCycles; }

  /** Sets the type of window to be used. */
  void setWindowType(rsWindowFunction::WindowType newType)
  { 
    rsAssert(isWindowTypeSupported(newType), "Desired window type not supported");
    windowType = newType; 
  }

  /** Window sidelobe-rejection parameter in dB, applicable only when windowType is dolphChebychev.
  40 dB is hamming-like, 60 dB is blackman-like. */
  void setSidelobeRejection(T rejection) { sidelobeRejection = rejection; }
  // the number of cycles per block should be increased when the rejection is increased (i think, 
  // this is because of the wider mainlobe-width) - todo: figure out, how exactly - also allow all
  // integers, not just powers of two, for allowing the desired functional dependency to be 
  // realized more accurately


  /** Sets up, whether or not inharmonic partials should be expected. If this is set to true (and 
  the number cycles per block is > 1), the algorithm tries to find the actual partial frequency by
  searching in some neighbourhood of the expected harmonic "slot". If false, it will just read out
  the spectrum at the expected exact harmonic frequency. */
  void setAllowInharmonics(bool allow) { allowInharmonics = allow; }


  /** Sets the relative width inside which we search for a spectral peak in the vicinity of an
  expceted harmonic. The absolute width in bins should be proportional to the window's mainlobe 
  width and the zero-padding factor. This function sets the proportionality constant 
  (default: 1). */
  void setSpectralPeakSearchWidth(T newWidth) { peakSearchWidth = newWidth; }

  /** When multiple cycles are used and inharmonic partials are allowed, we must search for the 
  mainlobe peak in the neighbourhood of the expected harmonic frequency. This search also involves 
  a decision whether or not a found peak/lobe is to be considered a mainlobe (we want avoid picking
  up sidelobes as partials). The decision is made based on the relative width of the measured lobe
  with respect to the window's mainlobe-width - this function sets a proportionality factor. The 
  found lobe must be at least as wide as this factor times the mainlobe width of the window in 
  order to be considered a partial - if its narrower, we assume a sidelobe peak and discard the 
  peak. Default value is 0.75 (maybe 0.5 would be better? ...tests needed...) */
  void setMinPeakToMainlobeWidthRatio(T newWidth) { minPeakToMainlobeWidthRatio = newWidth; }
  // i think, the optimal value should depend on the ratio of mainlobe-width to the width of the 
  // widest sidelobe - if the widest sidelobe is 0.5 times as wide as the mainlobe, we should use a 
  // factor > 0.5 or something - maybe measure this ratio for the Chebychev window - it will 
  // probably depend on the length and certainly on the attenuation

  void setMinPeakToHarmonicWidthRatio(T newWidth) { minPeakToHarmonicWidthRatio = newWidth; }



  // void setTemporalOversampling(int newFactor)
  // ...produce intermediate datapoints between the already existing ones...

  // not yet used - todo: make these limits functional:
  //void setMinPartialIndex(int newIndex) { minPartialIndex = newIndex; }
  // does not yet work due to treating DC outside the loop in fillHarmonicData

  /** Sets the maximum harmonic index to be analyzed */
  void setMaxPartialIndex(int newIndex) { maxPartialIndex = newIndex; }

  /** Sets the maximum number of partials to be an (including DC) */
  //void setMaxNumPartials(int newMaximum) { maxNumPartials = newMaximum; }

  /** This option can be used to remove any harmonics that exceed the Nyquist limit, even if just 
  temporarily. The analysis may prodcue such frequencies due to the fact that the original audio is 
  stretched before analysis and the post-processing then shifts all frequencies up. However, if a 
  (short) cycle gets stretched out a lot, it will not actually contain any high-freq content above 
  the original Nyquist rate - so when a partial in the model goes above that limit, it will usually
  have zero amplitude. ...it turned out that this option may sometimes remove valid high-freq 
  content, so it's off by default ...more research needed....  */
  void removePotentiallyAliasingHarmonics(bool shouldRemove)
  {
    antiAlias = shouldRemove;
  }


  void setFreqsByPhaseDerivative(bool shouldRefine) 
  { 
    freqsByPhaseDerivative = shouldRefine; 
  }

  void setFreqPhaseConsistency(bool shouldBeConsistent)
  {
    freqsPhaseConsistent = shouldBeConsistent;
  }

  //---------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, iff the given desired window type is one of our supported window types. */
  bool isWindowTypeSupported(rsWindowFunction::WindowType desiredType)
  { 
    typedef rsWindowFunction::WindowType WT;
    WT dt = desiredType;
    return dt == WT::rectangular || dt == WT::hamming || dt == WT::blackman 
        || dt == WT::blackmanHarris || dt == WT::dolphChebychev;
  }
  // todo: support more types. For the non-parameterized types, this is trivial - just add it to 
  // the or-chain. But for other parameterized types (such as Kaiser), we need to expose the 
  // relevant parameter via a setter and have to do a switch in fillWindow because for different
  // window types, the parameter may have different meanings....

  /** If a parameterized window function is used, this returns the numeric value of the window's
  parameter. For a Dolph-Chebychev window, this parameter means the attenuation of the sidelobes
  in dB. */
  T getWindowParameter() const
  {
    return sidelobeRejection;
    // If we want to use other parameterized windows later, we may need a switch here - currently, 
    // the only supported parameterized window is the Dolph-Chebychev window
  }


  std::vector<T> getOriginalTimeStamps() { return tIn; }

  std::vector<T> getWarpedTimeStamps() { return tOut; }

  /** Returns a reference to the embedded rsCycleMarkFinder object that is used to find the 
  cycle-marks for the time-warping, to give client code access to its settings. */
  rsCycleMarkFinder<T>& getCycleFinder() { return cycleFinder; }

  //---------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Analyzes the given input sound and returns the sinusoidal model object that models the given 
  sound. */
  RAPT::rsSinusoidalModel<T> analyze(T* sampleData, int numSamples);


  /** \name Misc */

  /** Copies the content of the given signal block into the given transform buffer, appropriately
  padding with zeros, if necessary. Made static to enable it to be conveniently tested. */
  static void prepareBuffer(const std::vector<T>& signalBlock, std::vector<T>& trafoBuffer);


  bool useOldCode = false; 
  // temporary, for development - switch between new and old implementation of filling the blocks
  // i.e. between single-cycle and multi-cycle implementation

protected:


  /** The first step in the analysis algo is to pre-process the audio by flattening the pitch, 
  which is done by this function. The flattened signal will be stored in our member array y and 
  this function will also fill the member arrays tIn, tOut that contain the time warping map that
  has been used for flattening (sampled at instants of cycle marks). These arrays will be later 
  used for post processing the model data to account for the pitch-flattening. The boolean return
  value informs, if the process was successful (it will fail, if it can't find at least 2 cycle 
  marks). */
  bool flattenPitch(T* sampleData, int numSamples);

  /** Sets the maximum cycle length that was measured in the original, unflattened signal and does
  all the required adjustments of blockSize, trafoSize, etc. and re-allocates our buffers, and 
  re-creates the window and so on, if necessary. Called from flattenPitch during the 
  pre-processing step. */
  void setMaxMeasuredCycleLength(T maxLength);

  /** The second step in the analysis algo is to perform an FFT on each cycle of pitch-flattened
  signal. Because the pitch is now flat, each cycle has the same length (which was chosen to be
  a power of two, greater or equal to the length of the longest cycle in the input signal). This
  fills in the model data with (preliminary) values. */
  void analyzeHarmonicsOld(RAPT::rsSinusoidalModel<T>& mdl);

  void analyzeHarmonics(RAPT::rsSinusoidalModel<T>& mdl);
  // new version that supports multi-cycle blocks - under construction - when finished, the old
  // version may be deleted

  /** The third step in the analysis algo is to modify the time and frequency data to account for 
  the pitch flattening that was done in the first step. */
  void deFlattenPitch(RAPT::rsSinusoidalModel<T>& mdl);

  /** Removes those partials in the model that are above the original Nyquist freqeuncy. Such 
  partials may be generated by the analysis due to the fact that original signal is pitch-shifted
  down in the flattening step. They will typically have amplitude close to zero because, of 
  course, the original signal cannot contain these frequencies - except at the edges - there they
  may indeed have nonzero amplitudes if the original signal is harshly cut off there. Anyway - 
  zero amplitude or not - we don't want them in the model - their mere existence is an artifact
  of the downshifting step. */
  void removeAliasing(RAPT::rsSinusoidalModel<T>& mdl);
  // todo: prevent them from being produces in the first place - this will speed up the analysis
  // because we will have to analyze fewer partials

  /** Fills the initial datapoint at time zero and the final datapoint at time 
  (numSamples-1)/sampleRate - these are treated separately. */
  void handleEdges(RAPT::rsSinusoidalModel<T>& mdl);

  /** During most of our computational steps in the algo, we represent time in units of samples,
  but ultimately, rsSinusoidalModel wants to have time values in seconds, so this function is used
  as final step to convert all values. */
  void convertTimeUnit(RAPT::rsSinusoidalModel<T>& mdl);

  /** Refines the frequency estimates in the model, if the respective options are set to true (this
  step is optional). */
  void refineFrequencies(RAPT::rsSinusoidalModel<T>& mdl);

  /** Returns length of time-warping map (sampled at cycle marks). */
  int getMapLength() const { return (int) tIn.size(); }
  // rename to getTimeWarpMapLength
   
  /** Returns the number of FFT analysis frames. */
  int getNumFrames() const { return getMapLength()-1; }

  /** Returns the number of analyzed harmonics (including DC). */
  int getNumHarmonics() const 
  { 
    if(maxPartialIndex == -1)
      return cycleLength / 2;
    else
      return rsMin(cycleLength / 2, maxPartialIndex) + 1; // +1 correct?
  } 
  // rename to getNumPartials (maybe)
  // todo: maybe decide in advance, which harmonics will or will not alias and analyze only those
  // which won't as an optimization (instead of analyzing them all and then discarding the aliasing 
  // ones) - in this case, return the actual number of harmonics here

  /** Returns the FFT bin index for the partial with given harmonic index. This may differ from the
  harmonic index due to having multiple cycles within an analysis frame, using zero-padding before 
  the FFT and omitting low-frequency partials */
  //int getPartialBinIndex(int partialIndex) const 
  //{ return cyclesPerBlock*zeroPad*(partialIndex-minPartialIndex); } 


  /** Returns the number of datapoints (per partial) in the sinusoidal model. */
  int getNumDataPoints() const { return getNumFrames() + 2; } // + 2 for fade in/out frames

  /** Computes and returns an array of cycle-marks. */
  std::vector<T> findCycleMarks(T* x, int N);

  /** Returns the time instant that should be associated with the frame of given index 
  (before post-processing). The time stamp will be based on the content of tOut member array. */
  T getTimeStampForFrame(int frameIndex);
  // maybe rename to getWarpedTimeForFrame and have a corresponding getUnWarpedTimeForFrame
  // function (that uses tIn instead of tOut)...this will actually be used to *replace* the 
  // original time-stamp data later in the post-processing - so we may not actually need to 
  // compute the values before post-processing - but for conceptual clarity of the algorithm, 
  // we just do it for the time being (the effort is neglibilibe anyway)

  T getUnWarpedTimeStampForFrame(int frameIndex);
  // maybe rename to getUnWarpedSampleIndexForFrame

  /** Returns the extent by which we look to the left and and right for searching for a peak in the 
  magnitude spectrum from a given expected peak bin index k. If the returned value is w, we look at
  bins k-w...k+w. This width is proportional to the zero-padding factor and the window's mainlobe 
  width.
  ..should we include or exclude the boundaries?
  */
  int getSpectralPeakSearchWidth();

  /** Used internally to fill in the data in the model at the given frame-index based on the 
  current content of our "sig" buffer member variable. The time-stamp of the frame should be passed
  by the caller - frequency, magnitude and phase data are computed from the FFT of the sig 
  buffer. */
  void fillHarmonicData(RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T timeStamp);

  /** Given a vector of FFT magnitude values, this function returns a vector of bin indices that
  correspond to partials. */
  //std::vector<int> findPartialBins(const std::vector<T> magnitudes);

  /** Given a vector of spectral magnitude values, this function returns the index of a local 
  maximum in the neighbourhood of a given centerBin. The centerBin is typically the bin, where we
  expect a harmonic/partial to be, if it were perfectly harmonic. We search through all bins from 
  centerBin - halfSearchWidth to centerBin + halfSearchWidth. If no local maximum is found (or some
  other constraints that indicate a proper partial are not met), -1 is returned which indicates 
  that there is no partial in the neighbourhood of the expected harmonic. */
  int findPeakBinNear(std::vector<T>& v, int centerBin, int halfSearchWidth);
  // make vector const again (but that clashes with plotting -> make GNUPlotter take const arrays)



  bool isPeakPartial(std::vector<T>& v, int peakBin);


  /** Fills our window function array. */
  void fillWindow();




  //RAPT::rsPitchFlattener<T, T> flattener;
  rsFourierTransformerRadix2<T> trafo;
  rsCycleMarkFinder<T> cycleFinder;


  // block/transform buffer sizes:

  // old: 
  int cycleLength    = 0;  // length of one cycle in samples (after pitch-flattening)
  int cyclesPerBlock = 1;  // number of cycles per block/window, power of 2

  // new:
  //T   cycleLength    = T(0);  // length of one cycle in samples (after pitch-flattening)
  //T   cyclesPerBlock = T(1);  // number of cycles per block/window

  int blockSize      = 0;  // analysis block size == cycleLength * cyclesPerBlock
  int zeroPad        = 1;  // zero padding factor for FFT, power of 2
  int trafoSize      = 0;  // FFT size == blockSize * zeroPad

  
  typedef rsWindowFunction::WindowType WindowType;
  WindowType windowType = WindowType::rectangular;
  // todo: switch to the Dolph-Chebychev window and let the user adjust its parameter - either in 
  // terms of the sidelobe-rejection (in dB) or in mainlobe-width (in bins? or Hz? or some other 
  // normalized unit?). there are some rules, how the number-of-cycles per block should be related
  // to the mainlobe-width - wider mainlobes (i.e. higher sidelobe-rejection) require more cycles 
  // per block - document these things...
  T sidelobeRejection = 40; // sidelobe-rejection parameter in dB, applicable only when windowType 
                            // is dolphChebychev, 40 dB is hamming-like, 60 blackman-like


  bool antiAlias   = false;

  //bool freqReAssignment = false;   // doesn't seem to make sense, uses up phase information
  bool freqsByPhaseDerivative = false;
  bool freqsPhaseConsistent = false;


  bool allowInharmonics = true;
  bool parabolicInterpolation = true;   // relevant only allowInharmonics = true
  bool phaseInterpolation = true;       // relevant only when parabolicInterpolation = true


  T minPeakToMainlobeWidthRatio = T(0.75);
  T minPeakToHarmonicWidthRatio = T(0.75);
  T peakSearchWidth = T(1); 
    // relative width in which we search for spectral peaks around an expected harmonic, adjusted 
    // according to windowType and zeroPad

  T sampleRate = 1;
  T sincLength = 512.0;  // length of sinc-interpolator for time-warping


  // limits to the range of analyzed partials:
  //int minPartialIndex = 0; 
  //int maxPartialIndex = std::numeric_limits<int>::max();
  int maxPartialIndex = -1;  // -1 encodes "no upper limit" - maybe use std::optional
  //int maxNumPartials = std::numeric_limits<int>::max();



  std::vector<T> y;     // pre-processed (time-warped) signal



  // buffers:
  std::vector<T> tIn, tOut;      // time-warping map (sampled at cycle marks)
  // maybe rename to "to", "tw" for original and warped

  std::vector<T> sig, wnd, sigPadded, mag, phs;  // block of signal, window, magnitude and phase

};