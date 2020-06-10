

//=================================================================================================

template<class T>
int rsSinusoidalAnalyzer<T>::getRequiredBlockSize(rsWindowFunction::WindowType type, T df, T fs, 
  bool oddSize)
{
  T B = rsWindowFunction::getMainLobeWidth(type, 0.0, 0); // later maybe pass a window parameter
  int size = (int) ceil(B*fs/df);
  if(oddSize && rsIsEven(size))   // maybe remove, if linear interpolation is used for the phase
    size += 1;
  return size;
  // Because the phase estimation is much more accurate for odd block sizes (because the phase 
  // spectrum of window functions is flat over the mainlobe for odd block sizes, for even sizes it 
  // has a linear trend), this will always return an odd number so you may use the returned value 
  // as your blockSize without further ado. or: maybe, if linear interpolation is used for phase 
  // estimation, we can get equally good results with even block-sizes? -> try it!
}

template<class T>
T rsSinusoidalAnalyzer<T>::getRequiredThreshold(rsWindowFunction::WindowType type, T margin)
{
  return rsWindowFunction::getSideLobeLevel(type, 0.0) + margin;
  // todo: maybe if we constrain the peak-conditions further, such that a peak must be higher than
  // 2,3,4,... of its neighbours to the left and right, we can use lower thresholds without picking 
  // up sidelobes and maybe this number should depend on the mainlobe width of the window 
  // ...hmm...but those peaks will be buried in sidelobes from other partials anyway
}

template<class T>
size_t rsSinusoidalAnalyzer<T>::findBestMatchingTrack(T freq, 
  std::vector<RAPT::rsSinusoidalPartial<double>>& tracks, 
  const std::vector<bool>& trackContinued) const
{
  T dfMin = RS_INF(T);  
  size_t bestIndex = 0;
  for(size_t i = 0; i < tracks.size(); i++) {
    T trackFreq = tracks[i].getEndFreq();
    T df = rsAbs(freq - trackFreq);
    if(df <= dfMin && trackContinued[i] == false) { // look only at tracks that are not already continued
      dfMin = df;
      bestIndex = i;
    }
  }
  if(dfMin < getMaxFreqDelta(freq))
    return bestIndex;
  else
    return tracks.size();
}

template<class T>
void rsSinusoidalAnalyzer<T>::continuePartialTracks1(
  std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeaks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& aliveTracks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& deadTracks) const
{
  // initializations:
  typedef std::pair<size_t, size_t> IndexPair;
  std::vector<IndexPair> continuations;               // 1st: activeTracks, 2nd: newPeaks
  std::vector<size_t> deaths;                         // index into activeTracks
  std::vector<size_t> births;                         // index into newPeaks
  size_t pkIdx, trkIdx;                               // peak index, track index
  size_t numActiveTracks = aliveTracks.size();
  std::vector<bool> trackContinued(numActiveTracks);  // true, if track has already been continued
  for(trkIdx = 0; trkIdx < numActiveTracks; trkIdx++)
    trackContinued[trkIdx] = false;

  // Figure out which tracks have to be continued (and if so, which new peak data should be 
  // appended), which tracks have to be discontinued ("death") and for which peaks a new track
  // has to be created ("birth"):

  // loop over the new peaks to figure out birthes and continuations:
  for(pkIdx = 0; pkIdx < newPeaks.size(); pkIdx++) {
    trkIdx = findBestMatchingTrack(newPeaks[pkIdx].freq, aliveTracks, trackContinued); // looks only in those that are not already continued
    if(trkIdx == aliveTracks.size())     // no match found, so a new track is born
      births.push_back(pkIdx);
    else {
      continuations.push_back(IndexPair(trkIdx, pkIdx));
      trackContinued[trkIdx] = true;
    }
  }
  // Actually, it would be better to run the loop over the active tracks and find a match for
  // each track among the new peaks. This would allow to search for a frequency that is not exactly 
  // the last value in the respective track but one that is obtained by extrapolation of the 
  // frequency trajectory so far (linear prediction or polynomial extrapolation). Also, the search 
  // could take advantage of the fact that the peak array is sorted by frequency (to use binary 
  // instead of linear search). This would imply that we would have to work with a boolean 
  // "peakUsed" array instead of "trackContinued". Maybe keep both variants in a prototype 
  // implementation. ...well, actually, it would be possible to use extrapolation with the current
  // loop over the peaks as well - we'd just have to replace
  // T trackFreq = tracks[i].getEndFreq();
  // by
  // T trackFreq = getTrackFreq(tracks[i])
  // where getTrackFreq would be a function that performs the extrapolation - but still the 
  // advantage of binary serach would remain - moreover, looping over the tracks is the way it's
  // described in the literature as well. ..maybe try to cook up a situation where this actually 
  // makes a difference - i think, it makes a difference, if there are two competing possible 
  // matchings - in this function, the one with the lower track-index is chosen, in the other, the 
  // one with the loew peak index is chosen - actually, in this case, the one with the smaller 
  // frequency difference should be chosen - how can we achieve this? maybe first collect all 
  // candidate matchings and if a peak or track index appears twice (or more often), select one?
  // that would let us get rid of the bool-arrays and the difference between the two versions of 
  // the algo disappears, so we may keep only one -> check, how Serra's SMS tools do it

  // loop over the "trackContinued" array to figure out deaths:
  for(trkIdx = 0; trkIdx < numActiveTracks; trkIdx++) {
    if(trackContinued[trkIdx] == false)
      deaths.push_back(trkIdx);
  }

  // we have figured out the deaths, births and continuations - apply them:
  applyContinuations(newPeaks, aliveTracks, deadTracks, births, deaths, continuations);
}

template<class T>
size_t rsSinusoidalAnalyzer<T>::findBestMatchingPeak(T freq, 
  std::vector<RAPT::rsInstantaneousSineParams<T>>& peaks,
  const std::vector<bool>& peakUsed) const
{
  T dfMin = RS_INF(T);  
  size_t bestIndex = 0;
  for(size_t i = 0; i < peaks.size(); i++) {
    T peakFreq = peaks[i].freq;
    T df = rsAbs(freq - peakFreq);
    if(df <= dfMin && peakUsed[i] == false) { // look only at peaks that are not already used up
      dfMin = df;
      bestIndex = i;
    }
  }
  if(dfMin < getMaxFreqDelta(freq))
    return bestIndex;
  else
    return peaks.size();
}

template<class T>
void rsSinusoidalAnalyzer<T>::continuePartialTracks0(
  std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeaks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& aliveTracks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& deadTracks) const
{
  // initializations:
  typedef std::pair<size_t, size_t> IndexPair;
  std::vector<IndexPair> continuations;               // 1st: aliveTracks, 2nd: newPeaks
  std::vector<size_t> deaths;                         // index into aliveTracks
  std::vector<size_t> births;                         // index into newPeaks
  size_t pkIdx, trkIdx;                               // peak index, track index
  size_t numNewPeaks = newPeaks.size();
  std::vector<bool> peakUsed(numNewPeaks);            // true, if peak is already used up
  for(pkIdx = 0; pkIdx < numNewPeaks; pkIdx++)
    peakUsed[pkIdx] = false;

  // loop over the alive tracks to figure out deaths and continuations:
  for(trkIdx = 0; trkIdx < aliveTracks.size(); trkIdx++) {
    pkIdx = findBestMatchingPeak(aliveTracks[trkIdx].getEndFreq(), newPeaks, peakUsed); // looks only at peaks that are not already used up
    if(pkIdx == newPeaks.size())        // no match found, so the track dies
      deaths.push_back(trkIdx);
    else {
      continuations.push_back(IndexPair(trkIdx, pkIdx));
      peakUsed[pkIdx] = true;
    }
  }

  // loop over the "peakUsed" array to figure out births:
  for(pkIdx = 0; pkIdx < numNewPeaks; pkIdx++) {
    if(peakUsed[pkIdx] == false)
      births.push_back(pkIdx);
  }

  // we have figured out the deaths, births and continuations - apply them:
  applyContinuations(newPeaks, aliveTracks, deadTracks, births, deaths, continuations);
}

template<class T>
void rsSinusoidalAnalyzer<T>::applyContinuations(
  std::vector<RAPT::rsInstantaneousSineParams<T>>& newPeaks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& aliveTracks,
  std::vector<RAPT::rsSinusoidalPartial<T>>& deadTracks,
  std::vector<size_t>& births, std::vector<size_t>& deaths,
  std::vector<std::pair<size_t, size_t>>& continuations) const
{
  size_t pkIdx;   // peak index
  size_t trkIdx;  // track index
  int i;          // signed because we use a >= 0 comparison in the deaths loop - 
                  // why actually? - because in the loop, we remove them from the end of the 
                  // aliveTracks array to not mess up the indices inside the loop

  // continue matched tracks with new peaks:
  for(i = 0; i < continuations.size(); i++) {
    trkIdx = continuations[i].first;
    pkIdx  = continuations[i].second;
    aliveTracks[trkIdx].appendDataPoint(newPeaks[pkIdx]);
  }

  // maybe we should assert that deaths is sorted ascending - because when we remove tracks from
  // the aliveTracks array, we must do so from behind in order to not mess up the array indexes
  // inside the loop...

  rsSinusoidalPartial<T> track;
  rsInstantaneousSineParams<T> params;

  // kill discontinued tracks by appending an additonal "fade-out" datapoint at the end of the 
  // track and moving them from the aliveTracks to the deadTracks array:
  if(deaths.size() > 0) {  // can we get rid of this if? i think so
    for(i = (int) deaths.size()-1; i >= 0; i--) {
      trkIdx = deaths[i];
      track  = aliveTracks[trkIdx];
      if(fadeOutTime > T(0))
        track.applyFadeOut(fadeOutTime);
      rsAppend(deadTracks, track);
      rsRemove(aliveTracks, trkIdx); 
    }
  }

  // create new tracks by creating a fresh track from the peaks that should give birth and also 
  // prepend a "fade-in" datapoint:
  for(i = 0; i < births.size(); i++) {
    pkIdx  = births[i];
    params = newPeaks[pkIdx];
    RAPT::rsSinusoidalPartial<T> newTrack;
    newTrack.appendDataPoint(params);
    if(fadeInTime > T(0)) 
      newTrack.applyFadeIn(fadeInTime);
    aliveTracks.push_back(newTrack);
  }
}

template<class T>
rsMatrix<std::complex<T>> rsSinusoidalAnalyzer<T>::getComplexSpectrogram(
  T* sampleData, int numSamples)
{
  return sp.getComplexSpectrogram(sampleData, numSamples);
}

template<class T>
RAPT::rsSinusoidalModel<T> rsSinusoidalAnalyzer<T>::analyzeSpectrogram(
  const RAPT::rsMatrix<std::complex<T>>& stft, T sampleRate)
{
  // Initializations:
  typedef RAPT::rsInstantaneousSineParams<double> InstParams;
  typedef RAPT::rsSinusoidalPartial<double> Partial;  // rename to SineTrack
  int numFrames  = stft.getNumRows();
  int numBins    = stft.getNumColumns();
  int firstFrame = 0; 
  int lastFrame  = numFrames-1;
  int frameIndex  = firstFrame;
  T binDelta   = sampleRate / sp.getFftSize();
  T frameDelta = sp.getHopSize() / sampleRate;


  rsMatrix<T> mag = matrixMagnitudes(stft);
  rsMatrix<T> phs = matrixPhases(stft);

  // ...maybe plot the spectrogram here...
  //plotPhasogram(numFrames, numBins, phs.getDataPointer(), sampleRate, sp.getHopSize());

  std::vector<Partial> activeTracks, finishedTracks;
  std::vector<InstParams> instPeakParams;

  //bool plotShortTimeSpectra = true;

  // needed for plotting only:
  std::vector<T> freqs(numBins);
  for(int i = 0; i < numBins; i++)
  {
    //freqs[i] = i * sampleRate / sp.getFftSize(); // plot against actual frequency
    freqs[i] = i;  // plot against bin idex
  }

  // loop over the frames:
  while(frameIndex != lastFrame) {

    T time  = frameIndex * frameDelta;
    T* pMag = mag.getRowPointer(frameIndex);
    T* pPhs = phs.getRowPointer(frameIndex);
    const std::complex<T>* pCmp = stft.getRowPointerConst(frameIndex);  // pointer to complex short-time spectrum


    //plotData(numBins, pMag);
    //plotDecibels(numBins, &freqs[0], pMag); // for development
    //plotData(numBins, &freqs[0], pMag); // for development
    //plotData(numBins, &freqs[0], pPhs);

    // find spectral peaks:
    std::vector<int> peaks = peakIndices(pMag, numBins, magThreshold);

    // determine exact peak frequencies, amplitudes and phases:
    instPeakParams.resize(peaks.size());
    for(size_t i = 0; i < peaks.size(); i++) {
      T peakBin, peakAmp;
      spectralMaximumPositionAndValue(pMag, peaks[i], &peakBin, &peakAmp);

      //T peakFreq  = peakBin*sampleRate/sp.getFftSize();

      // maybe later let the user select the phase-interpolation method between no-interpolation, 
      // linear and maybe others (maybe interpolate in the complex domain and then extract phase):
      T peakPhase = pPhs[peaks[i]];      // maybe use interpolation later

      //T peakPhase = RAPT::rsArrayTools::interpolatedValueAt(pPhs, numBins, peakBin);
      // simple linear interpolation is wrong here - it returns a totally wrong value when the 
      // phases of the bins are close to -pi and pi or vice versa - we need wrapped interpolation
      // ...but really unit-test the wrapped interpolation first

      // wrapped linear interpolation:
      //int k = peaks[i];
      //peakPhase = rsInterpolateWrapped(pPhs[k], pPhs[k+1], peakBin-floor(peakBin), -PI, PI);
      // strangely, with this interpolation, we may get louder residuals (w Hamming, zp=2) - maybe 
      // make it optional and do more tests
      // i think, it's because we need to dsitinguish two cases: 
      // peakBin >= peaks[i], peakBin < peaks[i]
      // todo: use this:
      // peakPhase = interpolatePhase(pPhs, peakBin);

      // maybe also try to interpolate the real and imaginary parts and the extract the phase and
      // do tests, which one is best


      instPeakParams[i].time  = time;
      instPeakParams[i].freq  = peakBin*sampleRate/sp.getFftSize(); // use getBinFrequency(peakBin);
      instPeakParams[i].gain  = peakAmp;
      instPeakParams[i].phase = peakPhase;
    }

    // peak continuation, birth or death - todo: dispatch between tracking algorithms:
    continuePartialTracks0(instPeakParams, activeTracks, finishedTracks);
    frameIndex += 1;
  }

  // create and return the model data structure:
  rsSinusoidalModel<T> model;
  model.addPartials(finishedTracks);       // add the finished tracks to the model
  if(fadeOutTime > T(0))                   // "finalize" active tracks by applying fade-out...
    for(size_t i = 0; i < activeTracks.size(); i++)
      activeTracks[i].applyFadeOut(fadeOutTime);
  model.addPartials(activeTracks);         // ...and add them, too
  cleanUpModel(model);                     // remove spurious tracks, merge, etc.
  return model;


  // algorithm:

  // -scan through this spectrogram from left to right
  // -for each frame m, do:
  //  -find the spectral peaks in the frame
  //  -for each peak k in the frame, do:
  //   -figure out its exact frequency, amplitude and phase (and maybe later time, if re-assignment 
  //    is used)
  //   -try to find a match in the set of active partials:
  //    -if a match is found, append the new datapoint to that matched partial
  //    -if no match is found, start a new partial 
  //   -if there are active partials that have not been used up by this matching procedure,
  //    they end now

  // -todo: post-processing: after this is all done, go over the data a second time to see if we 
  //  can merge partials and/or remove partials whose length is below a user-defined threshold

  // -find peak frequency by (parabolic) interpolation (or maybe cubic?)
  // -evaluate complex amplitude at peak-freq
  // -compute magnitude phase from interpolated complex amplitude and store in one of the partials

  // references:
  // http://www.cerlsoundgroup.org/Loris/
  // https://ccrma.stanford.edu/~juan/ATS_manual.html
  // http://clam-project.org/
}

template<class T>
std::vector<int> rsSinusoidalAnalyzer<T>::peakIndices(T* x, int N, T threshToMax)
{
  T max = RAPT::rsArrayTools::maxValue(x, N);
  std::vector<int> peaks;
  for(int i = 1; i < N-1; i++)
  {
    T dbg = x[i];
    if(x[i] > x[i-1] && x[i] > x[i+1] && x[i] > threshToMax*max)
      peaks.push_back(i);
  }
  return peaks;
} 

template<class T> 
void rsSinusoidalAnalyzer<T>::spectralMaximumPositionAndValue(T *x, int k, T* pos, T* val)
{
  // find coeffs of parabolic interpolant (maybe factor out, so we can plot the parabola):
  T lowAmp = 0.0000000001; // -200 dB - to prevent log-of-zero
  T a[3], y[3];
  y[0] = rsAmpToDbWithCheck(x[k-1], lowAmp);   // left
  y[1] = rsAmpToDbWithCheck(x[k],   lowAmp);   // mid
  y[2] = rsAmpToDbWithCheck(x[k+1], lowAmp);   // right
  rsPolynomial<T>::fitQuadratic_m1_0_1(a, y);

  // find maximum position and evaluate parabola there (unless parabola is degenerate):
  if(a[2] != T(0)) { 
    T d  = rsPolynomial<T>::quadraticExtremumPosition(a);
    *pos = k + d;
    *val = rsDbToAmp(a[0] + (a[1] + a[2]*d)*d);
  } 
  else {
    *pos = T(k);
    *val = x[k];
  }
}

template<class T> 
T rsSinusoidalAnalyzer<T>::interpolatePhase(T* phs, T pos)
{
  int k = (int) floor(pos);
  T frac = pos - T(k);
  return rsInterpolateWrapped(phs[k], phs[k+1], frac, -PI, PI);
}

template<class T>
rsSinusoidalModel<T> rsSinusoidalAnalyzer<T>::analyze(
  T* sampleData, int numSamples, T sampleRate)
{
  // -maybe pre-process the input signal by flattening the pitch and make the period coincide with
  //  a fraction of the analysis frame size 
  // -we should keep the time warping map around and when done with the analysis, use it to 
  //  re-map the time-instants in the model
  // -we should also keep the corresponding (instantaneous) frequency scaler map around and apply
  //  the inverse to the freq data in the model
  //  ...but later - first, let's see how far we get without such a pre-processing and try to make 
  //  the algorithm work well under this (suboptimal) condition, too

  rsMatrix<std::complex<T>> stft = getComplexSpectrogram(sampleData, numSamples);
  return analyzeSpectrogram(stft, sampleRate);
}

template<class T>
void rsSinusoidalAnalyzer<T>::cleanUpModel(rsSinusoidalModel<T>& model) const
{
  // todo: merge partials (before removing spurious ones)

  // remove spurious partials (maybe factor out):
  for(int i = (int)model.getNumPartials()-1; i >= 0; i--)  
    if(model.getPartial(i).getLength() <= minTrackLength) 
      model.removePartial(i);
  // We use <= instead of < to remove zero-length tracks, even when minTrackLength == 0, so we 
  // allow only for tracks that are *strictly* longer than 0 when minTrackLength == 0. A track with
  // exactly zero length makes no sense, so it should not occur - each partial must have a start and 
  // end and a finite time interval in between

  // de-bias the frequency estimates in the remaining, non-spurious partials:
  if(forceFreqPhaseConsistency)
    for(int i = 0; i < model.getNumPartials(); i++)
    {
      rsSinusoidalProcessor<T>::makeFreqsConsistentWithPhases(model.getModifiablePartialRef(i));

      //makeFreqsConsistentWithPhases(model.getModifiablePartialRef(i));
      //model.getModifiablePartialRef(i).makeFreqsConsistentWithPhases();
      // todo: have a function model.makeFreqsConsistentWithPhases(); that wraps this loop over
      // the partals
    }

  // maybe have an option to delete redundant datapoints, i.e. datapoints that remain constant
  // over a long time, such as 1000Hz, 1000Hz, 1000Hz, ...., 1000Hz, 1000Hz - just keep the outer
  // two or maybe the outer 4 (because of the interpolation) economizeModel or something
  // or simplifyModel, simplifyPartial
}

//template class rsSinusoidalAnalyzer<double>;


/*

Ideas:

Peak-Tracking:
-maybe instead of looping over the tracks or peaks and greedily picking the closest match in the 
 respective other array, we should first collect all candidate matches and then decide which gets
 assigned to which based on the whole picture
-maybe a single peak can actually correspond to two tracks (if two tracks cross each other in the
 current frame) and maybe the peak's amplitude should then be split between the two tracks

Window:
-use a Dolph-Chebychev window and let the user only select the amplitude threshold for partials
 -the window will then be automatically tuned such the sidelobes are some reasonable margin below 
  the detection threshold (the margin has to figured out empirically - maybe 15 or 20 dB)
 -this, in turn, determines the mainlobe width of the window which in turn determines the required
  block-size
 ...soo in the high-level setup, the user just selects frequency resolution amplitude threshold 
 (and possibly the margin - but maybe not) and anything else can then be set up automatically
 ...nevertheless, provide the low-level setup anyway
 -maybe the allowed frequency deviation between two frames (which should be proportional to the 
  hopSize) can be set up in a way that is independent from the hopSize - maybe in Hz/s
  10 Hz/s means the frequency is allowed to change by 10 Hz in one second
 -but also have a look at the phase-response - make sure, the window doesn't mess up phase 
  measurement

-when a Chebychev window is used, we may make the block-size a function of the threshold (because
 the width of the mainlobe gets smaller with higher thresholds, we may choose smaller block-sizes
-use an adaptive (time-varying) amplitude threshold - high at the beginning (like -20..-30) - we
 can use smaller block-sizes in transient portions of the sound and therefore have better 
 time-resolution there. the idea is that at the beginning, all partials still have high amplitude,
 so a higher threshold is justified - later in the decay phase, we can use lower thresholds 
 (-40..-80) to not loose track of the more highly damped partials (...this applies mostly to
 decaying sounds, like guitar - for steady sounds, like flute, the situation is different)
-maybe we should first figure out, where the transients are via an onset-detector


// how do we best compute the instantaneous phase - linear interpolation?
int binInt  = peaks[i];
T   binFrac = peakBin-binInt;
if(binFrac < 0) {
binInt  -= 1;
binFrac  = 1-binFrac;
}
//T amp2 = (1-binFrac)*pMag[binInt] + binFrac*pMag[binInt+1]; // another way to compute amplitude - vary bad
//std::complex<T> binCmpVal = (1-binFrac)*pCmp[binInt] + binFrac*pCmp[binInt+1]; 
//T binMag = abs(binCmpVal);
//T binPhs = arg(binCmpVal);
// complex value at peak-bin
T phs0 = pPhs[binInt];   // just for debug
T phs1 = pPhs[binInt+1]; // dito
peakPhase = rsInterpolateWrapped(pPhs[binInt], pPhs[binInt+1], binFrac, -PI, PI);
// hmm...it doesn't seem to make sense to interpolate the phase between bins like that - it 
// does not seem to be a smooth function - it jumps erratically between bins - maybe we need 
// a completely different way to obtain the instantaneus phase...maybe by comparing values in
// different frames? or maybe by correlating the function with a complex sine that has 
// exactly the frequency? maybe the Goertzel algo can be used?:
// https://en.wikipedia.org/wiki/Goertzel_algorithm
// https://www.embedded.com/design/configurable-systems/4006427/A-DSP-algorithm-for-frequency-analysis
// https://www.st.com/content/ccc/resource/technical/document/design_tip/group0/20/06/95/0b/c3/8d/4a/7b/DM00446805/files/DM00446805.pdf/jcr:content/translations/en.DM00446805.pdf
// https://stackoverflow.com/questions/13499852/scipy-fourier-transform-of-a-few-selected-frequencies
// https://stackoverflow.com/questions/13499852/scipy-fourier-transform-of-a-few-selected-frequencies
// https://stackoverflow.com/questions/11579367/implementation-of-goertzel-algorithm-in-c
//...if we go down that route, we may also use the obtained amplitude
// for refining our amplitude estimate...ok - for now, just leave the instantaneous phase 
// measurement at zero

https://ccrma.stanford.edu/~jos/sasp/Phase_Interpolation_Peak.html


see here for window functions:
https://en.wikipedia.org/wiki/Window_function
implement blackman-harris and maybe blackman-nutall, dolph-chebychev - i think, low sidelobes are
important for sinusoidal parameter estimation...maybe gaussian for frequency estimation? or use 
re-assignment?
*/


/*
Ideas for transformations:
-insert phase anomalies into partials - maybe at regular intervals, offset the phase by a certain
 dp (delta-phi) - this should impose some sort of periodic pattern, maybe a sort of robotization
 ...but maybe that time-interval of the phase anomalies could be different for different partials
 "Phase-Bumper"


*/