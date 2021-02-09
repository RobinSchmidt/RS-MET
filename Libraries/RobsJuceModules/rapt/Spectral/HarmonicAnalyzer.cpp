
template<class T>
rsHarmonicAnalyzer<T>::rsHarmonicAnalyzer() : cycleFinder(sampleRate)
{
  typedef rsFourierTransformerRadix2<T> FT;
  trafo.setDirection(FT::FORWARD);
  trafo.setNormalizationMode(FT::NORMALIZE_ON_FORWARD_TRAFO);
}

template<class T>
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N)
{
  RAPT::rsSinusoidalModel<T> mdl;   // init empty model
  if(flattenPitch(x, N) == false)   // pre-process audio (flatten pitch), sets up the blockSize
    return mdl;                     // return empty model if pre-processing has failed

  if(useOldCode)
    analyzeHarmonicsOld(mdl);       // todo: remove this switch eventually
  else
    analyzeHarmonics(mdl);          // create model from pitch-flattened signal (now in member y)

  deFlattenPitch(mdl);              // post process model data to account for flattening
  if(antiAlias)
    removeAliasing(mdl);            // remove freqs above orignal nyquist freq
  handleEdges(mdl);                 // add fade-in/out datapoints
  convertTimeUnit(mdl);             // convert from samples to seconds
  refineFrequencies(mdl);           // refines freq estimates, if desired

  //rosic::writeToMonoWaveFile("PitchFlattened.wav", &y[0], (int)y.size(), (int)sampleRate);
  // move to rapt - rapt is a lower layer than rosic an we are not supposed to call rosic functions
  // inside rapt functions...maybe rename to rsWavWrite

  return mdl;

  // todo: clean up the model: remove partials that are consistently above the nyquist freq - due 
  // to the stretching, we may get frequencies almost up to the sample-rate (we downshift at most 
  // by an octave (with respect to the longest original cycle) - so the model may end up having 
  // twice the number of frequencies as it should have...(although their amplitudes are really low)
  // ...maybe also introduce an amplitude threshold and remove partials that are consistently below
  // that threshold

  // maybe to fix the spurious on/off switching of certain harmonics, we should "repair" the 
  // partials that show this effect - detect the gaps and fill them by interpolating. ...but that's 
  // somehow cheating and should be used only as last resort - first, we should try to refine the 
  // decision criteria in order to avoid such a situation in the first place - try to cook up 
  // signals that show this artifact and investigate, how to adjust the decision criteria

  // Maybe during analysis, we should create a "report" that the user may inquire after an analysis
  // containing information about what may have went wrong (for example, the number of detected
  // gaps in the repair, etc.)
}

template<class T>
bool rsHarmonicAnalyzer<T>::flattenPitch(T* x, int Nx)
{
  typedef RAPT::rsArrayTools AR;
  typedef std::vector<T> Vec;

  // Find cycle marks and assign FFT blockSize:
  Vec cycleMarks = findCycleMarks(x, Nx);        // cycle marks
  //rsPlotSignalWithMarkers(x, Nx, &cycleMarks[0], (int) cycleMarks.size());

  rsAssert(cycleMarks.size() >= 2);              // something went wrong in the cycel mark finder
  if(cycleMarks.size() < 2) return false;        // report failure

  Vec cycleLengths = rsDifference(cycleMarks);   // cycle lengths
  T maxLength = rsMaxValue(cycleLengths);
  //maxLength   = rsMax(maxLength, cycleMarks[0]);             // delta between 0 and 1st mark
  //maxLength   = rsMax(maxLength, (Nx-1)-rsLast(cycleMarks)); // delta between end and last mark
  rsAssert(maxLength >= 2);                      // at least 2 samples per cycle
  if(maxLength < 2) return false;                // report failure

  // determine analysis cycle length, i.e. the fixed length to which all cycles will be stretched:
  setMaxMeasuredCycleLength(maxLength); // does also some buffer-re-allocation 
  //rsPlotVector(cycleLengths);
  //rsPlotVector(sampleRate/cycleLengths);

  // Create the mapping function for the time instants of the cycle marks:
  int mapLength = (int) cycleMarks.size() + 2;  // +2 for t = 0 and t = N-1
  tIn.resize( mapLength);
  tOut.resize(mapLength);
  tIn[0] = tOut[0] = 0;    // time-origin is zero for both, original and stretched signal

  // The first marker is mapped to an instant, such that the initial partial cycle is stretched by
  // the same amount as the first full cycle (the cycle between 1st and 2nd marker):
  tIn[1]  = cycleMarks[0];
  tOut[1] = cycleMarks[0] * cycleLength / cycleLengths[0];
  tOut[1] = round(tOut[1]);

  // All cycles between the initial partial cycle and final partial cycle are stretched to the same 
  // fixed length:
  for(int i = 2; i < mapLength-1; i++) {
    tIn[i]  = cycleMarks[i-1];
    tOut[i] = tOut[i-1] + cycleLength; }

  // The end time instant is mapped such that the final partial cycle is stretched by the same 
  // amount as the last full cycle:
  T tailLength = (Nx-1) - rsLast(cycleMarks);
  tIn [mapLength-1] = Nx-1;
  tOut[mapLength-1] = tOut[mapLength-2] + tailLength * cycleLength / rsLast(cycleLengths);
  tOut[mapLength-1] = round(tOut[mapLength-1]);

  // We have created the time warping map, sampled at the cycle-marks. For applying the 
  // warping, we interpolate it up to sample rate via linear interpolation:
  int Ny = (int) rsLast(tOut) + 1; // length of stretched signal and warping map
  Vec t(Ny), w(Ny);                // interpolated time axis and warping map
  AR::fillWithIndex(&t[0], Ny);
  RAPT::resampleNonUniformLinear(&tOut[0], &tIn[0], mapLength, &t[0], &w[0], Ny);

  // Do the time-warping:
  y.resize(Ny);
  rsTimeWarper<T, T>::timeWarpSinc(x, Nx, &y[0], &w[0], Ny, sincLength);
  return true;
}

template<class T>
void rsHarmonicAnalyzer<T>::setMaxMeasuredCycleLength(T maxLength)
{
  // old - allows ony powers of two for cyclesPerBlock:
  cycleLength = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));
  blockSize   = cyclesPerBlock * cycleLength;

  // new - allows arbitrary values for cyclesPerBlock (should be integer, though):
  // ...something to do...
  // blockSize = nextPowerOfTwo(cyclesPerBlock * maxLength);  // should work for float inputs(?)
  // cycleLength = T(blockSize) / T(cyclesPerBlock);


  sig.resize(blockSize); 
  wnd.resize(blockSize);
  fillWindow();
  trafoSize = zeroPad * blockSize;
  sigPadded.resize(trafoSize);
  mag.resize(trafoSize); 
  phs.resize(trafoSize);
  trafo.setBlockSize(trafoSize);
}

// can this be deleted sometime soon?
template<class T>
void rsHarmonicAnalyzer<T>::analyzeHarmonicsOld(RAPT::rsSinusoidalModel<T>& mdl)
{
  // maybe, if we want to keep the simpler algo around, we should assert that numCyclesPÜerBlock==1
  // and the window is rectangular

  // Initialize the model (create all datapoints, to filled with actual data later):
  mdl.init(getNumHarmonics(), getNumDataPoints());

  // The initial partial cycle is pre-padded with zeros:
  int n0 = 0;                          // first sample (from y-array) in current frame
  int m  = 0;                          // frame index
  int K  = blockSize;      // this is now potentially different from cycle-length...we need to upadate code below
  int L  = (int) tOut[1];              // length of initial partial cycle
  rsAssert(L >= 0 && L <= K);
  typedef RAPT::rsArrayTools AR;
  AR::fillWithZeros(&sig[0], K-L);
  if(L > 0)
    AR::copy(&y[n0], &sig[K-L], L);
  //rsPlotVector(sig);
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));

  // The inner cycles/frames are taken as is:
  L = K;                               // length of inner cycles
  for(m = 1; m < getNumFrames()-1; m++) {
    n0 = (int) tOut[m];                   // ...why not round?
    AR::copy(&y[n0], &sig[0], L);

    //rsPlotVector(sig);
    //// plot 2nd-to-last (debug):
    //if(m == getNumFrames()-2)
    //  plotVector(sig);

    fillHarmonicData(mdl, m, getTimeStampForFrame(m));
  }

  // The final partial cycle is post-padded with zeros:
  n0 = (int) tOut[m]; 
  L = int(tOut[tOut.size()-1] - tOut[tOut.size()-2]) + 1;
  rsAssert(L >= 0 && L <= K);
  AR::copy(&y[n0], &sig[0], L);
  if(L < K)  // is this correct?
    AR::fillWithZeros(&sig[L], K-L);
  //rsPlotVector(sig);  // debug
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));


  // todo: double-check all index computations against off-by-one errors, verify time-indices
}

template<class T>
void rsHarmonicAnalyzer<T>::analyzeHarmonics(RAPT::rsSinusoidalModel<T>& mdl)
{
  // Initialize the model (create all datapoints, to filled with actual data later):
  mdl.init(getNumHarmonics(), getNumDataPoints());
  int over = (blockSize - cycleLength) / 2; // amount of overhanging of block with respect to cycle
  for(int m = 0; m < getNumFrames(); m++)
  {
    int cycleStart = (int) tOut[m];
    int cycleEnd   = (int) tOut[m+1];  // safe: tOut.size() == getnumFrames()+1
    int blockStart = cycleStart - over;
    int blockEnd   = cycleEnd   + over;

    // adjustments for the first and last (improper) cycle:
    int length = blockEnd-blockStart; 
    if(length != blockSize) {
      rsAssert(m == 0 || m == getNumFrames()-1); // should only happen in first or last frame
      int delta = blockSize - length;
      if(m == 0) blockStart -= delta;            // check if this is correct
      else       blockEnd += delta;
      length = blockEnd-blockStart; }            // only relevant for debugging
    rsAssert(blockEnd-blockStart == blockSize);

    // copy section from y into sig, apply window and extract spectral data:
    rsArrayTools::copySection(&y[0], (int) y.size(), &sig[0], blockStart, blockSize);
    for(size_t n = 0; n < sig.size(); n++)
      sig[n] *= wnd[n];
    //rsPlotVector(sig);
    fillHarmonicData(mdl, m, getTimeStampForFrame(m));
      // maybe we should pass the "delta" from above and use it to adjust the phase values like
      // phase += delta*omega or something
  }
}



template<class T>
void rsHarmonicAnalyzer<T>::deFlattenPitch(RAPT::rsSinusoidalModel<T>& mdl)
{
  // todo: get rid of these vectors - compute values on the fly inside the loop
  std::vector<T> lw = rsDifference(tOut);  // warped lengths of cycles
  std::vector<T> lu = rsDifference(tIn);   // unwarped lengths of cycles
  for(int m = 0; m < getNumFrames(); m++) {
    T tw = getTimeStampForFrame(m);         // warped time
    T tu = getUnWarpedTimeStampForFrame(m); // unwarped time
    T r  = lw[m] / lu[m];                   // stretching ratio applied to frame m
    for(int k = 0; k < getNumHarmonics(); k++) {
      mdl.getDataRef(k, m+1).time  = tu;    // m+1 because datapoint-index is frame-index + 1
      mdl.getDataRef(k, m+1).freq *= r;     // due to initial "fade-in" datapoint at time zero
    }
  }
}

template<class T>
void rsHarmonicAnalyzer<T>::removeAliasing(RAPT::rsSinusoidalModel<T>& mdl)
{
  // maybe we should have an option to select, if only those partials are removed that are above
  // sampleRate/2 all the time or also those that exceed the Nyquist limit only temporarily
  // maybe call the option setAntAliasMode() options: based-on-min-freq, based-on-max-freq

  bool allTheTime = false;   // make user option

  // Find the index of the partial, above which all higher ones may be removed (todo: use binary 
  // instead of linear search):
  int maxIndexToRetain = getNumHarmonics() - 1;
  while(maxIndexToRetain > 0) {
    if( mdl.getPartial(maxIndexToRetain).willAlias(sampleRate, allTheTime) )
      maxIndexToRetain--;
    else
      break;
  }

  // ...and remove all partials above the found index:
  mdl.removePartialsAbove(maxIndexToRetain);
}

template<class T>
void rsHarmonicAnalyzer<T>::handleEdges(RAPT::rsSinusoidalModel<T>& mdl)
{
  // now we must fill in the data at the very first and very last datapoint index to get a 
  // fade-in/out:
  int k;
  int lastDataIndex = (int) mdl.getPartial(0).getNumDataPoints()-1;
  int i = lastDataIndex-1;  // where we read the freq and phase from
  //T endTime = rsLast(tOut);
  T endTime = rsLast(tIn);
  rsInstantaneousSineParams<T> params;
  for(k = 0; k < (int) mdl.getNumPartials(); k++) 
  {
    // fill first datapoint:
    params  = mdl.getPartial(k).getDataPoint(1);  
    T dt    = (T(0) - params.time) / sampleRate;  // time delta between 2nd and 1st datapoint
    T freq  = params.freq;
    T phase = params.phase + 2*PI*freq * dt;
    phase   = rsWrapToInterval(phase, -PI, PI);
    mdl.setData(k, 0, T(0), freq, T(0), phase);  

    // fill last datapoint:
    params = mdl.getPartial(k).getDataPoint(i);
    dt     = (endTime - params.time) / sampleRate; // time delta between 2nd-to-last and last datapoint
    freq   = params.freq;
    phase  = params.phase + 2*PI*freq * dt;
    phase  = rsWrapToInterval(phase, -PI, PI);
    mdl.setData(k, i+1, endTime, freq, T(0), phase);
    //mdl.setData(k, i+1, endTime, freq, params.gain, phase);  // test - nope - may make end artifacts worse
    int dummy = 0;
    // todo: check, if last datapoint is correct - there is an artifact in the resynthesized signal
    // at the end
  }

  int dummy = 0;

  // the distance of the very first marker from the time origin t=0 should probably used for 
  // determining the start phase - don't assume an additional "ghost" marker at t=0 - instead, let
  // the sinusoid start at zero amplitude, frequency determined by the distance between 1st and 2nd
  // marker and phase appropriate to the frequency and time-value of the 1st marker (i.e. if the 
  // first marker is at 25 and the second is at 125, assume a cycle length of 100 and start phase
  // of -90° (a quarter period) - for higher harmonics, take into account the phase-measurement
  // at first marker (for the fundamental, that phase is zero by construction ...but only, if we 
  // use f0 zero-crossings for the cycle-mark finder...hmmm....
}

template<class T>
void rsHarmonicAnalyzer<T>::convertTimeUnit(RAPT::rsSinusoidalModel<T>& mdl)
{
  for(size_t hi = 0; hi < mdl.getNumPartials(); hi++)
    for(int di = 0; di < getNumDataPoints(); di++)
      mdl.getDataRef((int)hi, di).time /= sampleRate;
}

template<class T>
void rsHarmonicAnalyzer<T>::refineFrequencies(RAPT::rsSinusoidalModel<T>& mdl)
{
  if(freqsByPhaseDerivative)
    rsSinusoidalProcessor<T>::refineFreqsViaPhaseDerivative(mdl);

  if(freqsPhaseConsistent)
    mdl.makeFreqsConsistentWithPhases(); // use rsSinusoidalProcessor
}

template<class T>
std::vector<T> rsHarmonicAnalyzer<T>::findCycleMarks(T* x, int N)
{
  // To ensure that initial and final section are really partial cycles (as opposed to a full
  // cycle plus something extra), we prepend and/or append artificial cycle marks in these cases.
  // The positions of the artificial marks are set from the length of the first or last cycle
  // respectively:
  std::vector<T> cm = cycleFinder.findCycleMarks(x, N);
  if(cm.size() >= 2) {
    T L = cm[1] - cm[0];
    while(cm[0] > L)
      rsPrepend(cm, cm[0]-L);
    L = cm[cm.size()-1] - cm[cm.size()-2];
    while((N-1) - cm[cm.size()-1] > L)
      rsAppend(cm, cm[cm.size()-1]+L);
  }
  //rsPlotSignalWithMarkers(x, N, &cm[0], (int) cm.size());
  return cm;
}

template<class T>
T rsHarmonicAnalyzer<T>::getTimeStampForFrame(int m)
{
  // maybe use just tOut[m], depending on a user setting which determines whether the time origin
  // for each frame should be at start or center of the frame - this setting should then also 
  // affect, whether or not we shift the signal buffer by half of its length before FFT)

  //return sampleRate * (tOut[m] + T(0.5)*(tOut[m+1]-tOut[m]));

  // hmm...maybe it's not a good idea to include the sampleRate factor at this stage for numerical
  // reasons - maybe keep the time in samples at this stage - eventually, the model wants its
  // time stamps in seconds - but we probably should include the sample-rate factor at the very end
  // of the algo

  return tOut[m] + T(0.5)*(tOut[m+1]-tOut[m]);
}

template<class T>
T rsHarmonicAnalyzer<T>::getUnWarpedTimeStampForFrame(int m)
{
  return tIn[m] + T(0.5)*(tIn[m+1]-tIn[m]);
}

template<class T>
int rsHarmonicAnalyzer<T>::getSpectralPeakSearchWidth()
{
  //T peakSearchWidth = T(1); // maybe make user parameter later
  //T mainlobeWidth = rsWindowFunction::getMainLobeWidth(windowType, T(0));  // old
  //T mainlobeWidth = rsWindowFunction::getMainLobeWidth(windowType, getWindowParameter());
  T mainlobeWidth = rsWindowFunction::getMainLobeWidth(
    windowType, getWindowParameter(), blockSize);

  //T s = T(0.5);  
  // old: 1.0 - works for hamming, blackman needs around 0.5 (numCycles=4) - when the mainlobe is
  // too wide, the actual maximum may occur at the border of the window and then the partial gets
  // discarded - (i think) we may remedy this by using more cycles per block - but the next larger
  // value is 8 which is a bit large -> todo: allow for non-powers-of-two - maybe the search width
  // should be some function of the mainlobe-width and number-of-cycles (and maybe a user 
  // parameter)

  return (int) round(T(0.5)*peakSearchWidth*zeroPad*mainlobeWidth);
  // 0.5, because we search to boht sides by the distance, so this function actually return half
  // of the serach-width
  

  //return (int) ceil(T(0.5)*peakSearchWidth*zeroPad*mainlobeWidth) + 1;
  //return (int) ceil(T(0.5)*peakSearchWidth*zeroPad*mainlobeWidth) - 1;

  // maybe this peakserachWidth should also be related to the minPeakWidth - if it's too wide (!), 
  // we may miss peaks because the max-values happen at the end of the search interval

  //return (int) round(T(0.5)*peakSearchWidth*zeroPad*mainlobeWidth);
}

template<class T>
void rsHarmonicAnalyzer<T>::fillHarmonicData(
  RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T time)
{
  //rsPlotVector(sig);
  prepareBuffer(sig, sigPadded);
  trafo.getRealSignalMagnitudesAndPhases(&sigPadded[0], &mag[0], &phs[0]);  // perform FFT

  //if(frameIndex == 231 || frameIndex == 232 )  // 2nd harmonic in flute switches from off to on (bm-window)
  //  rsPlotSpectrum(mag, T(0), T(-150), true);

  //if(frameIndex == getNumFrames()/2)
  //  rsPlotSpectrum(mag, T(0), T(-150), true); // freq axis wrong, if we pass the sampleRate

  //if(frameIndex >= 10)
  //  rsPlotSpectrum(mag, sampleRate, T(-200));

  // extract model data from FFT result:
  int dataIndex   = frameIndex + 1;        // +1 because of the fade-in datapoint
  int numBins     = trafo.getBlockSize();  // number of FFT bins
  int numPartials = getNumHarmonics();     // number of (pseudo) harmonics
  int kHarm;                               // bin index where partial/harmonic is expected

  if(numPartials == 0)
    return; 
    // fix crash - but probably, a situation like this should be avoided at a higher level         
    // of the algo - i.e. make sure that all cycles have a length >= minLength or something

  if(zeroPad == 1 && cyclesPerBlock == 1) { 
    // old version (before multi-cycle and zero-padding support):
    for(kHarm = 0; kHarm < numPartials; kHarm++) {
      T freq = trafo.binIndexToFrequency(kHarm, numBins, sampleRate);
      mdl.setData(kHarm, dataIndex, time, freq, T(2)*mag[kHarm], phs[kHarm]);
    }
  }
  else {
    // new, general version - supporting multiple cycles per block and/or zero-padding:
    int kPeak;                             // bin index where partial is actually found
    int kPeakOld = -1;                     // kPeak from previous iteration
    int w2 = getSpectralPeakSearchWidth(); // we look for a peak in the range: kHarm +- w2
    T freq, gain, phase, peakBin;


    // handle DC separately:  
    //mdl.setData(0, dataIndex, time, T(0), T(2*zeroPad)*mag[0], phs[0]); // this was wrong

    // this is still unfinished!
    mdl.setData(0, dataIndex, time, T(0), T(zeroPad)*mag[0], phs[0]); // handle DC separately
    phs[0] = T(0); // is this correct? the FFT analyzer uses the convention to put the Nyquist-freq
                   // analysis value there...is that actually correct? check that out!
      // the handling of DC is not yet good - we may produce negative amplitudes here - but maybe, 
      // we should allow negative amplitudes generally in the sinusoidal model?


    for(int h = 1; h < numPartials; h++) {

      kHarm = cyclesPerBlock*zeroPad*h;  // bin index where partial/harmonic is expected
      //kHarm = getPartialBinIndex();        // may be used when we use minPartialIndex

      kPeak = kHarm;                       // preliminary - refined below
      //kPeak = -1;                          // preliminary - assigned below

      if(allowInharmonics)                 // search for peaks near expected harmonics
      {
        kPeak = findPeakBinNear(mag, kHarm, w2);
        if(kPeak == -1 || kPeak == kPeakOld) {
          // no peak found or the same peak was found a second time (the latter case may occur, if 
          // there's a partial halfway between expected harmonic frequencies - then we take only
          // the first one seriously and discard the second) - use kHarm for the frequency, zero 
          // for the amplitude:
          freq = trafo.binIndexToFrequency(kHarm, numBins, sampleRate);
          gain = T(0);
          phase = phs[kHarm];
        } 
        else {
          if(parabolicInterpolation) {
            rsSinusoidalAnalyzer<T>::spectralMaximumPositionAndValue(
              &mag[0], kPeak, &peakBin, &gain);
            freq  = peakBin * sampleRate / numBins;
            gain *= T(2*zeroPad);
            if(phaseInterpolation)
              phase = rsSinusoidalAnalyzer<T>::interpolatePhase(&phs[0], peakBin);
            else
              phase = phs[kPeak];
          } 
          else {
            freq  = trafo.binIndexToFrequency(kPeak, numBins, sampleRate);
            gain  = T(2*zeroPad)*mag[kPeak];
            phase = phs[kPeak];
          }
        }
        mdl.setData(h, dataIndex, time, freq, gain, phase);
      }
      else {
        freq = trafo.binIndexToFrequency(kHarm, numBins, sampleRate);
        gain = T(2*zeroPad)*mag[kHarm]; 
        mdl.setData(h, dataIndex, time, freq, gain, phs[kHarm]);
      }
      kPeakOld = kPeak; // remember it to make sure to not record the same peak twice
    }
  }
}

/*
template<class T>
std::vector<int> rsHarmonicAnalyzer<T>::findPartialBins(const std::vector<T> mag)
{
  std::vector<int> peakBins;
  int w2 = getSpectralPeakSearchWidth();  // search width to the left and right
  int numPartials = getNumHarmonics();    // number of (pseudo) harmonics

  return peakBins;
}
*/

// move to rsArrayTools:
template<class T>
int numPeaks(T* x, int N)
{
  int np = 0;
  for(int i = 1; i < N-1; i++)
    if(x[i] >= x[i-1] && x[i] >= x[i+1]) // with >=, this finds the plateaus, too - maybe rename function
      np++;
  return np;
}

template<class T>
int rsHarmonicAnalyzer<T>::findPeakBinNear(std::vector<T>& v, int kCenter, int w2)
{
  //return kCenter; // test
  // maybe always return a valid index and an success/error-code in a second output

  //bool dbgIsHarmonic = kCenter == 32 || kCenter == 960 || kCenter == 992;
  //// for the two sines at 200Hz/6100Hz

  if(w2 == 0)
    return kCenter;
  int kLeft  = rsMax(kCenter - w2, 0);
  int kRight = rsMin(kCenter + w2, (int) v.size()-1);
  int length = kRight - kLeft + 1;

  //if(dbgIsHarmonic)
  //{
  //  rsPlotSpectrum(toVector(&v[kLeft], length), 0.0, -120.0, true); // for decibels
  //  //rsPlotArray(&v[kLeft], length); // plot segment where we search for a peak
  //}

  if(w2 == 1) {
    if(v[kCenter] >= v[kLeft] && v[kCenter] >= v[kRight])
      return kCenter;
    else
      return -1;
  }

  // If there's a partial near kCenter, we expect to see an unimodal distribution of values inside 
  // the search window because our search window equals the mainlobe width. This way, we may 
  // distinguish mainlobes from sidelobes - sidelobes are narrower which leads to multiple local 
  // maxima in the search window.
  // ..hmm...but what if the maximum
  // is shifted to the side - shouldn't we use *half* of the mainlobe-width? make plots and check
  //int nPeaks = numPeaks(&v[kLeft], length);
  //if(nPeaks != 1)
  //  return -1;
  // maybe make enforcing this unimodality condition optional


  int kMax = rsArrayTools::maxIndex(&v[kLeft], length) + kLeft;  // index of maximum
  //return kMax; // test - crashes bcs parabolic peak interpolation returns out-of range value
  if(kMax == kLeft || kMax == kRight) // *not* ensured already by nPeaks == 1: there could be a
    return -1;                        // bump in the middle but the side could still be higher
  else {
    if(isPeakPartial(v, kMax))
      return kMax;
    else
      return -1;
  }

  // todo: maybe impose an additional threshold constraint - maybe relative to the absolute maximum
  // of the whole spectrum
  // maybe, if a local maximum is found, find the two local minima that surround it and check if 
  // their distance is >= peak-search width

  // maybe, we should not just return -1 in the case when we don't find a peak that meets all 
  // constraints but instead return some negative number that indicates what exactly happened (i.e.
  // which test failed) - and then instead of recording a zero amplitude at kCenter, preliminarily
  // record that negative number and later clean these values up in a post-processing step
  // ...and in cases of short gaps, fill them by interpolation from the surrounding data
  // ...maybe also don't allow peaks to instantly switch on/off from one datapoint to the next

  // maybe let the user switch between various decision strategies - also include plain old
  // amplitude thresholds
}

// maybe rename to: isPeakWideEnough
template<class T>
bool rsHarmonicAnalyzer<T>::isPeakPartial(std::vector<T>& v, int peakBin)
{
  //return true;   // test

  // minWidth can be precomputed:
  //T mainlobeWidth = zeroPad * rsWindowFunction::getMainLobeWidth(windowType, T(0)); // old
  T mainlobeWidth = zeroPad * rsWindowFunction::getMainLobeWidth(
    windowType, getWindowParameter(), blockSize);

  T minWidth1     = minPeakToMainlobeWidthRatio*mainlobeWidth;

  T harmonicWidth = zeroPad * cyclesPerBlock;
  T minWidth2     = minPeakToHarmonicWidthRatio*harmonicWidth;

  int minWidth    = (int) round( rsMin(minWidth1, minWidth2) );
    // ...this needs some serious consideration - if it's really a good idea to do it like this...





  // find indices of local minima that surround the local maximum at peakBin and the width between
  // these two local minima:
  int kR = peakBin + 1;             // index of right local minimum
  while(kR <= (int)v.size()-2) {
    if(v[kR] <= v[kR-1] && v[kR] <= v[kR+1])
      break;
    if(kR - peakBin - 1 > minWidth)   // return early if the left minimum is very faar away
      return true;
    kR++;
  }
  int kL = peakBin - 1;
  while(kL >= 1) {
    if(v[kL] <= v[kL-1] && v[kL] <= v[kL+1])
      break;
    if(peakBin - kL - 1 > minWidth)
      return true;
    kL--;
  }
  int width = kR - kL + 1;

  // compare the computed/measured lobe width against a minimum allowed width - the width must be 
  // larger than that in order to consider the peak a partial/mainlobe:
  if(width >= minWidth)
    return true;
  else
    return false;
}


template<class T>
void rsHarmonicAnalyzer<T>::prepareBuffer(const std::vector<T>& sig, std::vector<T>& buf)
{
  size_t K2 = sig.size() / 2;
  size_t M  = buf.size();
  size_t i;
  for(i = 0;  i < K2;   i++) buf[i] = sig[i+K2];     // first section is 2nd half of sig
  for(i = K2; i < M-K2; i++) buf[i] = 0;             // middle section is zero padding
  for(i = 0;  i < K2;   i++) buf[M-K2+i] = sig[i];   // last section is 1st half of sig
  // It may seem, that just swapping left and right half of the buffer would lead to phase 
  // measurements that are off by half a sample because the center of an even-length buffer falls
  // on a half-integer - but: our datapoints are actually also placed at the half-integers, so in 
  // the end, it works out correctly.
}

template<class T>
void rsHarmonicAnalyzer<T>::fillWindow()
{
  rsWindowFunction::createWindow(
    &wnd[0], (int) wnd.size(), windowType, true, getWindowParameter());
}


/*

Ideas:

Currently, we stretch every cycle of the input signal to some fixed cycleLength (the power of two
greater or equal to the length of the longest cycle in the input). The number nc of cycles per 
block must be a power of two because we want the blockSize to be a power of two for the FFT 
friendliness.. However, having nc to be restricted to powers of two is quite restrictive, so here
is an idea to lift that restriction while still using power-of-2 blocksizes: Instead of stretching 
every cycle separately to a given cycleLength, we may stretch groups of cycles to the desired 
blockSize. Consider an example that has cycle lengths like:

  98,101,102,99,97,100,102,101,99,97,100,102

the current implementation would stretch all of them indiviudally to a length of 128. If nc=2, we 
would take 2 of these length 128 cycles for every block. Instead we, could group the cycles in 
pairs:

  98,101, 102,99,  97,100,  102,101,  99,97,  100,102
    199     201     197       203      196      202    -> stretch all these pairs to length 256

and stretch every such pair of cycles to length 256. The advantage is that it readily generalizes 
to using any number of cycles for such a group and using as blocksize the power of two that is 
greater or equal to the length of the longest group. For example, with nc=3:

  98,101,102, 99,97,100, 102,101,99, 97,100,102
     301         296         302        299            -> stretch all these triples to length 512

...wait - the pitch-flattening is not restricted to stretch every cycle to an integer number of 
samples. If we want to fit 3 cycles into a block of 512 samples, we could just stretch every cycle 
to 512/3 = 170.666... samples. I think, the required changes are in flattenPitch:

  cycleLength = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));
  setCycleLength(cycleLength); 

should be replaced by:

  blockSize = nextPowerOfTwo(cyclesPerBlock * maxLength);  // should work for float inputs(?)
  cycleLength = T(blockSize) / T(cyclesPerBlock);
  setCycleLength(cycleLength); 

where the type of the member cycleLength should be changed from int to T. setCycleLength will 
(again) set the blockSize - that line should probably be removed then. Or we do all assignments in
setCycleLength which should then take the measured maxLength as parameter. Maybe a (protected) 
function setMaxMeasuredCycleLength should be introduced which does everything....
Maybe cyclesPerBlock could then be type T, too - i.e. non-integer. i don't know if it's useful to 
use a non-integer number of cycles per block - probably not, but maybe for strongly inharmonic 
signals, but we could experiment with that...


Do temporal oversampling: take datapoints not only every cycle but every half-cycle or 
quarter-cycle. It was observed that the error signal is largest halfway between the datapoints, 
which is not surprising. So, placing new datapoints at exactly those positions where the error is
currently greatest, should effectively reduce the error signal.


Another idea to estimate the instantaneous frequncies: Do the same block based measurements as now 
additionally with one sample offset to the left and to the right and estimate the frequency by 
numerically differentiating the instantaneous phase: f = dp/dt ~= (p(t+dt) - p(t-dt)) / (2*dt) 
where dt is the sampling period. We could also estimate the 2nd derivative of the phase and store 
that in an additional data field (called "glide" or "slide" or something) that defaults to zero. 
That could be used as target value for quintic phase interpolation. Maybe one sample offset is too
little - maybe experiment with 2,3, etc. maybe up to half a cycle. There are some questions how 
exactly to incorporate the frequency information that is already available by the current method 
(and similar information that will be available at the two additional datapoints)...maybe some sort
of "best-fit" strategy (whatever that means in this context) could make sense?



ToDo: 

Figure out, if there's any advantage in any circumstances to use Blackman or Hamming or other 
windows as opposed to the Dolph-Chebychev - i think, the Dolph-Cheby window is actually optimal for 
this kind of sinusoidal analysis, so we may not need any others anymore and could throw them out at
some point to clean up the interface...but who knows - for the time being, i'll leave them in, just 
in case...maybe at a later stage, the usage of other windows can be retained in an experimental 
subclass somewhere outside the main library


*/