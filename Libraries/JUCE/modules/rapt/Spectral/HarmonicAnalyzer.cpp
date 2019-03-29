
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
  if(flattenPitch(x, N) == false)   // pre-process audio (flatten pitch), sets the blockSize
    return mdl;                     // return empty model if pre-processing has failed

  // todo: remove this switch eventually:
  if(useOldCode)
    analyzeHarmonics(mdl);            // create model from pitch-flattened signal (now in member y)
  else
    analyzeHarmonics2(mdl);

  deFlattenPitch(mdl);              // post process model data to account for flattening
  if(antiAlias)
    removeAliasing(mdl);            // remove freqs above orignal nyquist freq
  handleEdges(mdl);                 // add fade-in/out datapoints
  convertTimeUnit(mdl);             // convert from samples to seconds
  refineFrequencies(mdl);           // refines freq estimates, if desired

  rosic::writeToMonoWaveFile("PitchFlattened.wav", &y[0], (int)y.size(), (int)sampleRate);
  // move to rapt - rapt is a lower layer than rosic an we are not supposed to call rosic functions
  // inside rapt functions...maybe rename to rsWavWrite

  return mdl;

  // todo: clean up the model: remove partials that are consistently above the nyquist freq - due 
  // to the stretching, we may get frequencies almost up to the sample-rate (we downshift at most 
  // by an octave (with respect to the longest original cycle) - so the model may end up having 
  // twice the number of frequencies as it should have...(although their amplitudes are really low)
  // ...maybe also introduce an amplitude threshold and remove partials that are consistently below
  // that threshold
}

template<class T>
bool rsHarmonicAnalyzer<T>::flattenPitch(T* x, int Nx)
{
  typedef RAPT::rsArray AR;
  typedef std::vector<T> Vec;

  // Find cycle marks and assign FFT blockSize:
  Vec cycleMarks = findCycleMarks(x, Nx);        // cycle marks
  //plotSignalWithMarkers(x, Nx, &cycleMarks[0], (int) cycleMarks.size());
  if(cycleMarks.size() < 2)
    return false;                                // report failure
  Vec cycleLengths = rsDifference(cycleMarks);   // cycle lengths
  T maxLength = rsMax(cycleLengths);
  //maxLength   = rsMax(maxLength, cycleMarks[0]);             // delta between 0 and 1st mark
  //maxLength   = rsMax(maxLength, (Nx-1)-rsLast(cycleMarks)); // delta between end and last mark
  cycleLength = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));
  setCycleLength(cycleLength);  // does also some buffer-re-allocation
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
    tOut[i] = tOut[i-1] + cycleLength;
  }

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
void rsHarmonicAnalyzer<T>::analyzeHarmonics(RAPT::rsSinusoidalModel<T>& mdl)
{
  // Initialize the model (create all datapoints, to filled with actual data later):
  mdl.init(getNumHarmonics(), getNumDataPoints());

  // The initial partial cycle is pre-padded with zeros:
  int n0 = 0;                          // first sample (from y-array) in current frame
  int m  = 0;                          // frame index
  int K  = blockSize;      // this is now potentially different from cycle-length...we need to upadate code below
  int L  = (int) tOut[1];              // length of initial partial cycle
  rsAssert(L >= 0 && L <= K);
  typedef RAPT::rsArray AR;
  AR::fillWithZeros(&sig[0], K-L);
  if(L > 0)
    AR::copyBuffer(&y[n0], &sig[K-L], L);
  //rsPlotVector(sig);
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));

  // The inner cycles/frames are taken as is:
  L = K;                               // length of inner cycles
  for(m = 1; m < getNumFrames()-1; m++) {
    n0 = (int) tOut[m];                   // ...why not round?
    AR::copyBuffer(&y[n0], &sig[0], L);

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
  AR::copyBuffer(&y[n0], &sig[0], L);
  if(L < K)  // is this correct?
    AR::fillWithZeros(&sig[L], K-L);
  //rsPlotVector(sig);  // debug
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));


  // todo: double-check all index computations against off-by-one errors, verify time-indices
}

template<class T>
void rsHarmonicAnalyzer<T>::analyzeHarmonics2(RAPT::rsSinusoidalModel<T>& mdl)
{
  // Initialize the model (create all datapoints, to filled with actual data later):
  mdl.init(getNumHarmonics(), getNumDataPoints());

  typedef RAPT::rsArray AR;

  //int numFrames = getNumFrames();  //
  int over = (blockSize - cycleLength) / 2; // amount of overhanging of block with respect to cycle

  for(int m = 0; m < getNumFrames(); m++)
  {
    int cycleStart = (int) tOut[m];
    int cycleEnd   = (int) tOut[m+1];  // safe: tOut.size() == getnumFrames()+1
    int blockStart = cycleStart - over;
    int blockEnd   = cycleEnd   + over;

    // but for the first and last improper cycle, this is wrong...
    int length = blockEnd-blockStart; 
    if(length != blockSize)
    {
      // ...we must do something extra in these special cases...
      rsAssert(m == 0 || m == getNumFrames()-1); // should only happen in first or last frame
      int delta = blockSize - length;

      // check if this is correct:
      if(m == 0)
        blockStart -= delta;
      else
        blockEnd += delta;

      length = blockEnd-blockStart;  // update - only relevant for debug
    }
    rsAssert(blockEnd-blockStart == blockSize);

    // copy section from y into sig, apply window and extract spectral data:
    AR::copySection(&y[0], (int) y.size(), &sig[0], blockStart, blockSize);
    for(size_t n = 0; n < sig.size(); n++)
      sig[n] *= wnd[n];
    //rsPlotVector(sig);
    fillHarmonicData(mdl, m, getTimeStampForFrame(m));
      // maybe we should pass the "delta" from above and use it to adjust the phase values like
      // phase += delta*omega or something

    int dummy = 0;
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
  for(k = 0; k < mdl.getNumPartials(); k++) 
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
  for(int hi = 0; hi < mdl.getNumPartials(); hi++)
    for(int di = 0; di < getNumDataPoints(); di++)
      mdl.getDataRef(hi, di).time /= sampleRate;
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
void rsHarmonicAnalyzer<T>::setCycleLength(int newLength)
{
  cycleLength = newLength;
  blockSize   = cyclesPerBlock * cycleLength;
  sig.resize(blockSize); 
  wnd.resize(blockSize);
  fillWindow();
  trafoSize = zeroPad * blockSize;
  sigPadded.resize(trafoSize);
  mag.resize(trafoSize); 
  phs.resize(trafoSize);
  trafo.setBlockSize(trafoSize);
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
  T peakSearchWidth = T(1); // maybe make user parameter later
  T mainlobeWidth   = rsWindowFunction::getMainLobeWidth(windowType, T(0));
  return (int) round(T(0.5)*peakSearchWidth*zeroPad*mainlobeWidth);
}

template<class T>
void rsHarmonicAnalyzer<T>::fillHarmonicData(
  RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T time)
{
  //rsPlotVector(sig);
  prepareBuffer(sig, sigPadded);
  trafo.getRealSignalMagnitudesAndPhases(&sigPadded[0], &mag[0], &phs[0]);  // perform FFT


  if(frameIndex == getNumFrames()/2)
    rsPlotSpectrum(mag, T(0), T(-100), true); // freq axis wrong, if we pass the sampleRate
    //rsPlotVector(mag);

  //if(frameIndex >= 10) {
  //  rsPlotSpectrum(mag, sampleRate, T(-200));
  //  //rsPlotVector(rsAmpToDb(mag, T(-200))); 
  //}


  // extract model data from FFT result:
  int dataIndex   = frameIndex + 1;        // +1 because of the fade-in datapoint
  int numBins     = trafo.getBlockSize();  // number of FFT bins
  int numPartials = getNumHarmonics();     // number of (pseudo) harmonics - new
  int k;                                   // bin index
  T freq, gain, phase;

  int w2 = getSpectralPeakSearchWidth(); // to be used later

  if(zeroPad == 1 && cyclesPerBlock == 1) { // old version before multi-cycle and zero-padding
    for(k = 0; k < numPartials; k++) {
      freq = trafo.binIndexToFrequency(k, numBins, sampleRate);
      mdl.setData(k, dataIndex, time, freq, T(2)*mag[k], phs[k]);
    }
  }
  else
  {
    // we will need parabolic interpolation to find better frequency (and amplitude) measurements 
    // and for interpolating phase-data, we need to be careful about wrapping issues
    mdl.setData(0, dataIndex, time, T(0), T(2*zeroPad)*mag[0], phs[0]); // handle DC separately
    for(int h = 1; h < numPartials; h++) {

      k = cyclesPerBlock*zeroPad*h; // bin-index where we expect the peak for the (pseudo) harmonic

      bool parabolicInterpolation = true;  // make user option
      // with the tremolo-sine, not doing this is actually better, - why? maybe the search range 
      // for a maximum is too large and we pick up sidelobes of the window? ...maybe we should look
      // for a maximum *or* a minimum?

      if(parabolicInterpolation)
      {
        // not yet finished

        //int kPeak = findPeakBinNear(mag, k, zeroPad);   // old

        //int kPeak = findPeakBinNear(mag, k, zeroPad/cyclesPerBlock); 
        // best for tremolo-sine, but then it misses the inharmonic for the 200Hz/6100Hz example
        // ...we may need some more sophisticated way to find the relevant bin-index

        int kPeak = findPeakBinNearOld(mag, k, zeroPad*cyclesPerBlock); 
        // worst for tremolo-sine, best for two sines 200Hz/6100Hz
        // why should the cyclesPerBlock make the search range wider - if anything, it should make 
        // it narrower - but with the two versions above, we miss the 6100 partial
        // -> try to figure out most meaningful search range and/or let the user select it

        if(kPeak == -1) {
          // no peak found: use k for the frequency, zero for the amplitude
          freq = trafo.binIndexToFrequency(k, numBins, sampleRate);
          gain = T(0);
          phase = phs[k];
        } else {
          // preliminary:
          freq  = trafo.binIndexToFrequency(kPeak, numBins, sampleRate);
          gain  = T(2*zeroPad)*mag[kPeak]; // preliminary - compute parabola maximum
          phase = phs[kPeak];
          // todo: find exact frequency and amplitude by parabolic interpolation:
        }
        mdl.setData(h, dataIndex, time, freq, gain, phase);
      }
      else {
        freq = trafo.binIndexToFrequency(k, numBins, sampleRate);
        gain = T(2*zeroPad)*mag[k]; 
        mdl.setData(h, dataIndex, time, freq, gain, phs[k]);
      }
    }
  }


  int dummy = 0;
}

// move to rsArray:
template<class T>
int numPeaks(T* x, int N)
{
  int np = 0;
  for(int i = 1; i < N-1; i++)
    if(x[i] >= x[i-1] && x[i] >= x[i+1])
      np++;
  return np;
}

template<class T>
int rsHarmonicAnalyzer<T>::findPeakBinNear(std::vector<T>& v, int kCenter, int w2)
{
  if(w2 == 0)
    return kCenter;
  int kLeft  = rsMax(kCenter - w2, 0);
  int kRight = rsMin(kCenter + w2, (int) v.size()-1);
  int length = kRight - kLeft + 1;
  //rsPlotArray(&v[kLeft], length); // plot segement where we search for a peak
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
  int nPeaks = numPeaks(&v[kLeft], length);
  if(nPeaks != 1)
    return -1;
  // maybe make enforcing this unimodality condition optional

  // find index of maximum:
  int kMax = rsArray::maxIndex(&v[kLeft], length) + kLeft; // check this

  if(kMax == kLeft || kMax == kRight) // isn't this ensured already by nPeaks == 1?
    return -1;
  else
    return kMax;
}

// old version - to be reomoved:
template<class T>
int rsHarmonicAnalyzer<T>::findPeakBinNearOld(std::vector<T>& v, int k, int w)
{
  if(w == 1)
    return k;
  rsAssert(rsIsEven(w));
  int w2 = w/2;
  rsAssert(k-w2 >= 0 && k+w2+1 < (int) v.size());


  //rsPlotArray(&v[k-w2], w+1); // plot segement where we search for a peak


  int nPeaks = numPeaks(&v[k-w2], w+1);
  //if(nPeaks != 1)
  //  return -1;  
  // there's no single peak - but a single peak is the condition that we have a partial in this
  // segment of the spectrum
  // this sometimes misses the peak at 6100 - why ...but only in the later frames...
  // maybe using a blackman window could help because its mainlobe is wider

  //int kMax = rsArray::maxIndex(&v[k-w2], w) + k - w2;  // integer index, old
  int kMax = rsArray::maxIndex(&v[k-w2], w+1) + k - w2;  // integer index, new
  //if(kMax == k-w2 || kMax == k+w2) {  // old
  if(kMax == k-w2 || kMax == k+w2+1) {  // new
    // there is no actual peak - the max-value is at the boundary - maybe we should do something
    // special in this case? ...maybe just return the center k?
    return -1;
    // note that in this case, there may be actually a minimum at k, so the parabolic amplitude 
    // computation may end up with a negative value...maybe clip amplitudes at zero from below
  } else {
    //rsPlotArray(&v[k-w2], w+1); // plot segement where we search for a peak
    return kMax; // preliminary - todo: parabolic interpolation (of log-values)
  }

  // -smaller search ranges are better for the tremolo-sine, because with a larger range, it may
  //  pick up the sidelobes of the window
  // -larger ranges better for the inharmonic pair (200/6100) because a too small range will miss
  //  the inharmonic partial
  // -can we find an algorithm that works well in both cases? -> make some plots


  // maybe, we should also check, if the found peak is higher than its two neighbours and maybe 
  // if we don't find an actual peak, we should return -1 to indicate that - maybe the outlying 
  // code can the decide to set the amplitude to zero

  // maybe we should also check, if the maximum is some threshold above the min or mean value
  // and/or whether the curve is unimodal, i.e. has exactly one local peak

  // this unimodality condition seems to make a lot of sense when we use zeroPad*cyclesPerBlock as 
  // search length because the mainlobe is wider than the sidelobes - sidelobes can be detected by
  // having a multiple local peaks
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

  // todo: in the first and last loop, multiply by window during copying, i.e.
  // sig[i+K2] becomes sig[i+K2] * wnd[i+K2] and sig[i] becomes sig[i] * wnd[i]
}

template<class T>
void rsHarmonicAnalyzer<T>::fillWindow()
{
  rsWindowFunction::createWindow(&wnd[0], (int) wnd.size(), windowType, true);
  //rsWindowFunction::createWindow(&wnd[0], (int) wnd.size(), windowType, false);
}


//template class rsHarmonicAnalyzer<double>;