

// move elsewhere...
template<class T>
void plotSignalWithMarkers(T* signal, int signalLength, T* markers, int numMarkers)
{
  std::vector<T> zeros(numMarkers);    // y values for plotting (all zero)
  RAPT::rsArray::fillWithZeros(&zeros[0], numMarkers);
  GNUPlotter plt;
  plt.addDataArrays(signalLength, signal);
  plt.addDataArrays(numMarkers,   markers, &zeros[0]);
  plt.setGraphStyles("lines", "points");
  plt.setPixelSize(1000, 300);
  plt.plot();
}

// move to rapt:
template<class T>
std::vector<T> rsDifference(const std::vector<T> x)
{
  if(x.size() < 2)
    return std::vector<T>();  // result is empty
  std::vector<T> d(x.size()-1);
  for(size_t i = 0; i < d.size(); i++)
    d[i] = x[i+1] - x[i];
  return d;
}

//template<class T>
//T rsMinValue(T 

template<class T>
T rsMax(const std::vector<T> x)
{
  T max = std::numeric_limits<T>::min(); 
  // we should instead use -inf for double/float? -> make explicit specilizations

  for(size_t i = 0; i < x.size(); i++) {
    if(x[i] > max)
      max = x[i];
  }
  return max;
}


//-------------------------------------------------------------------------------------------------

template<class T>
rsHarmonicAnalyzer<T>::rsHarmonicAnalyzer()
{
  typedef rsFourierTransformerRadix2<T> FT;
  trafo.setDirection(FT::FORWARD);
  trafo.setNormalizationMode(FT::NORMALIZE_ON_FORWARD_TRAFO);
}

template<class T>
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N)
{
  RAPT::rsSinusoidalModel<T> mdl; // init empty model
  if(preProcess(x, N) == false)   // pre-process audio (flatten pitch)
    return mdl;                   // return empty model if pre-processing has failed
  analyzeHarmonics(mdl);          // create model from pitch-flattened signal (now in member y)
  postProcess(mdl);               // post process model data to account for flattening
  return mdl;

  //rosic::writeToMonoWaveFile("ModalPluckStretched.wav", &y[0], (int)y.size(), (int)sampleRate);

  // todo: clean up the model: remove partials that are consistently above the nyquist freq - due 
  // to the stretching, we may get frequencies almost up to the sample-rate (we downshift at most 
  // by an octave (with respect to the longest original cycle) - so the model my end up having 
  // twice the number of frequencies as it should have...(although their amplitudes are really low)
  // ...maybe also introduce an amplitude threshold and remove partials that are consistently below
  // that threshold
}

template<class T>
bool rsHarmonicAnalyzer<T>::preProcess(T* x, int Nx)
{
  typedef RAPT::rsArray AR;
  typedef std::vector<T> Vec;

  // Find cycle marks and assign FFT blockSize:
  Vec cycleMarks = findCycleMarks(x, Nx);        // cycle marks
  if(cycleMarks.size() < 2)
    return false;                                // report failure
  Vec cycleLengths = rsDifference(cycleMarks);   // cycle lengths
  T maxLength = rsMax(cycleLengths);
  blockSize = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));

  // Create the mapping function for the time instants of the cycle marks:
  int mapLength = (int) cycleMarks.size() + 2;  // +2 for t = 0 and t = N-1
  tIn.resize( mapLength);
  tOut.resize(mapLength);
  tIn[0] = tOut[0] = 0;    // time-origin is zero for both, original and stretched signal

  // The first marker is mapped to an instant, such that the initial partial cycle is stretched by
  // the same amount as the first full cycle (the cycle between 1st and 2nd marker):
  tIn[1]  = cycleMarks[0];
  tOut[1] = cycleMarks[0] * blockSize / cycleLengths[0];
  tOut[1] = round(tOut[1]);

  // All cycles between the initial partial cycle and final partial cycle are stretched to the same 
  // fixed length:
  for(int i = 2; i < mapLength-1; i++) {
    tIn[i]  = cycleMarks[i-1];
    tOut[i] = tOut[i-1] + blockSize;
  }

  // The end time instant is mapped such that the final partial cycle is stretched by the same 
  // amount as the last full cycle:
  T tailLength = (Nx-1) - rsLast(cycleMarks);
  tIn [mapLength-1] = Nx-1;
  tOut[mapLength-1] = tOut[mapLength-2] + tailLength * blockSize / rsLast(cycleLengths);
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
  typedef RAPT::rsArray AR;

  // Initialize the model (create all datapoints, to filled with actual data later):
  int numHarmonics  = blockSize / 2;        // number of partials
  int numFrames     = getMapLength() - 1;   // number of datapoints in each partial
  int numDataPoints = numFrames;            // maybe later use numFrames+2 for fade-in/out
  mdl.init(numHarmonics, numDataPoints);

  // Set up the fourier transformer object and block buffers:
  int K = blockSize;       // block-size (equals FFT size)
  trafo.setBlockSize(K);
  sig.resize(K); 
  mag.resize(K); 
  phs.resize(K);

  // The initial partial cycle is pre-padded with zeros:
  int n0 = 0;                          // first sample (from y-array) in current frame
  int m  = 0;                          // frame index
  int L  = (int) tOut[1];              // length of initial partial cycle
  AR::fillWithZeros(&sig[0], K-L);
  AR::copyBuffer(&y[n0], &sig[K-L], L);
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));

  // The inner cycles/frames are taken as is:
  L = K;                               // length of inner cycles
  for(m = 1; m < numFrames-1; m++) {
    n0 = (int) tOut[m];
    AR::copyBuffer(&y[n0], &sig[0], L);
    fillHarmonicData(mdl, m, getTimeStampForFrame(m));
  }

  // The final partial cycle is post-padded with zeros:
  n0 = (int) tOut[m]; 
  L = int(tOut[tOut.size()-1] - tOut[tOut.size()-2]);  // maybe +1? check against off-by-1
  AR::copyBuffer(&y[n0], &sig[0], L);
  AR::fillWithZeros(&sig[L], K-L);
  fillHarmonicData(mdl, m, getTimeStampForFrame(m));

  // todo: double-check all index computations against off-by-one errors, verify time-indices
  // maybe make some plots

  // the distance of the very first marker from the time origin t=0 should probably used for 
  // determining the start phase - don't assume an additional "ghost" marker at t=0 - instead, let
  // the sinusoid start at zero amplitude, frequency determined by the distance between 1st and 2nd
  // marker and phase appropriate to the frequency and time-value of the 1st marker (i.e. if the 
  // first marker is at 25 and the second is at 125, assume a cycle length of 100 and start phase
  // of -90° (a quarter period) - for higher harmonics, take into account the phase-measurement
  // at first marker (for the fundamental, that phase is zero by construction)
}

template<class T>
void rsHarmonicAnalyzer<T>::postProcess(RAPT::rsSinusoidalModel<T>& mdl)
{
  int numHarmonics  = blockSize / 2;        // number of partials
  int numFrames     = getMapLength() - 1;   // number of datapoints in each partial
  int numDataPoints = numFrames;            // maybe later use numFrames+2 for fade-in/out

  // todo: get rid of these vectors - compute values on the fly inside the loop
  std::vector<T> lw = rsDifference(tOut);  // warped lengths of cycles
  std::vector<T> lu = rsDifference(tIn);   // unwarped lengths of cycles
  for(int m = 0; m < numFrames; m++) {
    T tw = getTimeStampForFrame(m);         // warped time
    T tu = getUnWarpedTimeStampForFrame(m); // unwarped time
    T r  = lw[m] / lu[m];                   // stretching ratio applied to frame m
    for(int k = 0; k < numHarmonics; k++) {
      mdl.getDataRef(k, m).time  = tu;
      mdl.getDataRef(k, m).freq *= r;
    }
    int dummy = 0;
  }

  // convert time data from samples to seconds (multiply all time-stamps by 1/sampleRate):
  for(int hi = 0; hi < numHarmonics; hi++)
    for(int di = 0; di < numDataPoints; di++)
      mdl.getDataRef(hi, di).time /= sampleRate;
}


template<class T>
std::vector<T> rsHarmonicAnalyzer<T>::findCycleMarks(T* x, int N)
{
  T fl = 20;       // lower limit for fundamental (maybe let user set this up)
  T fu = 5000;     // upper limit for fundamental
  rsCycleMarkFinder<double> cmf(sampleRate, fl, fu);
  cmf.setSubSampleApproximationPrecision(2);  // 0: linear, 1: cubic, 2: quintic, ...
  cmf.setAlgorithm(rsCycleMarkFinder<double>::F0_ZERO_CROSSINGS);
  std::vector<T> cm = cmf.findCycleMarks(x, N);
  //plotSignalWithMarkers(x, N, &cm[0], (int) cm.size());
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
void rsHarmonicAnalyzer<T>::fillHarmonicData(
  RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T time)
{
  int dataIndex = frameIndex;
  // maybe we should use numDataPoints = numFrames+2 for prepending and appending a datapoint with
  // zero amplitude for fade-in/out (for each partial) - then we should init the model with
  // numDataPoints instead of numFrames and use dataIndex = frameIndex+1 here - the dataPoints
  // at indices 0 and numDataPoints-1 have to be treated separately..

  int K = trafo.getBlockSize();
  int numPartials = K/2; // maybe -1 for not including DC...or maybe just include
                         // DC into the model


  // we should shift the signal buffer by one half to move the time origin to the center
  // ...yes - at least, if we use the center of the frame for the time-stamp we need to do this
  // because only then it is all consistent...maybe we should have a function 
  // prepareSignalBlockBuffer(int frameIndex)

  // this is realy bad to do this here (it has an undocumented side-effect on the sig buffer - it 
  // doesn't matter, because the buffer is not used anymore after calling this function but still).
  // -> refactor and also optimize by passing a temp-buffer to the function:
  RAPT::rsArray::circularShift(&sig[0], K, -K/2);
  // optimize: pass tmpBuffer as workspace to avoid temporary memory allocation

  // for the first and last (partial) cycle that may not be correct...or is it? hmm...yes,
  // it could be
  // verify, if we should really shift by -K/2 (or if this ends up off-by-one)

 
  trafo.getRealSignalMagnitudesAndPhases(&sig[0], &mag[0], &phs[0]);       // perform FFT


  // extract model data from FFT result:
  for(int k = 0; k < numPartials; k++) {
    T freq = trafo.binIndexToFrequency(k, K, sampleRate);
    mdl.setData(k, dataIndex, time, freq, T(2)*mag[k], phs[k]);
  }
  int dummy = 0;
}

template class rsHarmonicAnalyzer<double>;