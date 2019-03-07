

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
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N, T fs)
{
  typedef RAPT::rsArray AR;

  RAPT::rsSinusoidalModel<T> mdl;


  //-----------------------------------------------------------------------------------------------
  // Step 1: - pre-processing (flatten pitch):
  // todo: factor out...maybe the function preProcess() should return a bool to indicate success

  typedef std::vector<T> Vec;
  Vec cycleMarks = findCycleMarks(x, N, fs);    // cycle marks
  if(cycleMarks.size() < 2)
    return mdl;
  Vec cycleLengths = rsDifference(cycleMarks);  // cycle lengths
  T maxLength = rsMax(cycleLengths);
  int targetLength = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));

  // create the mapping function for the time instants
  int mapLength = (int) cycleMarks.size() + 2;  // +2 for t = 0 and t = N-1

  Vec tIn(mapLength), tOut(mapLength);  
  // maybe these should be members (they are needed for the post-processing step, too - and we want
  // to factor out functions for pre- and post processing)

  // the time-origin is zero for both, original and stretched signal:
  tIn[0]  = tOut[0] = 0;

  // the first marker is mapped to an instant, such that the initial partial cycle is stretched by
  // the same amount as the first full cycle (the cycle between 1st and 2nd marker):
  tIn[1]  = cycleMarks[0];
  tOut[1] = cycleMarks[0] * targetLength / cycleLengths[0];
  tOut[1] = round(tOut[1]);

  // all cycles between the initial partial cycle and final partial cycle are stretched to the same 
  // fixed length
  for(int i = 2; i < mapLength-1; i++) {
    tIn[i]  = cycleMarks[i-1];
    tOut[i] = tOut[i-1] + targetLength;
  }

  // the end time instant is mapped such that the final partial cycle is stretched by the same 
  // amount as the last full cycle:
  double tailLength = (N-1) - rsLast(cycleMarks);
  tIn [mapLength-1] = N-1;
  tOut[mapLength-1] = tOut[mapLength-2] + tailLength * targetLength / rsLast(cycleLengths);
  tOut[mapLength-1] = round(tOut[mapLength-1]);

  Vec test; // for debug
  test = rsDifference(tOut);
  // elements should be all equal to targetLength except the first and the last (which should be
  // shorter than that) - ok - looks good

  // ok, we have created the time warping map, sampled at the cycle-marks, for applying the 
  // warping, we need to interpolate it up to sample rate - we use linear interpolation for that:
  int Ny = (int) rsLast(tOut) + 1; // length of stretched signal and warping map
  Vec t(Ny), w(Ny);                // interpolated time axis and warping map
  AR::fillWithIndex(&t[0], Ny);
  RAPT::resampleNonUniformLinear(&tOut[0], &tIn[0], mapLength, &t[0], &w[0], Ny);
  //test = rsDifference(w); // should be the readout-speed

  // do the time-warping:
  double sincLength = 64.0;        // length of sinc-interpolator
  Vec y(Ny);                       // stretched signal
  rsTimeWarper<T, T>::timeWarpSinc(x, N, &y[0], &w[0], Ny, sincLength);

  // pre-processing done: y contains the pre-processed (pitch-flattened) signal

  //rosic::writeToMonoWaveFile("StretchedModalPluck.wav", &y[0], Ny, (int)fs);

  //-----------------------------------------------------------------------------------------------
  // Step 2: - analyze harmonics in flattened signal and write data into sinusoidal model:

  // initialize the model (create all datapoints, to filled with actual data later):
  int numHarmonics = targetLength / 2;     // number of partials
  int numFrames    = mapLength;            // number of datapoints in each partial
  mdl.init(numHarmonics, numFrames);

  // set up the fourier transformer object and block buffers:
  int M = targetLength;    // block-size (equals FFT size) - maybe use K, and M for numFrames
  trafo.setBlockSize(M);
  sig.resize(M); 
  mag.resize(M); 
  phs.resize(M);

  // the initial partial cycle is pre-padded with zeros:
  //int n;
  int L = (int) tOut[1];                 // length of initial partial cycle
  AR::fillWithZeros(&sig[0], M-L);
  AR::copyBuffer(&y[0], &sig[M-L], L);
  fillHarmonicData(mdl, 0, T(0.5)*(tOut[1]-tOut[0]));
    // maybe make a function getTimeStamp(int i) - it should return T(0.5)*(tOut[i+1]-tOut[i]) 
    // or maybe just tOut[i], depending on user settings. for this, we need to make the tOut array
    // a member (which we need anyway for factoring out pre- and post processing steps)

  // the inner cycles/frames are taken as is:
  L = M;
  int m;  // frame index
  for(m = 1; m < numFrames-1; m++) {
    int frameStart = (int) tOut[m];
    AR::copyBuffer(&y[frameStart], &sig[0], L);
    fillHarmonicData(mdl, m, T(0.5)*(tOut[m]-tOut[m+1]));
  }


  int dummy = 0;


  // the final partial cycle is post-padded with zeros:
  // L = ...
  // ...







  //-----------------------------------------------------------------------------------------------
  // Step 3: - post process data in model to account for the flattening:




  // todo: maybe factor out the pre- and post-processing into a class:
  // rsFlatPitchPrePostProcessor
  //   std::vector<T> preProcessAudio(T* x, int N);
  //   postProcessModel(rsSinusoidalModel<T>& model);
  // in between these two calls, we do the harmonic extraction, the class keeps around the warping 
  // map after pre-processing the audio and uses it again in the post-processing step







  // the distance of the very first marker from the time origin t=0 should probably used for 
  // determining the start phase - don't assume an additional "ghost" marker at t=0 - instead, let
  // the sinusoid start at zero amplitude, frequency determined by the distance between 1st and 2nd
  // marker and phase appropriate to the frequency and time-value of the 1st marker (i.e. if the 
  // first marker is at 25 and the second is at 125, assume a cycle length of 100 and start phase
  // of -90° (a quarter period) - for higher harmonics, take into account the phase-measurement
  // at first marker (for the fundamental, that phase is zero by construction)






  return mdl;
}

template<class T>
std::vector<T> rsHarmonicAnalyzer<T>::findCycleMarks(T* x, int N, T fs)
{
  T fl = 20;       // lower limit for fundamental (maybe let user set this up)
  T fu = 5000;     // upper limit for fundamental
  rsCycleMarkFinder<double> cmf(fs, fl, fu);
  cmf.setSubSampleApproximationPrecision(2);  // 0: linear, 1: cubic, 2: quintic, ...
  cmf.setAlgorithm(rsCycleMarkFinder<double>::F0_ZERO_CROSSINGS);
  std::vector<T> cm = cmf.findCycleMarks(x, N);
  //plotSignalWithMarkers(x, N, &cm[0], (int) cm.size());
  return cm;
}

template<class T>
void rsHarmonicAnalyzer<T>::fillHarmonicData(
  RAPT::rsSinusoidalModel<T>& mdl, int frameIndex, T timeStamp)
{
  trafo.getRealSignalMagnitudesAndPhases(&sig[0], &mag[0], &phs[0]);
  int K = trafo.getBlockSize();
  int numPartials = K/2; // maybe -1 for not including DC...or maybe just include
                         // DC into the model

  for(int k = 0; k < numPartials; k++)
  {
    //T freq = trafo.binIndexToFrequency(k, K, sampleRate);

  }



  int dummy = 0;
}

template class rsHarmonicAnalyzer<double>;