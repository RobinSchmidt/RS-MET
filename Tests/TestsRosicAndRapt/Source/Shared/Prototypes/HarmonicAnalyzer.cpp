
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


template<class T>
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N, T fs)
{
  RAPT::rsSinusoidalModel<T> mdl;


  // pre-processing (flatten pitch):
  //RAPT::rsPitchFlattener<T, T> flattener;
  //...

  typedef std::vector<T> Vec;
  Vec cycleMarks = findCycleMarks(x, N, fs);    // cycle marks
  if(cycleMarks.size() < 2)
    return mdl;
  Vec cycleLengths = rsDifference(cycleMarks);  // cycle lengths
  T maxLength = rsMax(cycleLengths);
  int targetLength = RAPT::rsNextPowerOfTwo((int) ceil(maxLength));



  // todo: 
  // -create a time-warping map that maps the measured markers to their target instants
  //  -the very first marker should not be moved...how should we handle the first partial cycle?
  //   ...zero padding? extrapolation of freq and phase and start at zero amplitude? maybe let
  //   the user decide?
  //   -however this will be handled, we should probably not move the first marker..or wait
  //    ..we should probably scale its position by targetLength / cl[0]

  // create the mapping function for the time instants
  int mapLength = (int) cycleMarks.size();  // maybe we need +1 to map the N-1 th input-sample?
  Vec tIn(mapLength), tOut(mapLength);
  tIn[0]  = tOut[0] = 0;
  tIn[1]  = cycleMarks[0];
  tOut[1] = cycleMarks[0] * targetLength / cycleLengths[0];
  for(int i = 2; i < mapLength; i++) {
    tIn[i]  = cycleMarks[i-1];
    tOut[i] = tOut[i-1] + targetLength;
  }
  // verify this...actually, we want the tOut values to land on integers 2^k + offset where the 
  // offset comes from the initial partial cycle

  int dummy = 0;


  // the distance of the very first marker from the time origin t=0 should probably used for 
  // determining the start phase - don't assume an additional "ghost" marker at t=0 - instead, let
  // the sinusoid start at zero amplitude, frequency determined by the distance between 1st and 2nd
  // marker and phase appropriate to the frequency and time-value of the 1st marker (i.e. if the 
  // first marker is at 25 and the second is at 125, assume a cycle length of 100 and start phase
  // of -90° (a quarter period) - for higher harmonics, take into account the phase-measurement
  // at first marker (for the fundamental, that phase is zero by construction)

  // harmonic data extraction:

  RAPT::rsFourierTransformerRadix2<T> trafo;


  // post-processing (account for pitch flattenig):

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


template class rsHarmonicAnalyzer<double>;