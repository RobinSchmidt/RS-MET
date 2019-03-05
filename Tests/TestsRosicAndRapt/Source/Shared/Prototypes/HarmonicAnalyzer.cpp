
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


template<class T>
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N, T fs)
{
  // pre-processing (flatten pitch):
  //RAPT::rsPitchFlattener<T, T> flattener;
  //...

  typedef std::vector<T> Vec;
  Vec cm = findCycleMarks(x, N, fs);

  // todo: 
  // -find greatest distance dMax between adjacent cycle marks
  // -compute targetDistance dTgt = nextPowerOfTwo(dMax)
  // -create a time-warping map that maps the mesured markers to theri target instants


  // the distance of the very first marker from the time origin t=0 should probably used for 
  // determining the start phase - don't assume an additional "ghost" marker at t=0 - instead, let
  // the sinusoid start at zero amplitude, frequency determined by the distance between 1st and 2nd
  // marker and phase appropriate to the frequency and time-value of the 1st marker (i.e. if the 
  // first marker is at 25 and the second is at 125, assume a cycle length of 100 and start phase
  // of -90° (a quarter period) - for higher harmonics, take into account the phase-measurement
  // at first marker (for the fundamental, that phase is zero by construction)

  // harmonic data extraction:
  RAPT::rsSinusoidalModel<T> mdl;
  RAPT::rsFourierTransformerRadix2<T> trafo;


  // post-processing (account for pitch flattenig):

  return mdl;
}

template<class T>
std::vector<T> rsHarmonicAnalyzer<T>::findCycleMarks(T* x, int N, T fs)
{
  T fl = 20;       // lower limit for fundamental
  T fu = 5000;     // upper limit for fundamental
  rsCycleMarkFinder<double> cmf(fs, fl, fu);
  cmf.setSubSampleApproximationPrecision(2);  // 0: linear, 1: cubic, 2: quintic,  ...
  cmf.setAlgorithm(rsCycleMarkFinder<double>::F0_ZERO_CROSSINGS);
  std::vector<T> cm = cmf.findCycleMarks(x, N);
  //plotSignalWithMarkers(x, N, &cm[0], (int) cm.size());
  return cm;
}


template class rsHarmonicAnalyzer<double>;