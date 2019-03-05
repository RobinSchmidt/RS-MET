
// move elsewhere...
template<class T>
void plotSignalWithMarkers(T* signal, int signalLength, T* markers, int numMarkers)
{
  GNUPlotter plt;
  plt.addDataArrays(signalLength, signal);

  // add markers...

  plt.plot();

  //plotVector(cm);
}


template<class T>
RAPT::rsSinusoidalModel<T> rsHarmonicAnalyzer<T>::analyze(T* x, int N, T fs)
{
  // pre-processing (flatten pitch):
  //RAPT::rsPitchFlattener<T, T> flattener;
  //...

  typedef std::vector<T> Vec;
  T fl = 20;       // lower limit for fundamental
  T fu = 5000;     // upper limit for fundamental


  //// from removePitchModulation:
  //Vec fi(N);       // instantaneous frequency at each sample
  //rsInstantaneousFundamentalEstimatorD::measureInstantaneousFundamental(x, fi, N, fs, fl, fu);
  //// Compute desired readout speed for each sample which is given by the ratio of the desired 
  //// instantaneous frequency and the actual instantaneous frequency of the input signal (re-use the
  //// f-array for this)
  //for(int n = 0; n < N; n++)
  //  f[n] = targetFrequency / f[n];

  // find the cycle marks:
  rsCycleMarkFinder<double> cmf(fs, fl, fu);
  //cmf.setFundamentalRange(fl, fu);
  //cmf.setSampleRate(fs);
  cmf.setSubSampleApproximationPrecision(2);  // 0: linear, 1: cubic, ...
  cmf.setAlgorithm(rsCycleMarkFinder<double>::F0_ZERO_CROSSINGS);
  Vec cm = cmf.findCycleMarks(x, N);

  plotSignalWithMarkers(x, N, &cm[0], (int) cm.size());








  // harmonic data extraction:
  RAPT::rsSinusoidalModel<T> mdl;
  RAPT::rsFourierTransformerRadix2<T> trafo;


  // post-processing (account for pitch flattenig):




  return mdl;
}

template class rsHarmonicAnalyzer<double>;