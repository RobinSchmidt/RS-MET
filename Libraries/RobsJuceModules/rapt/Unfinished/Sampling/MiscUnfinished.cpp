//=================================================================================================

// x: input, y: output, N: length, processor: the filter/processor to be applied (needs to
// implement double getSample(double), P: padding length, numPasses: number of forward/backward
// passes maybe move to another file
template<class TSig, class TFlt> // signal and filter type
void rsApplyBiDirectionally(
  const TSig *x, TSig *y, int N, TFlt &processor, int P, int numPasses)
{
  // create a buffer containing the signal with (pre- and post) zero padding to allow the
  // processor/filter to ring out at the ends:
  int M = N+2*P;
  TSig *tmp = new TSig[M];
  rsArray::fillWithZeros(tmp, P);
  rsArray::copy(x, &tmp[P], N);
  rsArray::fillWithZeros(&tmp[P+N], P);

  // apply processor (multipass, bidirectionally):
  int p, n;
  for(p = 1; p <= numPasses; p++) {
    for(n = 0;   n <  M; n++) tmp[n] = processor.getSample(tmp[n]); // forward pass
    for(n = M-1; n >= 0; n--) tmp[n] = processor.getSample(tmp[n]); // backward pass
  }

  // copy result to output and clean up:
  rsArray::copy(&tmp[P], y, N);
  delete[] tmp;

  // todo: Actually, we don't really need pre- and post padding. Using just post padding and
  // applying first all forward passes and then all backward passes should give the same result.
  // ...but maybe it would show a different numeric roundoff behavior? -> test it - maybe make the
  // function a (static) member of rsBiDirectionalFilter
  // maybe it's better in practice to have all of the ring-out/warm-up section at the end of the
  // sample because samples may have a ahrd transient at the start but are already decayed down
  // toward the end such that at the end, the total amount of ringing effects is even less
}

template<class T>
std::vector<T> getPaddedSignal(const T* x, int N, int P)
{
  int M = N+2*P;
  std::vector<T> xp(M);
  rsArray::fillWithZeros(&xp[0], P);    // pre-padding
  rsArray::copy(x, &xp[P], N);    // actual signal
  rsArray::fillWithZeros(&xp[P+N], P);  // post-padding
  return xp;
}
// todo: use this function also for the function above - gets rid of code duplication

// non-uniform version of function above:
template<class TSig, class TTim, class TFlt> // signal, time and filter type
void rsApplyBiDirectionally(
  const TSig* x, const TTim* t, TSig* y, int N, TFlt& processor, int P, int numPasses)
{
  // signal array with pre- and post-padding:
  std::vector<TSig> tmp = getPaddedSignal(x, N, P);
  int M = (int) tmp.size(); // M = N+2*P

  // pre/post-padded time-delta array:
  std::vector<TTim> dt = getPaddedSignal(t, N, P);

  //size_t i = 0;
  //dt[0] = TTim(1); // test - gcc/win produces "assignment of read-only location" error - also in
             // the loops below - WTF? why is a vector element "read-only"?

  TTim dtAv = (t[N-1]-t[0]) / N; // average sampling interval
  int n;
  for(n = N+P-1; n >      P; n--) dt[n] = dt[n] - dt[n-1];
  for(n = 0;     n <=     P; n++) dt[n] = dtAv;   // don't do this before the 1st loop!
  for(n = N+P;   n <  N+2*P; n++) dt[n] = dtAv;

  // apply processor (multipass, bidirectionally):
  for(int p = 1; p <= numPasses; p++) {
    for(n = 0;   n <  M; n++) tmp[n] = processor.getSample(tmp[n], dt[n]);   // forward pass
    for(n = M-2; n >= 0; n--) tmp[n] = processor.getSample(tmp[n], dt[n+1]); // backward pass
  }

  // copy result to output:
  rsArray::copy(&tmp[P], y, N);
}

template<class T>
int rsBiDirectionalFilter::getPaddingLength(T bw, T fs)
{
  return rsCeilInt(10 * fs / bw);
  // factor 10 is ad hoc - experiment to find optimal factor, maybe the formula should include
  // the number of passes as well? and maybe also the order? ...ideally, it should be based on the
  // filter's ringing time
}

template<class TSig, class TPar>
void rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(const TSig *x, TSig *y, int N, TPar fc,
  TPar bw, TPar fs, int numPasses, TPar gc)
{
  // compute desired bandwidth for single-pass filter in octaves:
  bw *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, 1, gc);
  TPar bo = rsBandwidthConverter::absoluteBandwidthToOctaves(bw, fc);

  // create and set up the filter:
  rsStateVariableFilter<TSig, TPar> flt;  // maybe use a biquad later
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setBandwidth(bo);
  flt.setMode(rsStateVariableFilter<TSig, TPar>::BANDPASS_PEAK);

  // apply filter:
  int P = getPaddingLength(bw, fs);
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

template<class TSig, class TPar>
void rsBiDirectionalFilter::applyButterworthBandpassBwInHz(const TSig *x, TSig *y, int N, TPar fc,
  TPar bw, TPar fs, int order, int numPasses, TPar gc)
{
  // compute desired bandwidth for single-pass filter in octaves:
  bw *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);
  TPar bo = rsBandwidthConverter::absoluteBandwidthToOctaves(bw, fc);

  // create and set up the filter:
  rsEngineersFilter<TSig, TPar> flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setBandwidth(bo);
  flt.setMode(rsInfiniteImpulseResponseDesigner<TPar>::BANDPASS);
  flt.setApproximationMethod(rsPrototypeDesigner<TPar>::BUTTERWORTH);
  flt.setPrototypeOrder(order);

  // apply filter:
  int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

template<class TSig, class TPar>
void rsBiDirectionalFilter::applyButterworthLowpass(const TSig *x, TSig *y, int N, TPar fc,
  TPar fs, int order, int numPasses, TPar gc)
{
  // compute desired lowpass cutoff:
  fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);

  // create and set up the filter:
  rsEngineersFilter<TSig, TPar> flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setMode(rsInfiniteImpulseResponseDesigner<TPar>::LOWPASS);
  flt.setApproximationMethod(rsPrototypeDesigner<TPar>::BUTTERWORTH);
  flt.setPrototypeOrder(order);

  // apply filter:
  int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

template<class TSig, class TPar>
void rsBiDirectionalFilter::applyButterworthHighpass(const TSig *x, TSig *y, int N, TPar fc,
  TPar fs, int order, int numPasses, TPar gc)
{
  applyButterworthLowpass(x, y, N, fc, fs, order, numPasses, gc);
  rsArray::subtract(x, y, y, N); // works because we use a bidirectional (zero-phase) filter
  // but is this really equivalent to a proper bidirectional Butterworth highpass? i think so, but
  // this should to be verified....


  //// quick and dirty code duplication from applyButterworthLowpass - refactor

  //// compute desired lowpass cutoff:
  //fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);

  //// create and set up the filter:
  //rsEngineersFilter flt;
  //flt.setSampleRate(fs);
  //flt.setFrequency(fc);
  //flt.setMode(rsInfiniteImpulseResponseDesigner::HIGHPASS); // this line is the only difference to lowpass
  //flt.setApproximationMethod(rsPrototypeDesigner::BUTTERWORTH);
  //flt.setPrototypeOrder(order);

  //// apply filter:
  //int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  //rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

template<class TSig, class TTim, class TPar> // signal, time, parameter
void rsBiDirectionalFilter::applyButterworthLowpass(
  const TSig* x, const TTim* t, TSig* y, int N, TPar fc, int order, int numPasses, TPar gc)
{
  fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);
  rsNonUniformFilterIIR<TSig> flt;
  flt.setApproximationMethod(rsNonUniformFilterIIR<TSig>::ApproximationMethod::butterworth);
  flt.setFrequency(fc);
  flt.setOrder(order);

  //int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  // make such an estimation function for the non-uniform filter - it should probably return a
  // value in seconds
  int P = 1000;
  // ad-hoc - later we need to use something based a ringing time estimate - maybe something like
  // k * averageSampleRate / cutoff where k is a constant that depends on the ringing behavior of
  // the filter (Butterwort in this case - but we may later make it more general)

  rsApplyBiDirectionally(x, t, y, N, flt, P, numPasses);
}

template<class TSig, class TPar>
void rsBiDirectionalFilter::applyLowpass(const TSig *x, TSig *y, int N, TPar fc, TPar fs, int numPasses,
  TPar gc)
{
  // compute desired cutoff for single-pass filter:
  fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, 1, gc);

  // create and set up the filter:
  rsOnePoleFilter<TSig, TPar> flt;
  flt.setMode(rsOnePoleFilter<TSig, TPar>::LOWPASS_IIT); // impulse invariant - maybe switch to bilinear later
  flt.setSampleRate(fs);
  flt.setCutoff(fc);

  // apply filter:
  int P = getPaddingLength(fc, fs);
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

//=================================================================================================

template<class T>
bool rsZeroCrossingFinder::isUpwardCrossing(T *x, int n)
{
  if(x[n] < 0.0 && x[n+1] >= 0.0)
    return true;
  return false;

  // what about n < 0 or n >= N-2...this is unsafe
}

template<class T>
int rsZeroCrossingFinder::closestUpwardCrossingLeft(T *x, int N, int n)
{
  for(int i = n; i >= 0; i--)
    if(isUpwardCrossing(x, i))
      return i;
  return -1;
}

template<class T>
int rsZeroCrossingFinder::closestUpwardCrossingRight(T *x, int N, int n)
{
  for(int i = n; i < N-1; i++)
    if(isUpwardCrossing(x, i))
      return i;
  return -1;
}

template<class T>
int rsZeroCrossingFinder::closestUpwardCrossing(T *x, int N, int n)
{
  int cl = closestUpwardCrossingLeft( x, N, n);
  int cr = closestUpwardCrossingRight(x, N, n);
  if( cr-n < n-cl )
    return cr;
  else
    return cl;
}

template<class T>
T rsZeroCrossingFinder::upwardCrossingFrac(T *x, int N, int n, int p)
{
  T f = x[n]/(x[n]-x[n+1]);      // fractional part - init to zero of linear interpolant
  int q = rsMin(p, n, N-n-2);    // p, restricted to avoid access violation
  if(q > 0) {
    // refine linear zero estimate by Newton iteration on a higher order interpolating
    // polynomial using the zero of the linear interpolant as initial guess:
    T *a = new T[2*p+2];   // polynomial coefficients for interpolant (maybe use a static array)
    rsPolynomial<T>::interpolant(a, -q, 1, &x[n-q], 2*q+2);
    f = rsPolynomial<T>::rootNear(f, a, 2*q+1, 0.0, 1.0);
    delete[] a;
  }
  return f;
}

template<class T>
int rsZeroCrossingFinder::numUpwardCrossings(T *x, int N)
{
  int nz = 0;
  for(int n = 0; n < N-1; n++)
  {
    if( isUpwardCrossing(x, n) )
      nz++;
  }
  return nz;
}

template<class T>
std::vector<int> rsZeroCrossingFinder::upwardCrossingsInt(T *x, int N)
{
  int Nz = numUpwardCrossings(x, N);  // number of zero crossings in x
  int nz = 0;                         // index of the current zero crossing
  std::vector<int> z;
  z.resize(Nz);
  for(int n = 0; n < N-1; n++)
  {
    if(isUpwardCrossing(x, n))
    {
      z[nz] = n;
      nz++;
    }
  }
  return z;
}

template<class T>
std::vector<rsFractionalIndex> rsZeroCrossingFinder::upwardCrossingsIntFrac(T *x,
  int N, int p)
{
  std::vector<int> zi = upwardCrossingsInt(x, N);
  int Nz = (int)zi.size();
  std::vector<rsFractionalIndex> z;   // the positions of the zero crossings
  z.resize(Nz);
  for(int nz = 0; nz < Nz; nz++) {
    z[nz].intPart  = zi[nz];                                // store integer part
    z[nz].fracPart = upwardCrossingFrac(x, N, zi[nz], p);   // store fractional part
  }
  return z;
}

template<class T>
std::vector<T> rsZeroCrossingFinder::upwardCrossings(T *x, int N, int p)
{
  std::vector<rsFractionalIndex> z = upwardCrossingsIntFrac(x, N, p);
  std::vector<T> zd;
  zd.resize(z.size());
  for(unsigned int i = 0; i < z.size(); i++)
    zd[i] = z[i].intPart + z[i].fracPart; // combine integer and fractional part into a "double" value
  return zd;
}

template<class T>
std::vector<T> rsZeroCrossingFinder::bandpassedUpwardCrossings(T *x, int N, T fc,
  T bw, T fs, int np, int p)
{
  T *y = new T[N];
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y, N, fc, bw, fs, np);
  std::vector<T> z = upwardCrossings(y, N, p);
  delete[] y;
  return z;
}

//=================================================================================================

// move to rsArray, maybe get rid of using rsPolynomial<T>::fitQuadratic (formula simplifies
// because x-values are simple)
template<class T>
T rsMaxPosition(T* buffer, int N)
{
  int i = rsArray::maxIndex(buffer, N);

  if(i < 1 || i >= N-1)
    return i;

  T x[3] = { -1, 0, 1 };
  T a[3];
  rsPolynomial<T>::fitQuadratic(a, x, &buffer[i-1]);

  T offset = 0;
  if( abs(2*a[2]) >= abs(a[1]) ) // abs(offset) shall be <= 1
    offset = -a[1] / (2*a[2]);   // maximum of the parabola (zero-crossing of its slope)

  //return i; // for test
  return i + offset;  // also test
  //return i - offset;
  // figure out first, if it has to be + or - (i think +), then check, if the added constant in
  // the delta computation in line 302 should be 0.5, 1, 1.5 (or something else)
}

template<class T>
rsCycleMarkFinder<T>::rsCycleMarkFinder(T sampleRate, T minFundamental, T maxFundamental)
{
  fs   = sampleRate;
  fMin = minFundamental;
  fMax = maxFundamental;
}

template<class T>
void rsCycleMarkFinder<T>::refineCycleMarksByCorrelation(T *x, int N, std::vector<T>& cm, T f0)
{
  T* y = new T[N];
  if(correlationHighpass > 0)
    rsBiDirectionalFilter::applyButterworthHighpass(x, y, N, f0*correlationHighpass, fs, 4, 1);
  else
    rsArray::copy(x, y, N);

  // bug: i cannot refine the cycle marks in place inside the passed array because the refined mark
  // may overrun the unrefined, leading to a negative length in the process
  // try to devise an experiment that exposes this behavior (maybe with tow frequencies at
  // 100 and 200*1.618...

  int maxLength = 5000; // preliminary - use something based on the maximum time-delta between the
                        // cycle-marks in cm array

  int left = (int) cm[0];
  std::vector<T> cl(maxLength), cr(maxLength), corr(2*maxLength-1);
  for(unsigned int i = 1; i < cm.size(); i++)
  {
    int right      = (int) cm[i];
    int halfLength = rsFloorInt(correlationLength * 0.5 * (right-left));
    int length     = 2*halfLength;

    // prepare buffers for correlation computation:
    cl.resize(length); // maybe use raw arrays instead of vectors, avoid memory re-allocations
    cr.resize(length);
    corr.resize(2*length-1);
    rsArray::copySection(y, N, &cl[0],  left-halfLength, length);
    rsArray::copySection(y, N, &cr[0], right-halfLength, length);

    // use a matched filter to do the correlation,
    // see here: https://en.wikipedia.org/wiki/Matched_filter
    rsArray::reverse(&cl[0], length);                            // reversed left cycle is used as "impulse response"
    rsArray::convolve(&cr[0], length, &cl[0], length, &corr[0]); // OPTIMIZE: use FFT convolution
    T delta = rsMaxPosition(&corr[0], 2*length-1) - length + 1.5; // is +1.5 correct? was found by trial and error
    //delta = rsMaxPosition(&corr[0], 2*length-1) - length + 1.0; // test

    // overwrite array entry and update for next iteration:
    cm[i] = right + delta;
    left  = (int) cm[i];
  }

  delete[] y;
}

template<class T>
std::vector<T> rsCycleMarkFinder<T>::findCycleMarks(T *x, int N)
{
  if(algo == F0_ZERO_CROSSINGS)
    return findCycleMarksByFundamentalZeros(x, N);
  else if(algo == CYCLE_CORRELATION_OLD)         // deprecated...
    return findCycleMarksByCorrelationOld(x, N); // ...remove soon
  else
    return findCycleMarksByRefinement(x, N);

  // todo: maybe work with the more precise version that splits integer and fractional parts
  // of the zero-crossings

  rsError("Unknown algorithm setting");
  return std::vector<T>();
}

template<class T>
std::vector<T> rsCycleMarkFinder<T>::findCycleMarksByFundamentalZeros(T* x, int N)
{
  // Get initial estimate of fundamental by using an autocorrelation based algorithm at the center
  // of the input signal:
  T f0 = getFundamental(x, N);
  // here, we have a mutual dependency between rsInstantaneousFundamentalEstimator and
  // rsCycleMarkFinder - maybe break up by dragging estimateFundamentalAt out of the class

  T bw = bandPassWidth*f0; // absolute bandwidth
  T *y = new T[N];
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y, N, f0, bw, fs, bandpassSteepness);
  //rsPlotArrays(N, x, y); // plot original an bandpassed signals
  std::vector<T> z = rsZeroCrossingFinder::upwardCrossings(y, N, precision);
  delete[] y;
  return z;
}

template<class T>
std::vector<T> rsCycleMarkFinder<T>::findCycleMarksByRefinement(T* x, int N)
{
  // select a sample index n in the middle and estimate frequency there:
  //int nCenter = N/2;
  int nCenter = rsZeroCrossingFinder::closestUpwardCrossing(x, N, N/2);
  T f0 = getFundamental(x, N);
  T p = fs/f0; // curent period estimate (currently at nCenter)

  // create temporary signal to work with:
  std::vector<T> y(N);
  if(correlationHighpass > 0)
    rsBiDirectionalFilter::applyButterworthHighpass(
      x, &y[0], N, f0*correlationHighpass, fs, 4, 1);
  else
    rsArray::copy(x, &y[0], N);
  // rename correlationHighpass to refinementHighpass

  // nCenter serves as the initial cycle mark - from there, find the next one to the left by
  // correlating two segments of length p (or maybe 2*pn with windowing? - maybe experiment)
  std::vector<T> z;
  z.reserve((int) (2*ceil(N/p))); // estimated size needed is N/p, so twice that value should be more than enough

  // find cycle-marks to the left of nCenter by correlating two segments of length p (which is the
  // current estimate of the period):
  T right  = nCenter + rsZeroCrossingFinder::upwardCrossingFrac(&y[0], N, nCenter, precision);
  T left   = right - p;
  T length = right-left;
  z.push_back(right);  // nCenter = right serves as the initial cycle mark
  while(true)
  {
    left = refineCycleMark(&y[0], N, right, left); // new, refined left mark
    if(left == T(-1))
      break;
    z.push_back(left);
    length = right - left;
    right  = left;
    left   = right - length;
    if(left <= 0)
      break;
  }
  rsArray::reverse(&z[0], (int) z.size()); // bring marks in 1st half into ascending order

  // find cycle-marks to the right of nCenter (which is currently the last element in z):
  left  = z[z.size()-1];
  length = z[z.size()-1] - z[z.size()-2];
  right = left + length;
  while(true)
  {
    right = refineCycleMark(&y[0], N, left, right); // new, refined right mark
    if(right == T(-1))
      break;
    z.push_back(right);
    length = right - left;
    left   = right;
    right  = left + length;
    if(right >= N-1)
      break;
  }
  // hangs with Rhodes_F3.wav - still true?


  //rsPlotSignalWithMarkers(x, N, &z[0], (int) z.size());
  //rsPlotVector(rsDifference(z)); // plot the cycle-lengths

  return z;
}

template<class T>
std::vector<T> rsCycleMarkFinder<T>::findCycleMarksByCorrelationOld(T* x, int N)
{
  std::vector<T> z = findCycleMarksByFundamentalZeros(x, N); // initial estimates
  T f0 = getFundamental(x, N);
  refineCycleMarksByCorrelation(x, N, z, f0);
  return z;
}

template<class T>
T rsCycleMarkFinder<T>::getFundamental(T* x, int N)
{
  if(fundamental != T(0))
    return fundamental;
  else
  {
    int nCenter = rsZeroCrossingFinder::closestUpwardCrossing(x, N, N/2);
    return rsInstantaneousFundamentalEstimator<T>::
      estimateFundamentalAt(x, N, nCenter, fs, fMin, fMax);
  }
}


// introspection:

template<class T>
typename rsCycleMarkFinder<T>::ErrorMeasures
rsCycleMarkFinder<T>::getErrorMeasures(const std::vector<T>& cycleMarks, T period)
{
  int N = (int) cycleMarks.size();
  ErrorMeasures errors;
  errors.mean   = period - rsArray::meanDifference(&cycleMarks[0], N);
  errors.min    = T(0);
  errors.max    = T(0);
  errors.maxAbs = T(0);
  for(int i = 1; i < N; i++) {
    T estimate = cycleMarks[i] - cycleMarks[i-1];
    T error = period - estimate;
    if(error < errors.min)
      errors.min = error;
    if(error > errors.max)
      errors.max = error;
    if(abs(error) > errors.maxAbs)
      errors.maxAbs = abs(error);
  }

  // convert error measures to percent-of-true-value:
  T scale = T(100) / period;
  errors.mean   *= scale;
  errors.min    *= scale;
  errors.max    *= scale;
  errors.maxAbs *= scale;
  return errors;
}


// internal functions

template<class T>
T rsCycleMarkFinder<T>::periodErrorByCorrelation(T* x, int N, T left, T right)
{
  // a sort of crude, preliminary way to avoid producing garbage data or even crashes in
  // fade-in/out sections:
  T minPeriod = fs/fMax;
  if(right - left < minPeriod)
    return T(0);

  int iLeft  = rsRoundToInt(left);
  int iRight = rsRoundToInt(right);
  T error = periodErrorByCorrelation(x, N, iLeft, iRight);
  error  -= (right-left) - T(iRight-iLeft);
  return error;
}

template<class T>
T rsCycleMarkFinder<T>::periodErrorByCorrelation(T* x, int N, int left, int right)
{
  rsAssert(left  >= 0);
  rsAssert(right <  N);
  rsAssert(right >  left);
  int halfLength = rsFloorInt(correlationLength * 0.5 * (right-left));
  int length     = 2*halfLength;
  rsAssert(length >= 2);

  // prepare buffers for correlation computation:
  cl.resize(length);
  cr.resize(length);
  corr.resize(2*length-1);
  rsArray::copySection(x, N, &cl[0],  left-halfLength, length);
  rsArray::copySection(x, N, &cr[0], right-halfLength, length);
  applyWindow(&cl[0], length);
  applyWindow(&cr[0], length);

  // use a matched filter to get the cross-correlation sequence,
  // see here: https://en.wikipedia.org/wiki/Matched_filter
  // old:
  rsArray::reverse(&cl[0], length);                            // reversed left cycle is used as "impulse response"
  rsArray::convolve(&cr[0], length, &cl[0], length, &corr[0]); // OPTIMIZE: use FFT convolution
  //deBiasConvolutionResult(&corr[0], length);   // test

  // new:
  //rsCrossCorrelation(&cl[0], &cr[0], length, &corr[0], false);
  // doesn't seem to work -> make unit test -> ah, it produces a length N output instead of 2*N-1


#ifdef RS_DEBUG_PLOTTING
  // plot the signal chunks and correlation sequence:
  GNUPlotter plt; // #define DEBUG_PLOTTING in rapt.h to make it work
  rsArray::reverse(&cl[0], length);                     // reverse again for better visual match
  rsArray::scale(&corr[0], 2*length-1, 1./halfLength);  // to make the same scale as inputs
  plt.addDataArrays(length, &cl[0]);
  plt.addDataArrays(length, &cr[0]);
  plt.addDataArrays(2*length-1, &corr[0]);
  plt.plot();
#endif

  // from the cross-correlation, find the correlation lag for a best match:
  T maxPos  = rsMaxPosition(&corr[0], 2*length-1);
  T bestLag = maxPos + 1;       // found by trial and error - but why the +1?
  //T period  = bestLag / correlationLength;
  T period = bestLag * T(right-left) / T(length);
  // seems like dividing by correlationLength is better than dividing by the actual length-ratio
  // T(length) / T(right-left) which might be slightly different due to rounding to integers when
  // computing period (tried with a sine with period length of 45.5 cycles) ..or well, not seems
  // the exact ratio is better...makes more sense - but maybe more tests are needed

  return period - T(right-left);
}
// more ideas for improvement:
// -use de-biased cross-correlation sequence
// -fit quartic instead of quadratic for finding the subsample-precision maximum
// -apply a nonlinearity before fitting the parabola (it should be such that the found maximum
//  stays invariant with respect to scaling the signal - maybe log qualifies, maybe powers, too?)
// -use zero-crossings of the original signal near a correlation maximum (zeros are better than
//  just any value because they are (almost) invariant with respect to amplitude changes between
//  the cycles)

template<class T>
void rsCycleMarkFinder<T>::applyWindow(T* x, int N)
{
  T halfN = T(0.5) * T(N);
  for(int n = 0; n < N; n++)
  {
    x[n] *= rsWindowFunction::cosineSquared(T(n)-halfN, T(N));  // seems best among those that were tested
    //x[n] *= rsWindowFunction::rsRaisedCosineWindow(T(n)-halfN, T(N), 0.08);  // Hamming - produces large spikes of the error
    //x[n] *= rsWindowFunction::rsRaisedCosineWindow(T(n)-halfN, T(N), 0.07672);  // equiripple sidelobes - also spikes
    //x[n] *= rsWindowFunction::rsExactBlackmanWindow(T(n)-halfN, T(N));  // similar to rsCosineSquaredWindow
    // maybe try a window that has zeros for the first and second derivative at +-1
  }

  // todo: optimize (using some kind of rsWindowIterator class), flexibilize (allow client to
  // select window)
}

template<class T>
void rsCycleMarkFinder<T>::deBiasConvolutionResult(T* x, int N)
{
  for(int n = 0; n < N; n++)  {
    T scale = T(N) / T(N-n);
    x[N-1+n] *= scale;
    x[N-1-n] *= scale;
  }
}
// experimental - multiplies the result of a convolution with a sequence that compesates for the
// number of summed terms - the center value is due to a sum of N nonzero terms, left and right of
// it, only N-1 terms were summed, and so on

template<class T>
T rsCycleMarkFinder<T>::autoCorrelation(T* x, int N, int n1, int n2, int M)
{
  //return rsCrossCorrelation(&x[n1], rsMin(M, N-n1), &x[n2], rsMin(M, N-n2));
    // test - should produce same results as code below - nope - why?

  T xx = 0;
  T xy = 0;
  T yy = 0;
  int i, j;
  for(int n = 0; n < M; n++)
  {
    i = n1+n;
    j = n2+n;
    if(i < 0 || i >= N || j < 0 || j >= N)
      continue;
    // todo: optimize by precomputing nMin, nMax and get rid of the if inside the loop
    // do unit tests comparing against prototype with the "if"

    xx += x[i] * x[i];
    xy += x[i] * x[j];
    yy += x[j] * x[j];
  }
  //return xy;  // old - unnormalized sum-of-products
  if(xx == 0 || yy == 0)
    return 0;
  return xy / sqrt(xx*yy); // new - normalized autocorrelation - better results
}

template<class T>
T rsCycleMarkFinder<T>::bestMatchOffset(T* x, int N, int nFix, int nVar, int M)
{
  T left  = autoCorrelation(x, N, nFix, nVar-1, M);
  T mid   = autoCorrelation(x, N, nFix, nVar,   M);
  T right = autoCorrelation(x, N, nFix, nVar+1, M);
  while(right > mid)
  {
    nVar++;
    left  = mid;
    mid   = right;
    right = autoCorrelation(x, N, nFix, nVar+1, M);
  }
  while(left > mid)
  {
    nVar--;
    right = mid;
    mid   = left;
    left  = autoCorrelation(x, N, nFix, nVar-1, M);
  }
  //return T(nVar);  //preliminary

  // fit parabola
  T xp[3] = { -1, 0, 1 };
  T yp[3] = { left, mid, right };
  T a[3];
  rsPolynomial<T>::fitQuadratic(a, xp, yp);
  T offset = 0;


  if( abs(2*a[2]) >= abs(a[1]) ) // abs(offset) shall be <= 1
    offset = -a[1] / (2*a[2]);   // maximum of the parabola (zero-crossing of its slope)

  // catch cases where parabola becomes degenerate:
  if(rsAbs(a[2]) < RS_EPS(T))
    return T(nVar);
  // maybe i should investigate the issue some more - why does it become degenerate - should that
  // be possible to happen - for these investigations, just comment out this early return and
  // uncomment the debug-plotting code below:

  /*
  // debug:
  if(!rsIsFiniteNumber(offset)) {
    if(nVar > nFix)
      rsPlotArray(&x[nFix], nVar-nFix);
    else
      rsPlotArray(&x[nVar], nFix-nVar);
  }
  rsAssert(rsIsFiniteNumber(offset));
  // we hit this with the (long) rhodes sample - look into rsSinusoidalAnalyzer how to deal with it
  // it has to do with the parabolic interpolant becoming degenerate, i.e. a[2] == 0 ..actually,
  // in the case of the rhodes sample, a[1] is also zero - figure out, why the parabola becomes
  // degenerate in the first place ...it's because left, right, mid are all the same - but why are
  // they all the same?....
  */


  return T(nVar) + offset;
  //return T(nVar) - offset;  // test
}

template<class T>
T rsCycleMarkFinder<T>::bestMatchOffset(T* x, int N, T nFix, T nVar)
{
  // a sort of crude, preliminary way to avoid producing garbage data or even crashes in
  // fade-in/out sections:
  T minPeriod = fs/fMax;
  if(abs(nVar-nFix) < minPeriod)
    return T(0);

  int iFix = rsRoundToInt(nFix);
  int iVar = rsRoundToInt(nVar);
  int M    = (int) ceil(correlationLength*abs(iVar-iFix));

  T result = bestMatchOffset(x, N, iFix, iVar, M);
    // may have to add/or subtract nVar-iVar or something? otherwise may produce jitter?
    // -> experiment

  //result -= (nFix-iFix); // correct?
  result += (nFix-iFix); // correct?

  //rsAssert(rsIsFiniteNumber(result));

  return result;
}

template<class T>
T rsCycleMarkFinder<T>::refineCycleMark(T* x, int N, T anchor, T mark)
{
  //return bestMatchOffset(x, N, anchor, mark);  // rename to refineCrossCorr2 or sth.
  // preliminary - todo: switch between refinement algorithms

  switch(algo)
  {
  case CYCLE_CORRELATION: return bestMatchOffset(x, N, anchor, mark);
  case WINDOWED_CORRELATION:
  {
    if(anchor > mark)
      return mark - periodErrorByCorrelation(x, N, mark, anchor);
    else
      return mark + periodErrorByCorrelation(x, N, anchor, mark);
  }
  case ZERO_CROSSINGS:  return refineByZeroCrossing(x, N, anchor, mark);
  case CORRELATED_ZERO:
  {
    T tmp = bestMatchOffset(x, N, anchor, mark);
    return refineByZeroCrossing(x, N, anchor, tmp);
  }
  default:
  {
    rsError("Unknown algorithm");
    return mark;
  }
  }
}

template<class T>
T rsCycleMarkFinder<T>::refineByZeroCrossing(T* x, int N, T anchor, T mark)
{
  //T delta = mark - anchor; // for debug
  //rsAssert(abs(delta) > T(0));
  int n = rsRoundToInt(mark);
  n = rsZeroCrossingFinder::closestUpwardCrossing(x, N, n);
  if(n != -1)
  {
    //rsAssert(x[n] < 0 && x[n+1] >= 0);  // debug
    T newMark = n + rsZeroCrossingFinder::upwardCrossingFrac(x, N, n, precision);
    return newMark;
  }
  else
    return T(-1);  // -1 encodes failure to find a zero crossing
}

//=================================================================================================

// converts raw value from the cosine-generator into the window value - todo: write a class
// rsWindowFunctionIterator (subclass of rsSineIterator) that wraps that:
template<class T>
RS_INLINE T cosineToWindow(T c)
{
  //return T(1);  // for test: rectangular window

  // exact Blackman coefficients, modified such that we don't need to precompute c2=2*c*c-1:
  const T a0 = 6508.0/18608.0;
  const T a1 = 9240.0/18608.0;
  const T a2 = 2860.0/18608.0;
  return a0 + a1*c + a2*c*c;
}
// rename to cosineToBlackman
// todo: maybe write a general function that recursively computes cosines of successive multiples
// of an angle from the cosine of that angle. This recursion is based on the trigonometric
// identity:
// cos(a)*cos(b) = (cos(a-b) + cos(a+b))/2 which implies
// cos(w)*cos(w) = (cos(0)   + cos(2*w))/2, so cos(2*w) = 2*cos^2(w) - cos(0)
// the general recursion formula is:
// c[0] = cos(0*w) = 1
// c[1] = cos(1*w) = c
// c[n] = cos(n*w) = a*c[n-1]-c[n-2] with a = 2*cos(w) = 2*c
// from linear combinations of these successive cosine values, various windows can be created,
// see https://en.wikipedia.org/wiki/Window_function#Higher-order_generalized_cosine_windows






template<class TSig, class TPos>
RS_INLINE void sincInterpolatorLoop(int mMin, int mMax, TPos &tf, rsSineIterator<TPos> &sinIt,
  rsSineIterator<TPos> &wndIt, TSig &y, TSig *&x, int &ti, TPos &ws)
{
  TPos w;
  for(int m = mMin; m <= mMax; m++)
  {
    w   = sinIt.getValue() * cosineToWindow(wndIt.getValue()) / (m-tf);
    ws += w;
    y  += w * x[ti+m];
  }
}

template<class TSig, class TPos>
RS_INLINE void sincInterpolatorLoopNoStretch(int mMin, int mMax, TPos &tf, TPos &s,
  rsSineIterator<TPos> &wndIt, TSig &y, TSig *&x, int &ti, TPos &ws)
{
  TPos w;
  for(int m = mMin; m <= mMax; m++)
  {
    w  = s * cosineToWindow(wndIt.getValue())  / (m-tf);
    ws += w;
    y  += w * x[ti+m];
    s  *= -1.0;
  }
}
// todo: check, if types TSig/TPos are correct for all variables


template<class TSig, class TPos>
TSig rsResampler<TSig, TPos>::signalValueViaSincAt(TSig *x, int N, TPos t, TPos sincLength, TPos stretch)
{
  int L  = (int) floor(sincLength);
  int ti = (int) floor(t);         // integer part of t
  if( ti < 0 || ti >= N )
    return 0.0;
  TPos tf = t - ti;                // fractional part of t
  TSig y  = 0.0;                   // output value
  TPos s;                          // sine value
  TPos ws = 0.0;                   // sum of tap weights
  int mMin = -rsMin(L/2-1, ti);    // minimum shift
  int mMax = +rsMin(L/2, N-ti-1);  // maximum shift

  // optimized loop for stretch == 1.0 (used for downward transpositions):
  rsSineIterator<TPos> wndIt(2*PI/sincLength, 2*PI*(mMin-tf)/sincLength+PI/2);
  if( stretch == 1.0 )
  {
    s = sin(PI*(mMin-tf)/stretch) / PI;
    if( tf > RS_EPS(TPos) )
      sincInterpolatorLoopNoStretch(mMin, mMax, tf, s, wndIt, y, x, ti, ws);
    else
    {
      // split loop into a part for m < 0 and for m > 0, the case for m == 0 is handled separately
      // because it would lead to division-by-zero error:
      sincInterpolatorLoopNoStretch(mMin, -1, tf, s, wndIt, y, x, ti, ws);
      wndIt.getValue(); // for internal state-update
      ws += 1.0;
      y  += x[ti];
      sincInterpolatorLoopNoStretch(1, mMax, tf, s, wndIt, y, x, ti, ws);
    }
    return y/ws;
  }

  // general case with stretch (could handle special case above also, but less efficiently):
  rsSineIterator<TPos> sinIt(PI/stretch, PI*(mMin-tf)/stretch, 1.0/PI);
  if( tf > RS_EPS(TPos) )
    sincInterpolatorLoop(mMin, mMax, tf, sinIt, wndIt, y, x, ti, ws);
  else
  {
    sincInterpolatorLoop(mMin, -1, tf, sinIt, wndIt, y, x, ti, ws);
    sinIt.getValue();
    wndIt.getValue();
    ws += 1.0/stretch;
    y  += x[ti]/stretch;
    sincInterpolatorLoop(1, mMax, tf, sinIt, wndIt, y, x, ti, ws);
  }
  return y/ws;
}
// \todo: make a version of this function for periodic signals, maybe the loop should look like:
//  for(int m = sincLength/2; m <= sincLength/2; m++)
//    y += x[(ti+m)%N] * rsWindowedSinc(m-tf, sincLength, stretch);
// ...but we need to check, if this is correct. This version could be used for looped
// sample-playback (for example, in a sampler) or for wavetable-oscillators

template<class TSig, class TPos>
void rsResampler<TSig, TPos>::transposeLinear(TSig *x, int xN, TSig *y, int yN, TPos factor)
{
  int  nw;         // write position
  TPos nr = 0.0;   // read position
  int  nri;        // integer part of nr
  TPos nrf;        // fractional part of nr
  for(nw = 0; nw < yN; nw++)
  {
    nri = (int) floor(nr);

    if( nri+1 >= xN )
      break;  // end of input reached - todo: precompute the maximum nw and use it in the loop
              // condition to get rid of that conditional here

    nrf = nr - nri;
    y[nw] = (1.0-nrf)*x[nri] + nrf*x[nri+1];
    nr += factor;
  }
  rsArray::fillWithZeros(&y[nw], yN-nw); // fill tail with zeros
}

template<class TSig, class TPos>
void rsResampler<TSig, TPos>::transposeSinc(TSig *x, int xN, TSig *y, int yN, TPos factor,
  TPos sincLength, bool antiAlias)
{
  TPos stretch = 1.0;
  if( antiAlias == true )
    stretch = rsMax(1.0, factor);

  int  nw;         // write position
  TPos nr = 0.0;   // read position
  for(nw = 0; nw < yN; nw++)
  {
    y[nw] = signalValueViaSincAt(x, xN, nr, sincLength, stretch);
    nr += factor;
  }
  rsArray::fillWithZeros(&y[nw], yN-nw);
}

template<class TSig, class TPos>
void rsResampler<TSig, TPos>::shiftSinc(TSig *x, TSig *y, int N, TPos amount, TPos sincLength)
{
  if( x == y )
  {
    TSig *tmp = new TSig[N];
    shiftSinc(x, tmp, N, amount, sincLength);
    rsArray::copy(tmp, y, N);
    delete[] tmp;
  }
  else
  {
    for(int n = 0; n < N; n++)
      y[n] = rsResampler::signalValueViaSincAt(x, N, n-amount, sincLength, 1.0);
  }
}

//=================================================================================================

template<class TSig, class TPos>
void rsTimeWarper<TSig, TPos>::timeWarpSinc(TSig *x, int xN, TSig *y, TPos *w, int yN,
  TPos minSincLength, TPos maxLengthScaler, bool antiAlias)
{
  if( antiAlias == false )
  {
    for(int n = 0; n < yN; n++)
      y[n] = rsResampler<TSig, TPos>::signalValueViaSincAt(x, xN, w[n], minSincLength, 1.0);
  }
  else
  {
    TPos tOld = w[0];
    for(int n = 0; n < yN; n++)
    {
      TPos speed   = rsAbs(w[n]-tOld);
      TPos stretch = rsMax(1.0, speed);
      //TPos length  = minSincLength * rsMin(stretch, maxLengthScaler); // old, probably buggy
      TPos length  = rsMax(minSincLength, minSincLength*rsMin(stretch, maxLengthScaler));

      // maybe better, do just:
      // length = rsLimitToRange(minSincLength*stretch, minSincLength, maxSincLength);
      // where maxSincLength is precomputed as minSincLength*maxLengthScaler

      y[n] = rsResampler<TSig, TPos>::signalValueViaSincAt(x, xN, w[n], length, stretch);
      tOld = w[n];
    }
  }
}

template<class TSig, class TPos>
void rsTimeWarper<TSig, TPos>::invertMonotonousWarpMap(TPos *w, int N, TPos *wi)
{
  //bool cubic = true; // make this a parameter for the function

  int M = (int) ceil(w[N-1]); // length of inverse map
  int n, i, iOld = 0;
  for(n = 0; n < M; n++)
  {
    // find index i, such that ty[i] <= n and ty[i+1] >= n:
    for(i = iOld; i < N-1; i++)
    {
      if( w[i] <= n && w[i+1] >= n )
        break;
    }
    wi[n] = rsInterpolateLinear(w[i], w[i+1], (TPos)i, (TPos)(i+1), (TPos)n);


    /*
    if( cubic == true && i >= 1 && i < N-2 )
    {
      wi[n] = rsInterpolateCubicHermite(w[i-1], w[i], w[i+1], w[i+2], (TPos)(i-1), (TPos)i,
        (TPos)(i+1), (TPos)(i+2), (TPos)n);
    }
    else
      wi[n] = rsInterpolateLinear(w[i], w[i+1], (TPos)i, (TPos)(i+1), (TPos)n);
    */
    //wi[n] = rsInterpolateLinear(w[i], (TPos)i, w[i+1], (TPos)(i+1), (TPos)n);
      // Can be optimized: (i+1)-i == 1 - always, the function computes the value, but maybe it's
      // not worth it. On the other hand, maybe the quality could be improved by using cubic
      // hermite instead of linear interpolation.

    iOld = i;
  }
  // Remark: I have tried to use cubic hermite interpolation instead of linear, but it doesn't make
  // a difference quality wise (the spectral artifacts do not decrease measurably), so linear seems
  // to be appropriate. Moreover, with linear interpolation a forward/backward mapping roundtrip of
  // any time-value should be an indentity operation (up to roundoff error) whcih seems a desirable
  // property.
}

template<class TSig, class TPos>
int rsTimeWarper<TSig, TPos>::getPitchModulatedLength(TPos *r, int N)
{
  TPos s = 0.0;
  for(int n = 1; n < N; n++)
    s += 1.0 / r[n-1];
  return (int) ceil(s);
  // Remark: the last readout speed value r[N-1] is irrelevant because it would formally determine
  // the time-instant, where the (nonexistent) next-to-last sample x[N] would have to be written
}

template<class TSig, class TPos>
void rsTimeWarper<TSig, TPos>::applyPitchModulation(TSig *x, TPos *r, int N, TSig *y,
  TPos minSincLength, TPos maxLengthScaler, bool antiAlias)
{
  rsVariableSpeedPlayer<TSig, TPos> vsp;
  vsp.setInputAndSpeed(x, r, N);
  vsp.getOutput(y, minSincLength, maxLengthScaler, antiAlias);
}

//=================================================================================================

template<class TSig, class TPos>
rsVariableSpeedPlayer<TSig, TPos>::rsVariableSpeedPlayer()
{
  init();
}

template<class TSig, class TPos>
rsVariableSpeedPlayer<TSig, TPos>::~rsVariableSpeedPlayer()
{
  clear();
}

template<class TSig, class TPos>
void rsVariableSpeedPlayer<TSig, TPos>::setInputAndSpeed(TSig *input, TPos *r, int length)
{
  clear();

  x  = input;
  Nx = length;

  // Create inverse warping map (i.e., the map, that assigns for each input sample-index n in x a
  // (noninteger) time-instant t in y where the sample value x[n] should be written, such that
  // y(t) = x[n]:
  wi = new TPos[Nx];
  wi[0] = 0.0; // y(0.0) = x[0]
  for(int n = 1; n < Nx; n++)
    wi[n] = wi[n-1] + 1.0 / r[n-1];

  // Obtain desired warping map by inverting the wi map:
  Ny = (int) ceil(wi[Nx-1]);
  w  = new TPos[Ny];
  rsTimeWarper<TSig, TPos>::invertMonotonousWarpMap(wi, Nx, w);
}

template<class TSig, class TPos>
TPos rsVariableSpeedPlayer<TSig, TPos>::warpTime(TPos tx)
{
  return rsArray::interpolateClamped(wi, Nx, tx);
}

template<class TSig, class TPos>
TPos rsVariableSpeedPlayer<TSig, TPos>::unwarpTime(TPos ty)
{
  return rsArray::interpolateClamped(w, Ny, ty);
}

template<class TSig, class TPos>
void rsVariableSpeedPlayer<TSig, TPos>::getOutput(TSig *y, TPos minSincLength, TPos maxLengthScaler,
  bool antiAlias)
{
  rsTimeWarper<TSig, TPos>::timeWarpSinc(x, Nx, y, w, Ny, minSincLength, maxLengthScaler, antiAlias);
}

template<class TSig, class TPos>
std::vector<TSig> rsVariableSpeedPlayer<TSig, TPos>::getOutput(TPos minSincLength, TPos maxLengthScaler,
  bool antiAlias)
{
  std::vector<TSig> y(Ny);
  getOutput(&y[0], minSincLength, maxLengthScaler, antiAlias);
  return y;
}

template<class TSig, class TPos>
std::vector<TPos> rsVariableSpeedPlayer<TSig, TPos>::getTimeWarpMapXY()
{
  std::vector<TPos> map(Nx);
  rsArray::copy(wi, &map[0], Nx);
  return map;
}

template<class TSig, class TPos>
std::vector<TPos> rsVariableSpeedPlayer<TSig, TPos>::getTimeWarpMapYX()
{
  std::vector<TPos> map(Ny);
  rsArray::copy(w, &map[0], Ny);
  return map;
}

template<class TSig, class TPos>
std::vector<TPos> rsVariableSpeedPlayer<TSig, TPos>::invertSpeeds(std::vector<TPos>& speeds)
{
  rsVariableSpeedPlayer<TSig, TPos> vsp;
  vsp.setInputAndSpeed(nullptr, &speeds[0], (int)speeds.size());
  std::vector<TPos> map = vsp.getTimeWarpMapYX();
  int N = (int)map.size();
  for(int n = 0; n < N-1; n++)
    map[n] = 1.0 / (map[n+1]-map[n]);
  map[N-1] = 1;   // value irrelevant (not used in map computation)
  return map;
}

template<class TSig, class TPos>
std::vector<TSig> rsVariableSpeedPlayer<TSig, TPos>::applyPlaybackSpeed(std::vector<TSig>& x,
  std::vector<TPos>& s)
{
  rsAssert(x.size() == s.size());
  rsVariableSpeedPlayer<TSig, TPos> vsp;
  vsp.setInputAndSpeed(&x[0], &s[0], (int)x.size());
  return vsp.getOutput();
}

template<class TSig, class TPos>
void rsVariableSpeedPlayer<TSig, TPos>::init()
{
  x = w = wi = nullptr;
  Nx = Ny = 0;
}

template<class TSig, class TPos>
void rsVariableSpeedPlayer<TSig, TPos>::clear()
{
  delete[] wi;
  delete[] w;
  w = wi = nullptr;
}

//=================================================================================================

template<class TSig, class TPos>
void rsPitchFlattener<TSig, TPos>::setInput(TSig *x, TPos *f, int N, TPos ft)
{
  // x: input, f: instantaneous frequencies, N: length of x and f, ft: target frequency

  // set default target frequency, if 0 is passed:
  if(ft == 0.0)
    ft = rsArray::mean(f, N);

  // create temporary read-speed array:
  TPos *r = new TPos[N];
  for(int n = 0; n < N; n++)
    r[n] = ft / f[n];

  // setup inherited variable speed player and clean up temporary array:
  rsVariableSpeedPlayer<TSig, TPos>::setInputAndSpeed(x, r, N);
  delete[] r;
}

//=================================================================================================

template<class TSig, class TPos>
void rsPhaseLockedCrossfader<TSig, TPos>::setInputs(TSig *in1, TPos *f1, int len1, TSig *in2,
  TPos *f2, int len2, TPos ft)
{
  x1 = in1;
  x2 = in2;
  N1 = len1;
  N2 = len2;
  if(ft == 0.0)
    ft = 0.5 * (rsArray::mean(f1, N1) + rsArray::mean(f2, N2));
  pf1.setInput(x1, f1, N1, ft);
  pf2.setInput(x2, f2, N2, ft);
}

template<class TSig, class TPos>
void rsPhaseLockedCrossfader<TSig, TPos>::setFlattenedCrossfade(TPos start, TPos end, TPos shift)
{
  // adjust crossfade-start such that it is an integer with respect to x1 and crossfade-end to
  // be an integer with respect to x2 for seamless splicing:
  cs1   = rsRoundToInt((TPos)pf1.unwarpTime(start));
  start = pf1.warpTime(cs1);
  cs2   = pf2.unwarpTime(start-shift);
  ce2   = rsRoundToInt((TPos)pf2.unwarpTime(end-shift));
  end   = pf2.warpTime(ce2) + shift;
  ce1   = pf1.unwarpTime(end);
  this->shift = shift;
}

template<class TSig, class TPos>
int rsPhaseLockedCrossfader<TSig, TPos>::getCrossfadeOutputLength()
{
  TPos cl1 = ce1 - cs1;              // crossfade length in x1
  TPos cl2 = ce2 - cs2;              // crossfade length in x2
  TPos cly;                          // crossfade length in y
  //cly = 0.5*(cl1+cl2);               // arithmetic mean
  cly = 1 / (0.5*(1/cl1 + 1/cl2));     // harmonic mean
  return rsRoundToInt(cly);

  // Using arithmetic mean, the vibrato speed during crossfade will be between the vibrato speeds
  // of x1, x2. Using the harmonic mean, the output frequency will be the arithmetic mean between
  // the input frequencies of x1, x2. The geometric mean is in between. Maybe we can use a
  // generalized mean and expose the exponent as user parameter.
  // https://en.wikipedia.org/wiki/Generalized_mean
}

template<class TSig, class TPos>
std::vector<TSig> rsPhaseLockedCrossfader<TSig, TPos>::getFlattenedSignal1()
{
  return pf1.getOutput();
}

template<class TSig, class TPos>
std::vector<TSig> rsPhaseLockedCrossfader<TSig, TPos>::getFlattenedSignal2()
{
  return pf2.getOutput();
}

template<class TSig, class TPos>
std::vector<TSig> rsPhaseLockedCrossfader<TSig, TPos>::getOutput()
{
  // splice together heading section of x1, crossfade section and trailing section of x2:
  std::vector<TSig> yc = getCrossfadeSection();
  int L  = (int)yc.size();
  int s2 = ce2+1;             // start in x2
  int N  = cs1 + L + N2-s2;   // output length
  std::vector<TSig> y(N);
  int n;
  for(n = 0; n < cs1; n++)
    y[n] = x1[n];
  for(n = cs1; n < cs1+L; n++)
    y[n] = yc[n-cs1];
  for(n = cs1+L; n < N; n++)
    y[n] = x2[n-cs1-L+s2];
  return y;
}

template<class TSig, class TPos>
std::vector<TPos> rsPhaseLockedCrossfader<TSig, TPos>::getTimeWarpMapXY1()
{
  return pf1.getTimeWarpMapXY();
}

template<class TSig, class TPos>
std::vector<TPos> rsPhaseLockedCrossfader<TSig, TPos>::getTimeWarpMapYX1()
{
  return pf1.getTimeWarpMapYX();
}

template<class TSig, class TPos>
std::vector<TPos> rsPhaseLockedCrossfader<TSig, TPos>::getTimeWarpMapXY2()
{
  return pf2.getTimeWarpMapXY();
}

template<class TSig, class TPos>
std::vector<TPos> rsPhaseLockedCrossfader<TSig, TPos>::getTimeWarpMapYX2()
{
  return pf2.getTimeWarpMapYX();
}

template<class TSig, class TPos>
std::vector<TSig> rsPhaseLockedCrossfader<TSig, TPos>::getCrossfadeSection()
{
  computeReadoutTimes();
  int  L     = (int)t1.size();
  TSig scale = 1.0 / (L-1);
  TSig y1, y2, c;
  std::vector<TSig> yc(L);
  for(int n = 0; n < L; n++)
  {
    c  = scale * n;
    y1 = rsResampler<TSig, TPos>::signalValueViaSincAt(x1, N1, t1[n], 64, 1.0);
    y2 = rsResampler<TSig, TPos>::signalValueViaSincAt(x2, N2, t2[n], 64, 1.0);
    yc[n] = (1-c)*y1 + c*y2;
  }
  return yc;
  // maybe later include optional anti-aliasing
}

template<class TSig, class TPos>
void rsPhaseLockedCrossfader<TSig, TPos>::computeReadoutTimes()
{
  TPos cl1 = ce1 - cs1;                        // crossfade length in x1
  TPos cl2 = ce2 - cs2;                        // crossfade length in x2
  int  L   = getCrossfadeOutputLength();       // crossfade length in output
  TPos scl = 1.0 / (L-1);                      // scaler to map 0..L-1 to 0..1
  TPos t, c;                                   // t: normalized time 0..1, c(t)
  TPos a = 1, b = 0;                           // parameters of c(t) = a*t + b*t^2
  bool sweep = true;                           // maybe make this a member variable
  if(sweep == true)
  {
    TPos s0 = L/cl1;                           // slope of c(t) at t=0
    TPos s1 = L/cl2;                           // slope of c(t) at t=1
    a = s0;
    b = 0.5*(s1-a);
  }
  t1.resize(L);
  t2.resize(L);
  for(int n = 0; n < L; n++)
  {
    t = scl * n;
    c = a*t + b*t*t;
    TPos tx1 = pf1.warpTime(cs1+c*cl1);
    TPos tx2 = pf2.warpTime(cs2+c*cl2)+shift;
    TPos txw = (1-t)*tx1 + t*tx2;              // aggree on warped readout time
    t1[n] = pf1.unwarpTime(txw);
    t2[n] = pf2.unwarpTime(txw-shift);
  }

  // Notes:
  // The boolean "sweep" variable could at some point be exposed as user option. If sweep mode
  // is off and the two input signals have different fixed frequencies, we will see a fixed output
  // frequency during the crossfade which is intermediate between the two input frequencies. Which
  // frequency exactly that is, is determined by the return value of getCrossfadeOutputLength - if
  // it returns the harmonic mean between both input crossfade-lengths, the frequency during
  // crossfade will be the arithmetic mean between both input frequencies (and vice versa).
  // The computation of txw as weighted mean between tx1 and tx2: txw = (1-t)*tx1 + t*tx2 could
  // use also "c" instead of "t" as weight - or some combination of t and c and we could also try
  // weighted means different from the arithmetic mean. All of these modifications will lead to
  // different behavior when both input signals have vibratos with different speeds.
}

//=================================================================================================

template<class T>
void rsInstantaneousFundamentalEstimator<T>::estimateReliability(T *x, int N,
  const std::vector<T>& z, T *r)
{
  int Nz = (int) z.size(); // number of zero crossings in z
  int nz;                  // index of zero-crossing
  int sl, sr, el, er;      // start and end, left and right
  T c;                     // cross-correlation value

  sr = rsCeilInt( z[0]);
  er = rsFloorInt(z[1]);
  for(nz = 2; nz < Nz; nz++)
  {
    sl = sr;
    el = er;
    sr = el + 1;
    er = (int) z[nz];

    //T test = z[nz]; // for debug
      // seems like the last zero crossing is too close to the end of the signal such that the
      // zero crossing detection accesses invalid array indices - or something

    //c  = rsCrossCorrelation(&x[sl], el-sl+1, &x[sr], er-sr+1);
    c  = rsStretchedCrossCorrelation(&x[sl], el-sl+1, &x[sr], er-sr+1);
    rsArray::fillWithValue(&r[sl], el-sl+1, c);
  }

  // extend first and last reliability value to the start and end of the r-array:
  rsArray::fillWithValue(&r[el], N-el, c);
  sl = rsCeilInt(z[0]);
  rsArray::fillWithValue(r, sl, r[sl]);
}

template<class T>
void rsInstantaneousFundamentalEstimator<T>::measureInstantaneousFundamental(T *x, T *f,
  int N, T fs, T fMin, T fMax, T *r, int cycleMarkAlgo)
{
  rsCycleMarkFinder<T> cmf(fs, fMin, fMax); // todo: maybe set it up - or maybe have it a member and allow client code to set it up
  std::vector<T> z = cmf.findCycleMarks(x, N);

  // old:
  //std::vector<double> z = rsCycleMarkFinder::findCycleMarks(x, N, fs, fMin, fMax, cycleMarkAlgo, 3);
    // 3 is the precision - may be tweaked

  // Measure periods. They are given by the distance between successive zero-crossings. The time
  // instant at which we consider this period measurement to be effective is midway between the
  // two zero-crossings:
  int Nz = (int) z.size();
  T *p  = new T[Nz-1];
  T *tp = new T[Nz-1];
  int n;
  for(n = 0; n < Nz-1; n++)
  {
    p[n]  = z[n+1]-z[n];
    tp[n] = 0.5*(z[n+1]+z[n]);
  }

  // Interpolate period measurements up to samplerate and convert to frequencies
  T *tn = new T[N];
  rsArray::fillWithIndex(tn, N);
  rsInterpolateSpline(tp, p, Nz-1, tn, f, N, 1);
  for(n = 0; n < N; n++)
    f[n] = fs/f[n];

  // maybe a better way of interpolation would be to compute an instantaneous phase for each
  // zero-crossing/pitch-mark (each would have a value of phase[n] = 2*pi*n), interpolate these
  // phases, and obtain the instantaneous frequencies as difference - this would ensure a "phase
  // coherent interpolation" of the instantaneous frequencies (i.e. the instantaneous phase for a
  // resynthesized signal will not tend to drift depending on the specifics of the interpolation
  // scheme)
  // we should not really extrapolate the head and tail values with the first and last polynomial
  // repeating the last value will probably be better

  // If more actual data values are needed before the interpolation, we could find the zeros of
  // the derivatives which would give values halfway in between our current values. If the envelope
  // is exponential, the distance between maxima also equals the period (see Physics of Musical
  // Instruments, 2nd Ed, p.12)

  // assessment of reliability, if desired:
  if( r != nullptr )
    estimateReliability(x, N, z, r);

  // Cleanup:
  delete[] p;
  delete[] tp;
  delete[] tn;
}

template<class T>
T rsInstantaneousFundamentalEstimator<T>::estimateFundamentalAt(T *x, int N, int n,
  T fs, T fMin, T fMax)
{
  T pMax = fs/fMin; // maximum detectable period

  // Compute desired Length of autocorrelation sequence. It should be significantly longer than the
  // lag that corresponds to our maximum detectable period because the accuracy of the measured
  // values degrades for lags that are close to the end of the sequence due to less samples used in
  // the averaging. We need it to be odd so we can use a chunk that is centered around n.
  T k = 2.0;
  int L = (int) ceil(k*pMax);
  if( rsIsEven(L) )
    L += 1;

  // Get a chunk from the signal, typically centered around n, but the center is shifted at the
  // start and end of the input signal, such that we don't need to zero-extend the input
  // (conceptually):
  T *y = new T[L];                      // chunk from input signal
  int ns = rsMin(rsMax(0, n-L/2), N-L); // start of the chunk
  rsArray::copySection(x, N, y, ns, L);
  //rsPlotArray(y, L);

  // measure frequency:
  T r;  // reliability
  T f = rsAutoCorrelationPitchDetector<T>::estimateFundamental(y, L, fs, fMin, fMax, &r);

  // cleanup and return:
  delete[] y;
  return f;
}



//=================================================================================================

template<class T>
std::vector<int> rsPeakPicker<T>::getPeakIndices(const T* x, int N) const
{
  // todo: if smoothing is selected, create a pre-smoothed copy of x and operate on that

  std::vector<int> peaks;
  for(int i = 0; i < N; i++)
    if(isRelevantPeak(i, x, N))
      peaks.push_back(i);
  return peaks;
}

template<class T>
bool rsPeakPicker<T>::isRelevantPeak(int i, const T* x, int N) const
{
  return false;  // preliminary

  // this function is supposed to check all the criteria that must be met for a relevant peak
}

//=================================================================================================

//template<class T>
//void rsEnvelopeExtractor<T>::extractEnvelope(const T* x, int N, T* env)
//{
//
//}

template<class T>
void rsEnvelopeExtractor<T>::sineEnvelopeWithDeBeating(const T* x, int N, T* env)
{
  typedef std::vector<double> Vec;

  // get raw envelope:
  Vec envTime, envValue;
  getAmpEnvelope(x, N, envTime, envValue);

  // get envelope-of-envelope:
  Vec metaEnvTime, metaEnvValue;
  getMetaEnvelope(&envTime[0], &envValue[0], (int)envTime.size(), metaEnvTime, metaEnvValue, N);

  // get envelope-of-envelope as audio signal by interpolation:
  Vec t(N);
  rsArray::fillWithRangeLinear(&t[0], N, 0.0, N-1.0);
  interpolateEnvelope(&metaEnvTime[0], &metaEnvValue[0], (int)metaEnvTime.size(),
    &t[0], env, (int)t.size());


  // -maybe the bump can be avoided using a quartic interpolant
  // -and/or: let the env start at 0 and use a segement of lower order by not prescribing values
  //  for the derivative(s) at 0, same at the end, i.e. a quadratic
  // -maybe apply heuristics, when to start and/or end at zero - maybe if the time-delta between
  //  0 and first peak is below k*(time-delta between 1st and 2n peak) where k is some number < 1
  //  do not go to zero

  // smoothing (test):
  //rsBiDirectionalFilter::applyLowpass(env, env, N, 20.0, 44100.0, 4);
  // maybe instead of a filter, use an attack/release slew-rate limiter with zero attack in order
  // to pass through the actual peaks
}

template<class T>
void rsEnvelopeExtractor<T>::getMetaEnvelope(
  const T* rawEnvTime, const T* rawEnvValue, int rawEnvLength,
  std::vector<T>& metaEnvTime, std::vector<T>& metaEnvValue, T endTime)
{
  getPeaks(rawEnvTime, rawEnvValue, rawEnvLength, metaEnvTime, metaEnvValue);
  //rsAssert(rsArray::isSortedStrictlyAscending(&metaEnvTime[0], (int)metaEnvTime.size()));

  //GNUPlotter plt;
  //plt.addDataArrays((int) metaEnvTime.size(), &metaEnvTime[0], &metaEnvValue[0]);
  ////rsPlotVectorsXY(metaEnvTime, metaEnvValue); // debug

  T maxSpacing =   // this must be computed *before* calling setupEndValues!
    maxSpacingMultiplier * rsArray::maxDifference(&metaEnvTime[0], (int)metaEnvTime.size());


  setupEndValues(metaEnvTime, metaEnvValue, endTime);
  //rsAssert(rsArray::isSortedStrictlyAscending(&metaEnvTime[0], (int)metaEnvTime.size()));

  fillSparseAreas(rawEnvTime, rawEnvValue, rawEnvLength, metaEnvTime, metaEnvValue, maxSpacing);
  //rsAssert(rsArray::isSortedStrictlyAscending(&metaEnvTime[0], (int)metaEnvTime.size()));

  ////rsPlotVectorsXY(metaEnvTime, metaEnvValue); // debug
  //metaEnvValue = metaEnvValue; // little offset for visibility
  //plt.addDataArrays((int) metaEnvTime.size(), &metaEnvTime[0], &metaEnvValue[0]);
  //plt.plot();
}
// -maybe fillSparseAreas should be done before setupEndValues?
// -maybe, if there are less than 2 peaks, we should conclude that there is no beating present and
//  skip the de-beating process

template<class T>
void rsEnvelopeExtractor<T>::interpolateEnvelope(const T* envTimes, T* envValues, int envLength,
  const T* interpolatedTimes, T* interpolatedValues, int interpolatedLength)
{
  rsInterpolatingFunction<T, T> intFunc;
  intFunc.setMode(interpolationMode);
  intFunc.interpolate(envTimes, envValues, envLength,
    interpolatedTimes, interpolatedValues, interpolatedLength);
}

template<class T>
void rsEnvelopeExtractor<T>::connectPeaks(const T* envTimes, T* envValues, T* peakValues,
  int length)
{
  rsAssert(rsArray::isSortedStrictlyAscending(&envTimes[0], length));
  std::vector<T> metaEnvTime, metaEnvValue;
  getMetaEnvelope(envTimes, envValues, length, metaEnvTime, metaEnvValue, envTimes[length-1]);
  interpolateEnvelope(&metaEnvTime[0], &metaEnvValue[0], (int)metaEnvTime.size(),
    envTimes, peakValues, length);
  rsAssert(rsLast(metaEnvTime) == envTimes[length-1]);

  //GNUPlotter plt;
  //plt.addDataArrays(length, envTimes, envValues);
  //plt.addDataArrays((int) metaEnvTime.size(), &metaEnvTime[0], &metaEnvValue[0]);
  //plt.plot();
  ////rsPlotVectorsXY(metaEnvTime, metaEnvValue); // debug
  ////rsPlotArraysXY(length, envTimes, envValues); // debug
}

template<class T>
void rsEnvelopeExtractor<T>::setupEndValues(
  std::vector<T>& envTimes, std::vector<T>& envValues, T endTime)
{
  if(envTimes.size() < 2) {
    rsPrepend(envValues, T(0));
    rsPrepend(envTimes,  T(0));
    rsAppend( envValues, T(0));
    rsAppend( envTimes,  endTime);
    return;
  }

  T v;
  if(startMode == ZERO_END) {
    rsPrepend(envValues, T(0));
    rsPrepend(envTimes,  T(0));
  } else if(startMode == EXTRAPOLATE_END) {
    v = rsInterpolateLinear(envTimes[0], envTimes[1], envValues[0], envValues[1], 0.0);
    rsPrepend(envValues, rsMax(v, T(0)));
    rsPrepend(envTimes,  T(0));
  }

  if(endMode == ZERO_END) {
    rsAppend(envValues, T(0));
    rsAppend(envTimes,  endTime);
  } else if(endMode == EXTRAPOLATE_END) {
    int M = (int)envTimes.size()-1;
    v = rsInterpolateLinear(envTimes[M-1], envTimes[M], envValues[M-1], envValues[M], endTime);
    rsAppend(envValues, rsMax(v, T(0)));
    rsAppend(envTimes, endTime);
  }

  rsAssert(rsArray::isSortedStrictlyAscending(&envTimes[0], (int)envTimes.size()));
  // rsAssert(rsLast(envTimes) == endTime); // no - in free-end mode, it may be different
}

template<class T>
void rsEnvelopeExtractor<T>::fillSparseAreas(const T* rawEnvTime, const T* rawEnvValue, int rawEnvLength,
  std::vector<T>& metaEnvTime, std::vector<T>& metaEnvValue, T maxSpacing)
{
  if(maxSpacing == T(0))
  {
    rsCopyToVector(rawEnvTime,  rawEnvLength, metaEnvTime);
    rsCopyToVector(rawEnvValue, rawEnvLength, metaEnvValue);
    return;
  }

  /*
  // for debug:
  GNUPlotter plt1;
  plt1.addDataArrays(rawEnvLength, rawEnvTime, rawEnvValue);
  plt1.addDataArrays((int)metaEnvTime.size(), &metaEnvTime[0], &metaEnvValue[0]);
  plt1.plot();
  */

  std::vector<T> tmpTime, tmpValue;   // buffers for extra datapoints to be inserted
  for(size_t i = 1; i < metaEnvTime.size(); i++) {
    T t1 = metaEnvTime[i];
    T t0 = metaEnvTime[i-1];
    T dt = t1 - t0;
    if(dt > maxSpacing) { // we need to insert extra datapoints between i-1 and i

      int numExtraPoints = (int) floor(dt/maxSpacing);  // verify floor ...maybe use ceil?
      rsAssert(numExtraPoints >= 0);

      // it seems, maxSpacing is zero? :-O

      tmpTime.resize(numExtraPoints);
      tmpValue.resize(numExtraPoints);
      for(int j = 0; j < numExtraPoints; j++) {
        T t = t0 + (j+1) * (dt/(numExtraPoints+1));  // verify this formula
        int idx = rsArray::splitIndexClosest(rawEnvTime, rawEnvLength, t);
        tmpTime[j]  = rawEnvTime[idx];
        tmpValue[j] = rawEnvValue[idx];
      }
      rsInsert(metaEnvTime,  tmpTime,  i);
      rsInsert(metaEnvValue, tmpValue, i);
    }
  }

  /*
  // for debug:
  GNUPlotter plt2;
  plt2.addDataArrays(rawEnvLength, rawEnvTime, rawEnvValue);
  plt2.addDataArrays((int)metaEnvTime.size(), &metaEnvTime[0], &metaEnvValue[0]);
  plt2.plot();
  int dummy = 0;
  */
}

template<class T>
std::vector<size_t> rsEnvelopeExtractor<T>::findPeakIndices(const T* x, int N,
  bool includeFirst, bool includeLast)
{
  std::vector<size_t> peaks;

  bool includePlateaus = true; // todo: make user option

  if(N == 0)
    return peaks;

  if(includeFirst) {
    if(N > 1) {
      if(x[0] >= x[1])      // todo: switch between >= and > depending on includePlateaus
        peaks.push_back(0); }
    else
      peaks.push_back(0);  } // the one and only element is a "peak"

  if(includePlateaus) {
    for(int n = 1; n < N-1; n++) {
      if(RAPT::rsArray::isPeakOrPlateau(x, n))
        peaks.push_back(n); }}
  else {
    for(int n = 1; n < N-1; n++) {
      if(RAPT::rsArray::isPeak(x, n))
        peaks.push_back(n); }}

  if(includeLast){
    if(N > 1) {
      if(x[N-1] >= x[N-2])   // todo: switch between >= and > depending on includePlateaus
        peaks.push_back(N-1); }}

  return peaks;

  // make a unit test for this covering all special cases
}
// what about situations where there are several values of the same height, like
// 0 1 2 2 2 1 0 1 3 3 2 1 2 1
//       *          *      *    desired peaks
// ...actually, it seems like quadratic interpolation using 3 points seems unfair/asymmetric
// - probably it's better to use cubic interpolation and use two points to the left and two
// to the right to find the actual peak location
// or use >= as condition, then we would get
// 0 1 2 2 2 1 0 1 3 3 2 1 2 1
//     * * *       * *     *
// seems better for envelope extraction

template<class T>
void rsEnvelopeExtractor<T>::getAmpEnvelope(const T* x, int N,
  std::vector<T>& sampleTime, std::vector<T>& envValue)
{
  std::vector<T> xAbs(N);
  rsArray::applyFunction(&x[0], &xAbs[0], N, fabs);                  // absolute value
  std::vector<size_t> peakIndices = findPeakIndices(&xAbs[0], N, true, true); // peak indices

  // peak coordinates:
  size_t M = peakIndices.size();
  sampleTime.resize(M);
  envValue.resize(M);
  for(size_t m = 0; m < M; m++) {
    size_t n = peakIndices[m];
    sampleTime[m] = (T) n;
    envValue[m]   = xAbs[n];
    // todo: optionally refine to subsample-precision
  }
}

template<class T>
void rsEnvelopeExtractor<T>::getPeaks(const T *x, const T *y, int N,
  std::vector<T>& peaksX, std::vector<T>& peaksY)
{
  std::vector<size_t> peakIndices = findPeakIndices(y, N, false, false);
    // false because, we don't want to include the end-values, because they will be set
    // setupEndValues in getMetaEnvelope

  size_t M = peakIndices.size();
  peaksX.resize(M);
  peaksY.resize(M);
  for(size_t m = 0; m < M; m++) {
    size_t n = peakIndices[m];
    peaksX[m] = x[n];
    peaksY[m] = y[n];
    // todo: refine to subsample-precision
  }
}

//=================================================================================================

template<class T>
T rsExponentialEnvelopeMatcher<T>::getMatchOffset(const T* x1, int N1, const T* x2, int N2)
{
  // take decibel values of both input envelopes - this transforms the exponential decay into a 
  // linear one and when both decays are the same, the lines have the same slope:
  std::vector<T> xdB1, xdB2;  // envelopes converted to decibels
  std::vector<T> t1,   t2;    // time axes
  xdB1.reserve(N1);
  xdB2.reserve(N2);
  t1.reserve(N1);
  t2.reserve(N2);

  for(int n = initIgnore1; n < N1 - finalIgnore1; n++) {
    if(x1[n] >= ignoreThresh1) {
      t1.push_back(T(n));
      xdB1.push_back(rsAmpToDb(x1[n]));
    }
  }
  for(int n = initIgnore2; n < N2 - finalIgnore2; n++) {
    if(x2[n] >= ignoreThresh2) {
      t2.push_back(T(n));
      xdB2.push_back(rsAmpToDb(x2[n]));
    }
  }

  // find the regression lines y = a*t + b for both signals:
  T a1, b1, a2, b2;                // linear regression coeffs
  rsStatistics::linearRegression((int)t1.size(), &t1[0], &xdB1[0], a1, b1);
  rsStatistics::linearRegression((int)t2.size(), &t2[0], &xdB2[0], a2, b2);

  // compute, by how much we must shift x2 to match x1 at the given match level:
  T tm1, tm2;
  tm1 = (matchLevel - b1) / a1;  // time instant, where xdB1 crosses the matchLevel
  tm2 = (matchLevel - b2) / a2;  // same for xdB2
  return tm1 - tm2;              // this is the resulting desired shift
}
// factor out a function getRegressionCoeffs - can be used for plotting the regression lines

// idea: maybe to be more flexible with respect to the shape of the envelope, we should try 
// polynomial regression: https://en.wikipedia.org/wiki/Polynomial_regression instead of linear
// regression. that would amount to not necessarily assume an exponential envelope shape - instead,
// the shape would be estimated as well. this would complicate the situation in two ways - 1st, 
// obviously, we need an algorithm to find the polynomial coeffs, 2nd, we would have to find the 
// points, where the polynomials are equal to one another in the same way that we now find the 
// point where the lines cross - but that could give multiple solutions - how do we pick the right 
// one? ...or maybe we should not look for points where the curves cross but rather make them match
// in some sort of (to be suitably defined) least-squares sense? ...maybe the user could specify a
// "degree" parameter for the polynomial - with 1, we get the linear regression that we have now

// or maybe spline-regression:  https://data.princeton.edu/eco572/smoothing.pdf ...but for that, we
// would also need a procedure to place the knots

// if we allow for shapes other than exponential, the class should be renamed

//=================================================================================================

template<class T>
void rsSineQuadraturePart(T *x, T *y, int N, T f, T fs, bool backward)
{
  rsOnePoleFilter<T,T> flt;
  flt.setSampleRate(fs);
  flt.setCutoff(f);
  flt.setMode(rsOnePoleFilter<T,T>::ALLPASS_BLT);
  if( backward == true )
  {
    for(int n = N-1; n >= 0; n--)
      y[n] = flt.getSample(x[n]);
  }
  else
  {
    for(int n = 0; n < N; n++)
      y[n] = flt.getSample(x[n]);
  }
}

// used internally by rsSineEnvelopeViaQuadrature, rsSineEnvelopeViaAmpFormula
template<class T>
void rsSmoothSineEnvelope(T *y, int N, T f, T fs, T s)
{
  T a[3], b[3]; // filter coeffs
  if( s > 0.0 && f/s < 0.5*fs )
  {
    rsBiquadDesigner::calculateFirstOrderLowpassCoeffs(
      b[0], b[1], b[2], a[1], a[2], T(1.0/fs), f/s);
    rsArray::negate(a, a, 3);
    a[0] = 1.0;
    rsArray::filter(y, N, y, N, b, 1, a, 1);
    rsArray::reverse(y, N);
    rsArray::filter(y, N, y, N, b, 1, a, 1);
    rsArray::reverse(y, N);
     // do this with a 1st order filter as well, maybe use rsBiDirectionalFilter (lowpass version)
  }
}

template<class T>
void rsSineEnvelopeViaQuadrature(T *x, T *y, int N, T f, T fs, T s)
{
  // get 90° phase shifted version:
  rsSineQuadraturePart(x, y, N, f, fs, true);

  // obtain magnitude of complex phasor:
  for(int n = 0; n < N; n++)
    y[n] = sqrt(x[n]*x[n] + y[n]*y[n]);

  // optional smoothing:
  rsSmoothSineEnvelope(y, N, f, fs, s);
}

template<class T>
void rsSineEnvelopeViaAmpFormula(T *x, T *y, int N, T f, T fs, T s)
{
  T w = 2*PI*f/fs;
  T a, p;
  for(int n = 0; n < N-1; n++)
  {
    rsSineAmplitudeAndPhase(x[n], x[n+1], w, &a, &p);
    y[n] = a;
  }

  // for the last value, we look backwards:
  rsSineAmplitudeAndPhase(x[N-1], x[N-2], w, &a, &p);
  y[N-1] = a;

  // optional smoothing:
  rsSmoothSineEnvelope(y, N, f, fs, s);
}

template<class T>
void rsEnvelopedSine(T *y, int N, T f, T fs, T p, T *a)
{
  rsSineIterator<T> sinIt(2*PI*f/fs, p, 1.0);
  for(int n = 0; n < N; n++)
    y[n] = a[n] * sinIt.getValue();
}

template<class T>
double rsUnwrappedPhaseForCatchSweep(T p0, T pk, T wk, int k, int sweepDirection = 0)
{
  // compute phase that would be reached at k, if we were running a fixed frequency sine at w from
  // 0 to k:
  T pkw = p0 + k*wk;

  // find the unwrapped target phase pk at k - this is given by pk plus/minus some integer multiple
  // of 2*PI, such that we get most closely to pkw:
  T pku = pk;
  int n = 0;
  while( pku < pkw )
  {
    pku = pk + n*2*PI;
    n++;
  }
  if( sweepDirection > 0 || (sweepDirection == 0 && pku-pkw > PI) )
    pku -= 2*PI;

  return pku;
}

template<class T>
void rsEnvelopedPhaseCatchSweep(T *y, int k, T p0, T pk, T wk, T *a, int sweepDirection)
{
  // Notation: p0, pk, w0, wk: instantaneous phases and normalized radian frequencies at sample 0
  // and sample k, respectively, wn: inst. freq. at sample n = 0...k-1. We have
  // pk = p0 + sum_{n=0}^{k-1} wn and we assume a linear frequency sweep such that
  // wn = w0 + n * dw, where dw = (wk-w0)/k. Inserting this into the pk formula, we get:
  // pk = p0 + sum_{n=0}^{k-1} (w0 + n*dw) = p0 + k*w0 + sn*dw, where sn = sum_{n=0}^{k-1} n
  // This can be solved for w0 = (pk - p0 - sn*wk/k) / (k - sn/k)

  pk = rsUnwrappedPhaseForCatchSweep(p0, pk, wk, k, sweepDirection);
  T snk = (T) rsSum(0, k-1) / (T) k;
  T w0  = (pk-p0-snk*wk)/(k-snk); // normalized radian frequency at start
  T dw  = (wk-w0)/k;              // increment for normalized radian frequency
  T wn  = w0;                     // instantaneous normalized radian frequency
  T pn  = p0;                     // instantaneous phase
  for(int n = 0; n < k; n++)
  {
    y[n] = a[n] * sin(pn);
    wn   = w0 + n*dw;
    pn  += wn;
  }

  // todo: factor out a function to compute an array of w[n] - mainly for plotting the frequency
  // sweep and comparison between linear and parabolic sweep (to be implemented)
}

template<class T>
void rsPhaseCatchSweep(T *y, int k, T p0, T pk, T wk, int sweepDirection)
{
  rsArray::fillWithValue(y, k, 1.0);
  rsEnvelopedPhaseCatchSweep(y, k, p0, pk, wk, y, sweepDirection);
}

template<class T>
void rsEnvelopedPhaseCatchSine(T *y, int N, T f, T fs, T p0, T pk, int k, T *a, int sweepDirection)
{
  // 1st part (sweeping frequency):
  T w = 2*PI*f/fs;
  rsEnvelopedPhaseCatchSweep(y, k, p0, pk, w, a, sweepDirection);

  // 2nd part (fixed frequency):
  rsEnvelopedSine(&y[k], N-k, f, fs, pk, &a[k]);
}

template<class T>
void rsRecreateSine(T *x, T *y, int N, T fx, T fy, T fs, T p0, T smooth)
{
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);
  rsEnvelopedSine(y, N, fy, fs, p0, y);
}

template<class T>
void rsRecreateSineWithPhaseCatch(T *x, T *y, int N, T fx, T fy, T fs, T p0, T pk, int k,
  T smooth, int sd)
{
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);
  rsEnvelopedPhaseCatchSine(y, N, fy, fs, p0, pk, k, y, sd);
}

/*
template<class T>
void rsRecreateSine(T *x, T *y, int N, T fx, T fy, T fs, T p0, T smooth)
{
  // get the input sinusoid's envelope (we may use y as temporary buffer for this)
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);

  // generate the sine at the desired frequency and multiply it by the envelope:
  T w = 2*PI*fy/fs;
  T p = p0;
  for(int n = 0; n < N; n++)
  {
    y[n] *= sin(p);
    p    += w;
  }

  // optimize using rsSineOscillator

  // factor out the sine creation part - should take an amplitude array as input
}
*/

template<class T>
T getMaxShortTimeRMS(T* x, int N, int averagingLength)
{
  rsMovingAverage<T, T> ma;
  ma.setLengthInSamples(averagingLength);
  T maxRms = T(0);
  for(int n = 0; n < N; n++) {
    T rms = sqrt( ma.getSample(x[n]*x[n]) );
    if(rms > maxRms)
      maxRms = rms;
  }
  return maxRms;
}


template<class T>
T extremumViaLineIntersect(const T* x, const int N, const int k)  
{
  // make a similar function extremumViaParabola, move to rsArray

  // check, if our minimum point has two neighbours to each side:
  if(k < 2 || k > N-3) return k;  // can we do something better?

  // use the two neighbours of the minimum to the left and the two neighbours to the right to 
  // define two lines and find the x-coordinate of the intersection of these two lines - this is 
  // our estimate with subsample precision:
  T al, bl, ar, br;
  rsLine2D<T>::twoPointToExplicit(T(-1), x[k-1], T(-2), x[k-2], &al, &bl);
  rsLine2D<T>::twoPointToExplicit(T(+1), x[k+1], T(+2), x[k+2], &ar, &br);
  // this can be optimized - it divides by x2-x1 = 1 -> divisions are superfluous

  const T d = (br-bl) / (al-ar);  // delta
  return T(k) + rsClip(d, T(-1), T(+1));
}
// maybe a parabolic fit would be better - and maybe squared differences are better for this?
// perhaps, even an FFT based algo can be used for the squared differences (messing around with
// the two autocorrelations and the cross-correlation?)

template<class T>
T rsEnvelopeMatchOffset(const T* x, const int Nx, const T* y, const int Ny)
{
  std::vector<T> s(Nx);
  for(int k = 0; k < Nx; k++)
    s[k] = rsArray::meanOfAbsoluteDifferences(&x[k], &y[0], rsMin(Nx-k, Ny));

  //rsPlotVector(s);

  const int k = RAPT::rsArray::minIndex(&s[0], Nx); // index of minimum...
  return extremumViaLineIntersect(&s[0], Nx, k);    // ...refined to subsample precision
}

template<class T>
T rsEnvelopeMatchOffset(const T* x, const int Nx, const T* y, const int Ny, const int D)
{
  if(D == 1)
    return rsEnvelopeMatchOffset(x, Nx, y, Ny);
  else {
    const int NxD = Nx/D;
    const int NyD = Ny/D;
    std::vector<T> xd(NxD), yd(NyD);

    RAPT::rsArray::decimateViaMean(&x[0], Nx, &xd[0], D); 
    RAPT::rsArray::decimateViaMean(&y[0], Ny, &yd[0], D);

    //RAPT::rsArray::decimate(&x[0], Nx, &xd[0], D);  // use decimateViaMean
    //RAPT::rsArray::decimate(&y[0], Ny, &yd[0], D);

    //// debug:
    //T dbg = rsEnvelopeMatchOffset(&xd[0], NxD, &yd[0], NyD);
    //dbg *= D;
    //return dbg;

    // todo: maybe use log-of-envelope instead of raw envelope - maybe use a boolean parameter
    // or dB envelopes - and use a threshold, like: 
    // sum += abs(max(e1[n], thres) - max(e2[n], thresh))
    // or better: apply that threshold when computing the dB values
    // ...or maybe it's about time to wrap it all into a class - let user switch between lin/log
    // and absolute/squared difference, and parabolic/linear interpolation for subsample precision

    return T(D) * rsEnvelopeMatchOffset(&xd[0], NxD, &yd[0], NyD);
  }
}
// maybe merge this with the rsExponentialEnvelopeMatcher class (rename to rsEnvelopeMatcher), so we
// may also use its ingnore-facilities - needs functions 
//  setAlgorithm(absoluteDifference, squaredDifference, correlation, linearRegression, ...)
//  setDecimation, setInterpolation






template<class T>
void rsExpDecayParameters(T t1, T a1, T t2, T a2, T* A, T* tau)
{
  T dt =  t2 - t1;               // time difference
  T ra =  a2 / a1;               // amplitude ratio
  *tau = -dt / log(ra);          // time-constant of exponential decay
  *A   =  a1 / exp(-t1 / *tau);  // amplitude multiplier
}
// can we get rid of the call to exp? 

template<class T>
std::vector<T> rsExpDecayTail(int numFrames, const T* timeArray, const T* ampArray, 
  int matchIndex1, int matchIndex2, T sampleRate, T freq, T phase, int phaseMatchIndex, 
  int numSamples)
{
  rsAssert(matchIndex2 < numFrames);
  rsAssert(matchIndex1 < matchIndex2);

  // estimate exponential decay parameters:
  T A, tau;
  rsExpDecayParameters(timeArray[matchIndex1], ampArray[matchIndex1], 
    timeArray[matchIndex2], ampArray[matchIndex2], &A, &tau);

  // generate exponentially enveloped sinusoid:
  if(numSamples == -1) // use default length when user passes no desired length
    numSamples = (int) ceil(timeArray[numFrames-1] * sampleRate); 
  std::vector<T> x(numSamples);
  T ts = timeArray[phaseMatchIndex];  // time-instant for splicing
  T p0 = phase - 2*PI*freq*ts;        // start-phase
  T w  = 2*PI*freq/sampleRate;        // normalized radian frequency
  for(int n = 0; n < numSamples; n++)
    x[n]  = A * exp(-(n/sampleRate) / tau) * cos(w*n + p0);
  // we use the convention here the the phase is with respect to a cosine wave - this is consistent
  // with the sinusoidal modeling framework....but the modal synthesis stuff uses a sine...hmmm...
  // this should probably be treated consistently...

  // todo: optimize: use the exponential-decay filter instead of calling exp/cos explicitly

  return x;
}
