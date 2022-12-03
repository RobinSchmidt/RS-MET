// Construction/Destruction:

template<class T>
rsSpectrogramProcessor<T>::rsSpectrogramProcessor()
{
  //transformer.setNormalizationMode(rsFourierTransformerRadix2<T>::NEVER_NORMALIZE);
  // this is, how it should work

  transformer.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO);
  // and this is how it actually works - why do we need the normalizaion on iFFT? ...see comment in 
  // synthesize - there seems to be an additional unnormalization - but it's weird that it works 
  // anyway because the normalization just divides by the FFT size whereas this "unnormalization" 
  // multiplies by the sum of window values which is a different number, unless a rectangular 
  // window is used -> weird -> figure out

  setTrafoSize(512);
  setBlockSize(512);
}

template<class T>
rsSpectrogramProcessor<T>::~rsSpectrogramProcessor()
{

}

// Setup:

template<class T>
void rsSpectrogramProcessor<T>::setBlockSize(int newSize) 
{ 
  setBlockAndTrafoSize(newSize, trafoSize);
}

template<class T>
void rsSpectrogramProcessor<T>::setTrafoSize(int newSize) 
{ 
  setBlockAndTrafoSize(blockSize, newSize);
}

template<class T>
void rsSpectrogramProcessor<T>::setBlockAndTrafoSize(int newBlockSize, int newTrafoSize)
{
  rsAssert(newTrafoSize >= newBlockSize, "FFT size should be >= block size");
  if(newBlockSize != blockSize) {
    blockSize = newBlockSize;
    updateAnalysisWindow();
    updateSynthesisWindow();
  }
  if(newTrafoSize != trafoSize) {
    trafoSize = newTrafoSize;
    transformer.setBlockSize(trafoSize);
  }
}

// Inquiry:

template<class T>
int rsSpectrogramProcessor<T>::getNumFrames(int N, int H)
{
  int tmp = N+H-1;
  int F = tmp/H;   // integer division gives floor(tmp/H)
  if(tmp%H > 0)    // if there was a truncated fractional part in the floor-division... 
    F += 1;        // ...add 1 to get the ceil (otherwise floor == ceil)
  return F;        // F == ceil((N+H-1)/H)
}

template<class T>
T rsSpectrogramProcessor<T>::getWindowSum(T *wa, T *ws, int B, int H)
{
  T s = 0;
  for(int n = 0; n < B; n += H)
    s += wa[n] * ws[n];
  return s;
}

template<class T>
std::vector<T> rsSpectrogramProcessor<T>::getRoundTripModulation(int F)
{
  //  F: number of frames

  T*  wa = &analysisWindow[0];
  T*  ws = &synthesisWindow[0];
  int B  = blockSize;
  int H  = hopSize;
  int N = (F-1) * H + B/2;      // number of samples
  std::vector<T> y(N);          // modulation signal
  T *w = new T[B];              // product-window
  for(int n = 0; n < B; n++)
    w[n] = wa[n] * ws[n];
  for(int i = 0; i < F; i++)
    rsArrayTools::addInto(y.data(), N, w, B, i*H-B/2);
  delete[] w;
  return y;
}

// Processing:

// x: signal, N: number of samples, n: block center sample, X: complex short-time spectrum (output)
template<class T>
void rsSpectrogramProcessor<T>::shortTimeSpectrum(const T* x, int N, int n, std::complex<T> *X)
{
  prepareTrafoBuffer(x, N, n, X);
  fft(X, trafoSize);
}

// B: blockSize, M: trafoSize, L: amount of left-padding, R: amount of right-padding
template<class T>
void rsSpectrogramProcessor<T>::getLeftRightPaddingAmount(int B, int M, int* L, int* R)
{
  int P = M-B;        // total amount of padding
  if(rsIsEven(B)) {
    *L = P/2;
    *R = P - *L; }
  else {
    *R = P/2;
    *L = P - *R; }
}

template<class T>
void rsSpectrogramProcessor<T>::prepareTrafoBuffer(const T* x, int N, int n, std::complex<T> *X)
{
  T*  w = &analysisWindow[0];
  int B = blockSize;
  int M = trafoSize;
  int L, R;
  getLeftRightPaddingAmount(B, M, &L, &R);
  rsArrayTools::fillWithZeros( X,      L);       // left zero padding
  rsArrayTools::fillWithZeros(&X[M-R], R);       // right zero padding
  std::complex<T> *Xs = &X[L];              // pointer, from which we write into the X-array
  rsArrayTools::copySection(x, N, Xs, n-B/2, B); // copy signal section into FFT buffer
  for(int i = 0; i < B; i++)                // apply window
    Xs[i] *= w[i];
  swapForZeroPhase(X, M);
/*
// plot the windowed segment:
#ifdef RS_DEBUG_PLOTTING
  std::vector<T> dbg(M);
  for(int i = 0; i < M; i++)
    dbg[i] = real(X[i]);
  GNUPlotter plt;
  plt.plotArrays(M, &dbg[0]);
#endif
*/
}

template<class T>
rsMatrix<std::complex<T>> rsSpectrogramProcessor<T>::getComplexSpectrogram(const T* x, int N)
{
  // x: signal, N: number of samples

  T*  w = &analysisWindow[0];
  //int B = blockSize;
  int H = hopSize;
  int F = getNumFrames(N, H);                  // number of STFT frames
  int M = trafoSize;
  int K = M/2 + 1;                             // number of non-redundant bins
  rsMatrix<std::complex<T>> s(F, K);        // spectrogram (only positive frequency bins)
  //T a = 2 / rsArrayTools::sum(w, B);              // amplitude scaler
  T a = getAnalysisScaler();                   // amplitude scaler = 2 / rsArrayTools::sum(w, B)
  int n = 0;                                   // sample, where current block is centered
  std::complex<T> *X = new std::complex<T>[M]; // short-time spectrum centered at sample n
  for(int i = 0; i < F; i++)                   // loop over frames
  {
    // todo: maybe obtain amplitude compensation factors for the first few and last few frames 
    // which are not filled completely with samples leading to amplitude estimation errors
    // ...if we do this, the resynthesis should apply inverse factors



    shortTimeSpectrum(x, N, n, X);             // obtain STFT centered at n
    for(int j = 0; j < K; j++)                 // collect results for positive frequencies
      s(i, j) = a * X[j];
    n += H;                                    // advance n by the hop size
  }
  delete[] X;
  return s;
}

template<class T>
void rsSpectrogramProcessor<T>::lowpass(rsMatrix<std::complex<T>>& s, int hi)
{
  int numFrames = s.getNumRows();
  int numBins   = s.getNumColumns();
  for(int i = 0; i < numFrames; i++)
    for(int k = hi+1; k < numBins; k++)
      s(i, k) = T(0);
  // todo: maybe allow for a floating point cutoff value - the amplitude of the last bin would then
  // be something between 0 and 1
}

/*
// still buggy
template<class T>
void rsSpectrogramProcessor<T>::lowpass(rsMatrix<std::complex<T>>& s, T cutoffBin)
{
  int hi = (int) floor(cutoffBin);
  lowpass(s, hi);

  int numBins   = s.getNumColumns();
  if(hi+1 >= numBins) // do we actually need to check that edge case?
    return;
  scaleBin(s, hi + 1, cutoffBin - hi);
}
*/

template<class T>
void rsSpectrogramProcessor<T>::highpass(rsMatrix<std::complex<T>>& s, int lo)
{
  int numFrames = s.getNumRows();
  int numBins   = s.getNumColumns();
  for(int i = 0; i < numFrames; i++)
    for(int k = 0; k < lo; k++)
      s(i, k) = T(0);
}

/*
// still buggy
template<class T>
void rsSpectrogramProcessor<T>::highpass(rsMatrix<std::complex<T>>& s, T cutoffBin)
{
  int lo = (int) ceil(cutoffBin);
  highpass(s, lo);

  int numBins = s.getNumColumns();
  if(lo-1 < 0)  // do we actually need to check that edge case?
    return;
  scaleBin(s, lo, lo - cutoffBin);  // verify formula for scale factor
}
*/

template<class T>
void rsSpectrogramProcessor<T>::bandpass(rsMatrix<std::complex<T>>& s, int lo, int hi)
{
  lowpass( s, hi);
  highpass(s, lo);
}

template<class T>
std::vector<T> rsSpectrogramProcessor<T>::synthesize(const rsMatrix<std::complex<T>> &s)
{
  // s: spectrogram
  int B = blockSize;
  int H = hopSize;
  T*  wa = &analysisWindow[0];
  T*  ws = &synthesisWindow[0];
  std::vector<T> y = synthesizeRaw(s);
  T a = rsArrayTools::sum(wa, B) / 2;  // might this be the additional scaling because of which we
                                  // need to set NORMALIZE_ON_INVERSE_TRAFO?
  if(demodulateOutput == true) {
    std::vector<T> m = getRoundTripModulation(s.getNumRows());
    for(unsigned int n = 0; n < y.size(); n++)
      y[n] *= (a/m[n]);
    //rsPlotVector(m);
  }
  else
    rsArrayTools::scale(&y[0], (int) y.size(), a);
  return y;
}

template<class T>
rsMatrix<T> rsSpectrogramProcessor<T>::timeReassignment(T *x, int N,
  const rsMatrix<std::complex<T>> &s, T *wr, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wr: time-ramped window, B: blocksize,
  // H: hopsize

  rsMatrix<T> tr;

  // use the getComplexSpectrogram function to compute a spectrogram with the ramped window and then
  // apply the time reassignment formula using the original spectrogram s and the "ramped"
  // spectrogram to compute corresponding value of the time reassignment matrix tr
  // ...

  return tr;
}

template<class T>
rsMatrix<T> rsSpectrogramProcessor<T>::frequencyReassignment(T *x, int N,
  const rsMatrix<std::complex<T>> &s, T *wd, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wd: window derivative, B: blocksize,
  // H: hopsize

  rsMatrix<T> fr;

  // use the getComplexSpectrogram function to compute a spectrogram with the derivative window and 
  // then apply the frequency reassignment formula using the original spectrogram s and the 
  // "derivative" spectrogram to compute corresponding value of the frequency reassignment matrix 
  // fr

  return fr;
}

template<class T>
std::vector<T> rsSpectrogramProcessor<T>::synthesizeRaw(const rsMatrix<std::complex<T>> &s)
{
  // s: complex spectrogram

  int B  = blockSize;
  int H  = hopSize;
  T*  w  = &synthesisWindow[0];
  int F  = s.getNumRows();           // number of frames
  int K  = s.getNumColumns();        // number of (non-redundant) bins
  int N  = (F-1) * H + B/2;          // number of samples
  int M  = (K-1) * 2;                // FFT size
  int k0 = (M-B) / 2;                // read start in resynthesized grain before the
                                     // resynthesis window is applied
  std::vector<T> y(N);               // allocate signal
  std::vector<T> g(B);               // grain
  std::complex<T> *Y = new std::complex<T>[M];  // short-time spectrum

  int i, k;
  for(i = 0; i < F; i++)
  {
    // create symmetrized FFT spectrum and transform to time domain:
    Y[0] = s(i, 0);
    for(k = 1; k < K; k++) {
      Y[k]   = s(i, k);
      Y[M-k] = conj(Y[k]);
    }
    ifft(Y, M);
    swapForZeroPhase(Y, M);

    // apply synthesis-window:
    for(k = 0; k < B; k++)
      g[k] = Y[k0+k].real() * w[k];  // k0 != 0, if zero-padding was used

    // overlap/add into output signal:
    rsArrayTools::addInto(y.data(), N, g.data(), B, i*H-B/2);
  }

  delete[] Y;
  return y;
}



// Misc:

template<class T>
void rsSpectrogramProcessor<T>::updateAnalysisWindow()
{
  analysisWindow.resize(blockSize); // later: analysisBlockSize
  fillWindowArray(&analysisWindow[0], blockSize, analysisWindowType);
  // todo: create also the time-derivative and the time-ramped window for time/frequency 
  // reassignment later
}

template<class T>
void rsSpectrogramProcessor<T>::updateSynthesisWindow()
{
  synthesisWindow.resize(blockSize); // later: synthesisBlockSize
  fillWindowArray(&synthesisWindow[0], blockSize, synthesisWindowType);
}

template<class T>
void rsSpectrogramProcessor<T>::fillWindowArray(T* w, int length, rsWindowFunction::WindowType type)
{
  rsWindowFunction::createWindow(w, length, type, false); 
  // ...actually, this function is obsolete now...but maybe not when we later create derivative and
  // ramp windows for time/frequency reassignment
}

template<class T>
void rsSpectrogramProcessor<T>::swapForZeroPhase(std::complex<T>* X, int L)
{
  if(timeOriginAtWindowCenter)
    rsArrayTools::circularShift(X, L, -L/2);  
    // optimize: pass tmpBuffer as workspace to avoid temporary memory allocation
}

template<class T>
void rsSpectrogramProcessor<T>::fft(std::complex<T> *X, int M)
{
  transformer.setDirection(rsFourierTransformerRadix2<T>::FORWARD);
  transformer.setBlockSize(M);
  transformer.transformComplexBufferInPlace(X);

  //rsFFT(X, M);   // old
}

template<class T>
void rsSpectrogramProcessor<T>::ifft(std::complex<T> *X, int M)
{
  transformer.setDirection(rsFourierTransformerRadix2<T>::INVERSE);
  transformer.setBlockSize(M);
  transformer.transformComplexBufferInPlace(X);

  //rsIFFT(X, M);  // old - does scale by 1/M
}


/*
Notes:
-the identity analysis/resynthesis works, but only if the final demodulation step is executed
-when using H = B/4, it almost works even without the demodulation step, when using the 
 HANNING_ZN window for analysis and synthesis
 -we need H = B/4 instead of H = B/2 because the window is applied twice, once in the analysis
  and then again in the synthesis
 -at the start and end of the resynthesized signal, there's a little fade-in/fade-out of length 
  B/4 (i think), because the first window starts n = 0-B/2 - where for perfect overlap/add, it 
  should start at n = 0-3*B/4
 -the overlapped windows add up to 1.5 instead of 1.0
-maybe in practice it's indeed best to just do the final demodulation step but choose windows 
 where this step does not have to do very much (or almost nothing)
-maybe we could prepend and append dummy sections of length B before analysis and cut them off 
 after resynthesis to deal with the fade-in/out issues (that's not very elegant, though)




todo:

-the original phase vocoder relies on multiplication of the signal with various complex
 sinusoids ("heterodyning") and a filter bank - maybe additionally implement a true phase vocoder
 (see DAFX, Ch. 8 - Time-frequency processing)
-maybe that true phase vocoder may be extended by allowing the complex sinusoids to be not 
 necessarily harmonic and vary in frequency (and maybe also allowing the filters to have 
 time-variying bandwidths) - this could be used to track the frequency of the partial being 
 analyzed (but requires a preliminary knowlegde of the frequency trajectory of the partial=
 ...maybe we can first obtain preliminary (rough) frequency tracks by spectrogram processing and
 the refine them by a time-varying phase vocoder approach...and then use that data to refine 
 further etc. until it converges

See also:
https://en.wikipedia.org/wiki/Time%E2%80%93frequency_analysis
https://en.wikipedia.org/wiki/Uncertainty_principle#Signal_processing
https://en.wikipedia.org/wiki/Reassignment_method
https://en.wikipedia.org/wiki/Least-squares_spectral_analysis


*/