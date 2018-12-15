// Construction/Destruction:

template<class T>
rsSpectrogram<T>::rsSpectrogram()
{
  //init();
  setBlockSize(512);
}

template<class T>
rsSpectrogram<T>::~rsSpectrogram()
{

}

// Setup:


// Inquiry:

template<class T>
int rsSpectrogram<T>::getNumFrames(int N, int H)
{
  int tmp = N+H-1;
  int F = tmp/H;
  if(tmp%H > 0)
    F += 1;
  return F; // F == ceil((N+H-1)/H)
}

template<class T>
T rsSpectrogram<T>::getWindowSum(T *wa, T *ws, int B, int H)
{
  T s = 0;
  for(int n = 0; n < B; n += H)
    s += wa[n] * ws[n];
  return s;
}

// Processing:


template<class T>
rsMatrix<std::complex<T>> rsSpectrogram<T>::complexSpectrogram(const T *signal, int numSamples) const
{
  // preliminary - call static function - later, we should use an object of type 
  // rsFourierTransformerBluestein for more efficiency (it does some pre-computations)
  return complexSpectrogram(
    signal, numSamples, &analysisWindow[0], blockSize, hopSize, zeroPaddingFactor);
}

/*
template<class T>
void rsSpectrogram<T>::hanningWindowZN(T *w, int N)
{
  T s = 2*PI/N; // for a window that ends with zero: w[N-1]=0, this would be s=2*PI/(N-1)
  for(int n = 0; n < N; n++)
    w[n] = 0.5*(1-cos(s*n));
}
*/

// x: signal, N: number of samples, n: block center sample, w: window, B: blocksize, M: FFT size,
// X: complex short-time spectrum (output)
template<class T>
void rsSpectrogram<T>::shortTimeSpectrum(const T* x, int N, int n, const T* w,
  int B, int M, std::complex<T> *X)
{
  int pad = (M-B)/2;                        // amount of pre/post zero padding
  if(pad > 0)
  {
    rsArray::fillWithZeros(X, pad);         // pre padding
    rsArray::fillWithZeros(&X[M-pad], pad); // post padding
  }
  std::complex<T> *Xs = &X[pad];            // pointer, from which we write into the X-array
  rsArray::copySection(x, N, Xs, n-B/2, B); // copy signal section into FFT buffer (convert to complex)
  for(int i = 0; i < B; i++)                // apply window
    Xs[i] *= w[i];

  // shift buffer to make the window center the reference for phase values
  // ...

  rsFFT(X, M);                        // transform to frequency domain
}

template<class T>
rsMatrix<std::complex<T>> rsSpectrogram<T>::complexSpectrogram(const T* x, int N, const T* w,
  int B, int H, int P)
{
  // x: signal, N: number of samples, w: window, B: blocksize, H: hopsize, P: padding factor

  int F = getNumFrames(N, H);                  // number of STFT frames
  int M = B * P;                               // FFT size (maybe use L)
  int K = M/2 + 1;                             // number of non-redundant bins
  rsMatrix<std::complex<T>> s(F, K);           // spectrogram (only positive frequency bins)
  T a = 2 / rsArray::sum(w, B);                // amplitude scaler
  int n = 0;                                   // sample, where current block is centered
  std::complex<T> *X = new std::complex<T>[M]; // short-time spectrum centered at sample n
  for(int i = 0; i < F; i++)                   // loop over frames
  {
    shortTimeSpectrum(x, N, n, w, B, M, X);    // obtain STFT centered at n
    for(int j = 0; j < K; j++)                 // collect results for positive frequencies
      s(i, j) = a * X[j];
    n += H;                                    // advance n by the hop size
  }
  delete[] X;
  return s;
}

template<class T>
std::vector<T> rsSpectrogram<T>::synthesize(const rsMatrix<std::complex<T>> &s, T *ws,
  int B, int H, T *wa)
{
  // s: spectrogram, ws: synthesis-window, B: block size, H: hop size, wa: analysis window,
  std::vector<T> y = synthesizeRaw(s, ws, B, H);
  std::vector<T> m = getModulation(wa, ws, B, H, s.getNumRows());
  T a = rsArray::sum(wa, B) / 2;
  for(unsigned int n = 0; n < y.size(); n++)
    y[n] *= (a/m[n]);
  return y;
}

template<class T>
rsMatrix<T> rsSpectrogram<T>::timeReassignment(T *x, int N,
  const rsMatrix<std::complex<T>> &s, T *wr, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wr: time-ramped window, B: blocksize,
  // H: hopsize

  rsMatrix<T> tr;

  // use the complexSpectrogram function to compute a spectrogram with the ramped window and then
  // apply the time reassignment formula using the original spectrogram s and the "ramped"
  // spectrogram to compute corresponding value of the time reassignment matrix tr
  // ...

  return tr;
}

template<class T>
rsMatrix<T> rsSpectrogram<T>::frequencyReassignment(T *x, int N,
  const rsMatrix<std::complex<T>> &s, T *wd, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wd: window derivative, B: blocksize,
  // H: hopsize

  rsMatrix<T> fr;

  // use the complexSpectrogram function to compute a spectrogram with the ramped window and then
  // apply the frequency reassignment formula using the original spectrogram s and the "ramped"
  // spectrogram to compute corresponding value of the frequency reassignment matrix fr
  // ...

  return fr;
}

/*
// add length-L array y into length-N array x starting at n
template<class T>
void addInto(T *x, int N, T *y, int L, int n = 0)
{
  int r = 0;                // read start
  if(n < 0)
  {
    L += n;
    r -= n;
    n  = 0;
  }
  int d = n + L - N;        // number of overhanging values
  if(d > 0)
    L -= d;
  for(int i = 0; i < L; i++)
    x[n+i] += y[r+i];
}
// the index manipulation code can be factored out
*/

template<class T>
std::vector<T> rsSpectrogram<T>::synthesizeRaw(const rsMatrix<std::complex<T>> &s,
  T *w, int B, int H)
{
  // w: window, B: blocksize, H: hopsize, s: complex spectrogram

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
    for(k = 1; k < K; k++)
    {
      Y[k]   = s(i, k);
      Y[M-k] = conj(Y[k]);
    }
    rsIFFT(Y, M);

    // apply synthesis-window and overlap/add into output signal:
    for(k = 0; k < B; k++)
      g[k] = Y[k0+k].real() * w[k];  // k0 != 0, if zero-padding was used
    rsArray::addInto(y.data(), N, g.data(), B, i*H-B/2);
  }

  delete[] Y;
  return y;
}

template<class T>
std::vector<T> rsSpectrogram<T>::getModulation(T *wa, T *ws, int B, int H, int F)
{
  // wa: analysis window, ws: synthesis-window, B: block size, H: hop size, F: number of frames
  int N = (F-1) * H + B/2;      // number of samples
  std::vector<T> y(N);          // modulation signal
  T *w = new T[B];              // product-window
  for(int n = 0; n < B; n++)
    w[n] = wa[n] * ws[n];
  for(int i = 0; i < F; i++)
    rsArray::addInto(y.data(), N, w, B, i*H-B/2);
  delete[] w;
  return y;
}

// Misc:
/*
template<class T>
void rsSpectrogram<T>::init()
{
  //fs = 44100.0;   // samplerate
  //Nb = 512;       // block size
  //Nh = Nb/2;      // hop size
  //Nf = Nb*2;      // FFT size
}
*/

template<class T>
void rsSpectrogram<T>::updateAnalysisWindow()
{
  analysisWindow.resize(blockSize); // later: analysisBlockSize
  fillWindowArray(&analysisWindow[0], blockSize, analysisWindowType);
  // todo: create also the time-derivative and the time-ramped window for time/frequency 
  // reassignment later
}

template<class T>
void rsSpectrogram<T>::updateSynthesisWindow()
{
  synthesisWindow.resize(blockSize); // later: synthesisBlockSize
  fillWindowArray(&synthesisWindow[0], blockSize, synthesisWindowType);
}

template<class T>
void rsSpectrogram<T>::fillWindowArray(T* w, int length, int type)
{
  rsWindowFunction::createWindow(w, length, type);
  // ...actually, this function is obsolete now...but maybe not when we later create derivative and
  // ramp windows for time/frequency reassignment
}



/*
todo:

-the original phase vocoder relies on multiplication of the signal with various complex
 sinusoids ("heterodyning") and a filter bank - maybe additionally implement a true phase vocoder
 (see DAFX, Ch. 8 - Time-frequency processing)
-maybe that true phase vocoder may be extended by allowing the complex sinusoids to be not 
 necessarily harmonic and vary in frequency (and maybe also allowing the filters to have 
 time-variying bandwidths) - this could be used to track the frequency of the partial being 
 analyzed (but requires a preliminary knowlegde of the frequency trajectory of the partial=
 ...maybe we can firsat obtain preliminary (rough) frequency tracks by spectrogram processing and
 the refine them by a time-varying phase vocoder approach...and then use that data to refine 
 further etc. until it converges

*/