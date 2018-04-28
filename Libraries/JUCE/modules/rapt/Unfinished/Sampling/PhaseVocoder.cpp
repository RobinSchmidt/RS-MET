using namespace RSLib;

// Construction/Destruction:

rsPhaseVocoder::rsPhaseVocoder()
{
  init();

  int dummy = 0;
}

rsPhaseVocoder::~rsPhaseVocoder()
{

}

// Setup:


// Inquiry:

int rsPhaseVocoder::getNumFrames(int N, int H)
{
  int tmp = N+H-1;
  int F = tmp/H;
  if(tmp%H > 0)
    F += 1;
  return F; // F == ceil((N+H-1)/H)
}

double rsPhaseVocoder::getWindowSum(double *wa, double *ws, int B, int H)
{
  double s = 0;
  for(int n = 0; n < B; n += H)
    s += wa[n] * ws[n];
  return s;
}

// Processing:

void rsPhaseVocoder::hanningWindowZN(double *w, int N)
{
  double s = 2*PI/N; // for a window that ends with zero: w[N-1]=0, this would be s=2*PI/(N-1)
  for(int n = 0; n < N; n++)
    w[n] = 0.5*(1-cos(s*n));
}

// x: signal, N: number of samples, n: block center sample, w: window, B: blocksize, M: FFT size, 
// X: complex short-time spectrum (output)
void rsPhaseVocoder::shortTimeSpectrum(double *x, int N, int n, double *w, 
  int B, int M, rsComplexDbl *X)
{
  int pad = (M-B)/2;                 // amount of pre/post zero padding
  if(pad > 0)
  {
    rsFillWithZeros(X, pad);         // pre peadding
    rsFillWithZeros(&X[M-pad], pad); // post padding
  }
  rsComplexDbl *Xs = &X[pad];        // pointer, from which we write into the X-array
  rsCopySection(x, N, Xs, n-B/2, B); // copy signal section into FFT buffer (convert to complex)
  for(int i = 0; i < B; i++)         // apply window
    Xs[i] *= w[i];

  // shift buffer to make the window center the reference for phase values
  // ...

  rsFFT(X, M);                        // transform to frequency domain
}

rsMatrix<rsComplexDbl> rsPhaseVocoder::complexSpectrogram(double *x, int N, double *w, 
  int B, int H, int P)
{
  // x: signal, N: number of samples, w: window, B: blocksize, H: hopsize, P: padding factor

  int F = getNumFrames(N, H);               // number of STFT frames
  int M = B * P;                            // FFT size (maybe use L)
  int K = M/2 + 1;                          // number of non-redundant bins
  rsMatrix<rsComplexDbl> s(F, K);           // spectrogram (only positive frequency bins)
  double a = 2 / rsSum(w, B);               // amplitude scaler
  int n = 0;                                // sample, where current block is centered
  rsComplexDbl *X = new rsComplexDbl[M];    // short-time spectrum centered at sample n
  for(int i = 0; i < F; i++)                // loop over frames
  {
    shortTimeSpectrum(x, N, n, w, B, M, X); // obtain STFT centered at n
    for(int j = 0; j < K; j++)              // collect results for positive frequencies
      s(i, j) = a * X[j];
    n += H;                                 // advance n by the hop size
  }
  delete[] X;
  return s;
}

std::vector<double> rsPhaseVocoder::synthesize(const rsMatrix<rsComplex<double>> &s, double *ws, 
  int B, int H, double *wa)
{
  // s: spectrogram, ws: synthesis-window, B: block size, H: hop size, wa: analysis window,
  std::vector<double> y = synthesizeRaw(s, ws, B, H);
  std::vector<double> m = getModulation(wa, ws, B, H, s.getNumRows());
  double a = rsSum(wa, B) / 2;
  for(unsigned int n = 0; n < y.size(); n++)
    y[n] *= (a/m[n]);
  return y;
}

rsMatrix<double> rsPhaseVocoder::timeReassignment(double *x, int N,
  const rsMatrix<rsComplex<double>> &s, double *wr, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wr: time-ramped window, B: blocksize, 
  // H: hopsize

  rsMatrix<double> tr;

  // use the complexSpectrogram function to compute a spectrogram with the ramped window and then 
  // apply the time reassignment formula using the original spectrogram s and the "ramped" 
  // spectrogram to compute corresponding value of the time reassignment matrix tr
  // ...

  return tr;
}

rsMatrix<double> rsPhaseVocoder::frequencyReassignment(double *x, int N,
  const rsMatrix<rsComplex<double>> &s, double *wd, int B, int H)
{
  // x: signal, N: number of samples, s: complex spectrogram, wd: window derivative, B: blocksize, 
  // H: hopsize

  rsMatrix<double> fr;

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

std::vector<double> rsPhaseVocoder::synthesizeRaw(const rsMatrix<rsComplex<double>> &s, double *w, 
  int B, int H)
{
  // w: window, B: blocksize, H: hopsize, s: complex spectrogram

  int F  = s.getNumRows();                // number of frames
  int K  = s.getNumColumns();             // number of (non-redundant) bins
  int N  = (F-1) * H + B/2;               // number of samples
  int M  = (K-1) * 2;                     // FFT size
  int k0 = (M-B) / 2;                     // read start in resynthesized grain before the 
                                          // resynthesis window is applied
  std::vector<double> y(N);               // allocate signal
  std::vector<double> g(B);               // grain
  rsComplexDbl *Y = new rsComplexDbl[M];  // short-time spectrum

  int i, k;
  for(i = 0; i < F; i++)
  {
    // create symmetrized FFT spectrum and transform to time domain:
    Y[0] = s(i, 0);
    for(k = 1; k < K; k++)
    {
      Y[k]   = s(i, k);
      Y[M-k] = Y[k].getConjugate();
    }
    rsIFFT(Y, M);

    // apply synthesis-window and overlap/add into output signal:
    for(k = 0; k < B; k++)
      g[k] = Y[k0+k].re * w[k];  // k0 != 0, if zero-padding was used
    rsAddInto(y.data(), N, g.data(), B, i*H-B/2);
  }

  delete[] Y;
  return y;
}

std::vector<double> rsPhaseVocoder::getModulation(double *wa, double *ws, int B, int H, int F)
{
  // wa: analysis window, ws: synthesis-window, B: block size, H: hop size, F: number of frames
  int N = (F-1) * H + B/2;         // number of samples
  std::vector<double> y(N);        // modulation signal
  double *w = new double[B];       // product-window
  for(int n = 0; n < B; n++)
    w[n] = wa[n] * ws[n];             
  for(int i = 0; i < F; i++)
    rsAddInto(y.data(), N, w, B, i*H-B/2);
  delete[] w;
  return y;
}

// Misc:

void rsPhaseVocoder::init()
{
  //fs = 44100.0;   // samplerate
  //Nb = 512;       // block size
  //Nh = Nb/2;      // hop size
  //Nf = Nb*2;      // FFT size
}

