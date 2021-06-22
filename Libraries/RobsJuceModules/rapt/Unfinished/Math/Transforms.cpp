template<class T>
void smbFft(T *fftBuffer, long fftFrameSize, long sign)
{
  T wr, wi, arg, *p1, *p2, temp;
  T tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
  long i, bitm, j, le, le2, k, logN;

  logN = (long)(log((double)fftFrameSize)/log(2.)+.5); // pass this value as parameter

  // Bit-reversed ordering (?):
  for(i = 2; i < 2*fftFrameSize-2; i += 2) {
    for(bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
      if(i & bitm)
        j++;
      j <<= 1; }
    if(i < j) {
      p1      = fftBuffer+i;
      p2      = fftBuffer+j;
      temp    = *p1;
      *(p1++) = *p2;
      *(p2++) = temp;
      temp    = *p1;
      *p1     = *p2;
      *p2     = temp;  }}

  // The actual FFT:
  for(k = 0, le = 2; k < logN; k++) {
    le  <<= 1;
    le2   = le>>1;
    ur    = 1.0;
    ui    = 0.0;
    arg   = (T)PI/(le2>>1);
    wr    =      cos(arg);
    wi    = sign*sin(arg);
    for(j = 0; j < le2; j += 2) {  
      p1r = fftBuffer+j;
      p1i = p1r+1;
      p2r = p1r+le2;
      p2i = p2r+1;
      for(i = j; i < 2*fftFrameSize; i += le) {
        tr    = *p2r * ur - *p2i * ui;
        ti    = *p2r * ui + *p2i * ur;
        *p2r  = *p1r - tr;
        *p2i  = *p1i - ti;
        *p1r += tr;
        *p1i += ti;
        p1r  += le;
        p1i  += le;
        p2r  += le;
        p2i  += le;  }
      tr = ur*wr - ui*wi;
      ui = ur*wi + ui*wr;
      ur = tr;  }}
}

template<class T>
void rsDFT(std::complex<T> *x, int N) // maybe have a boolean "inverse" parameter
{
  std::complex<T> *X = new std::complex<T>[N];
  for(int k = 0; k < N; k++)
  {
    X[k] = 0;
    for(int n = 0; n < N; n++)
      X[k] += x[n]*exp(std::complex<T>(T(0), T(-2.0*PI*n*k/N)));
  }
  rsArrayTools::copy(X, x, N);
  delete[] X;
}

template<class T>
void rsFFT(std::complex<T> *a, int N)
{
  if(rsIsPowerOfTwo(N))
    rsRadix2FFT(a, N);
  else
  {
    rsError("Arbitrary length FFT not yet implemented");
    // rsBluesteinFFT(a, N);
  }
}

template<class T>
void rsIFFT(std::complex<T> *a, int N)
{
  // see "Understanding Digital Signal Processing", page 427:
  rsConjugate(a, N);
  rsFFT(a, N);
  rsConjugate(a, N);
  rsArrayTools::scale(a, N, T(1)/N);
  // \todo check, if this works for arbitrary (non power-of-two) N
  // shouldn't it be possible to just use a "sign" in the core FFT routine that just switches the
  // exponent of the twiddle factor to be positive or negative (->more efficient, we don't need to
  // iterate through the arrays - but i'm not sure, if it works with arbitrary lengths)
}

template<class T>
void rsMagnitudeAndPhase(T *signal, int N, T *magnitudes, T *phases)
{
  std::complex<T> *tmp = new std::complex<T>[N];

  rsArrayTools::convert(signal, tmp, N);
  rsFFT(tmp, N);

  int k;
  for(k = 0; k < N/2; k++)
    magnitudes[k] = abs(tmp[k]);
  // the bin at 0 (and N/2-1?) should be multiplied by 2?

  if(phases != NULL)
  {
    for(k = 0; k < N/2; k++)
      phases[k] = arg(tmp[k]);
  }
  delete[] tmp;
}

template<class T>
void rsRadix2FFT(std::complex<T> *a, int N)
{
  // The algorithm was ported from algorithm 4.2 in the book "Inside the FFT black box" and 
  // then simplified. All twiddle factors are computed on the fly via recursion. There's just one
  // single complex exponential (i.e. a sin/cos pair) evaluation in the whole routine. That's 
  // perhaps not very advisable from a numeric point of view but very nice mathematically. But for
  // number theoretic transform (i.e. FFT using the roots of unity of modular arithmetic), that's
  // no issue, so that seems desirable when the algo is adapted for that.

  int np = 1;          // NumOfProblems
  int h  = N/2;        // HalfSize -> distance between butterflied values?
  int jf;              // JFirst
  int jl;              // JLast
  std::complex<T> tmp; // Temp
  std::complex<T> wj;  // W (current twiddle factor), maybe rename to wjk or wkj
  std::complex<T> wm = exp(std::complex<T>(T(0), T(-2.0*PI/N))); 
  // Multiplier for twiddle factor, todo: use polar(), maybe using +2.0 instead of -2.0 will give 
  // the inverse trafo?

  // \todo use setRadiusAndAngle for initializing wm (more efficient), we also do not need to 
  // include the ComplexFunctions.inl

  // \todo maybe use long double precision for twiddle factors to avoid excessive error 
  // accumulation in the recursions

  while(h > 0)
  {
    for(int k = 0; k < np; k++)        // loop over the sub-FFTs
    {
      wj = std::complex<T>(1.0, 0.0);  // init twiddle factor for current sub-FFT
      jf = 2*k*h;                      // first index for k-th sub-FFT
      jl = jf+h-1;                     // last index for k-th sub-FFT
      for(int j = jf; j <= jl; j++)    // loop over the values
      {
        tmp    = a[j];
        a[j]   =     tmp+a[j+h];       // upper wing of Gentleman-Sande butterfly
        a[j+h] = wj*(tmp-a[j+h]);      // lower wing
        wj *= wm;                      // update twiddle factor for next iteration
      }
    }
    np *= 2;  // next stage has twice as much sub-FFTs
    h  /= 2;  // distance between butterflied values halves - maybe use bit-shift
    wm *= wm; // twiddle factor rotates twice as fast in next stage
  }

  rsArrayTools::orderBitReversed(a, N, (int)(rsLog2(N)+0.5)); // descramble outputs
}
// ToDo: adapt algo for finite fields - see:
// https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field
// it may have to take the n-th root of unity as argument...i guess that's our initial value for 
// wm? ...maybe then, the function can be take as is - the template parameter should not be the 
// real type T that underlies the complex type, but the complex type itself, which can then be
// replaced by something like rsModularInteger? -> test it to do a fast FFT convolution of two 
// sequences of modular integers...can we also devise a Bluestein-like algorithm for sequences
// of arbitrary length for modular integers?

/*
template<class T>
void rsBluesteinFFT(std::complex<T> *a, int N)
{
  rsError("Not yet implemented.");
}
*/

template<class T>
void rsTransforms::kronecker2x2(T* A, int N, T a, T b, T c, T d)
{
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 2*h) {
      for(int j = i; j < i+h; j++) {
        T x = A[j];
        T y = A[j+h];
        A[j]   = a*x + b*y;
        A[j+h] = c*x + d*y;  }}
    h *= 2;  }
}

template<class T>
void rsTransforms::hadamard(T* A, int N)
{
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 2*h) {
      for(int j = i; j < i+h; j++) {
        T x = A[j];
        T y = A[j+h];
        A[j]   = x + y;
        A[j+h] = x - y;  }}
    h *= 2;  }
}


//-------------------------------------------------------------------------------------------------

// Deprecated alias names:
template<class T> 
void rsFGHT(T* A, int N, T a, T b, T c, T d) { rsTransforms::kronecker2x2(A, N, a, b, c, d); }


//-------------------------------------------------------------------------------------------------
/*

ToDo:
-Wavelet transforms: Gabor, Daubechies, Walsh, ...
-Fourier, Hadamard, Walsh, Householder, Givens, Toeplitz (?)
-KarhunenLoeve/Eigenvector/PCA 
 http://fourier.eng.hmc.edu/e161/lectures/klt/node3.html
 https://link.springer.com/chapter/10.1007/978-3-540-72943-3_10
-Add a variants of the Kronecker transforms that use an additional workspace array and avoid the
 internal temp variables - then benchmark them against the in-place algo. Maybe, these out-of-place
 should first go into prototypes, and be added to the library only when it is found that they are
 faster

*/

