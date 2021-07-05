template<class T>
void smbFft(T *fftBuffer, long fftFrameSize, long sign)
{
  T wr, wi, arg, *p1, *p2, temp;
  T tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
  long i, bitm, j, le, le2, k, logN;

  logN = (long)(log((double)fftFrameSize)/log(2.)+.5); 
  // todo: pass this value as parameter or use a special rsLog2Int function, see comments below
  // rsLinearTransforms::fourierRadix2DIF

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

// move to prototypes - this is really not meant to be used anywhere in production code:
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
  {
    //rsRadix2FFT(a, N);
    rsLinearTransforms::fourierRadix2DIF(a, N);
  }
  else
  {
    rsError("Arbitrary length FFT not yet implemented");
    // rsBluesteinFFT(a, N);
  }
}

template<class T>
void rsIFFT(std::complex<T> *a, int N)
{
  rsLinearTransforms::fourierInvRadix2DIF(a, N);
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

/*
template<class T>
void rsBluesteinFFT(std::complex<T> *a, int N)
{
  rsError("Not yet implemented.");
}
*/

template<class T>
static void rsLinearTransforms::fourierRadix2DIF(T* a, int N, T W)
{
  rsAssert(rsIsPowerOfTwo(N), "N must be a power of 2");
  int n = 1;          // NumOfProblems
  int h = N/2;        // HalfSize -> distance between butterflied values?
  while(h > 0) {                       // loop over the problems(?)
    for(int k = 0; k < n; k++) {       // loop over the sub-FFTs
      T Wjk  = rsUnityValue(W);        // init twiddle factor W^(j*k) for current sub-FFT to 1
      int jf = 2*k*h;                  // first index for k-th sub-FFT, JFirst
      int jl = jf+h-1;                 // last index for k-th sub-FFT, JLast
      for(int j = jf; j <= jl; j++) {  // loop over the values
        T aj   = a[j];                 // temporary
        a[j]   =  aj + a[j+h];         // upper wing of Gentleman-Sande butterfly
        a[j+h] = (aj - a[j+h])*Wjk;    // lower wing
        Wjk *= W; }}                   // update twiddle factor for next iteration
    n *= 2;           // next stage has twice as many sub-FFTs
    h /= 2;           // distance between butterflied values halves
    W *= W; }         // twiddle factor rotates twice as fast in next stage, W_N = W_(N/2)?
  rsArrayTools::orderBitReversed(a, N, (int)(rsLog2(N)+0.5)); // descramble outputs
}
// ToDo:
// -use a special, optimized rsLog2Int function, using a check for the index of the highest bit:
//  https://stackoverflow.com/questions/671815/what-is-the-fastest-most-efficient-way-to-find-the-highest-set-bit-msb-in-an-i
//  Answer 4: "The Debruin technique should only be used when the input is already a power of two"
//  or: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
//  or: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
//      this has nice 2 functions for uint32 and uint64
// -maybe use bit-shifts instead of n *= 2, h /= 2, h = N/2

template<class T>
void rsLinearTransforms::hadamard(T* A, int N)
{
  rsAssert(rsIsPowerOfTwo(N), "N must be a power of 2");
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

template<class T>
void rsLinearTransforms::kronecker2x2(T* A, int N, T a, T b, T c, T d)
{
  rsAssert(rsIsPowerOfTwo(N), "N must be a power of 2");
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
void rsLinearTransforms::kronecker3x3(T* v, int N, const rsMatrix3x3<T>& A)
{
  rsAssert(rsIsPowerOfN(N, 3), "N must be a power of 3");
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 3*h) {
      for(int j = i; j < i+h; j++) {
        T x = v[j+0*h];
        T y = v[j+1*h];
        T z = v[j+2*h];
        v[j+0*h] = A(0,0) * x + A(0,1) * y + A(0,2) * z;
        v[j+1*h] = A(1,0) * x + A(1,1) * y + A(1,2) * z;
        v[j+2*h] = A(2,0) * x + A(2,1) * y + A(2,2) * z; }}
    h *= 3;  }
}

//-------------------------------------------------------------------------------------------------
/*

ToDo:
-Wavelet transforms: Gabor, Daubechies, Haar, ...
 https://en.wikipedia.org/wiki/Haar_wavelet
-Walsh, Householder, Givens, Toeplitz (?), Slant (uses sawtooths as basis 
 functions)..what about using rectangular waves? is that what the Hadamard trafo does? ..find a way
 to plot the basis vectors of the transforms - set one coeff one and all others zero and perform an 
 inverse transform - this should give the basis vector corresponding to the coeff that was set to 1
 https://en.wikipedia.org/wiki/Walsh_function
 https://en.wikipedia.org/wiki/Dyadic_transformation (sawtooth?)
-KarhunenLoeve/Eigenvector/PCA 
 http://fourier.eng.hmc.edu/e161/lectures/klt/node3.html
 https://link.springer.com/chapter/10.1007/978-3-540-72943-3_10
-Add variant of the Kronecker transforms that use an additional workspace array and avoid the
 internal temp variables - then benchmark them against the in-place algo. Maybe, these out-of-place
 should first go into prototypes, and be added to the library only when it is found that they are
 faster
-DCT, DST, Mellin
-test the FFT algo for finite fields by doing a fast FFT convolution of two sequences of modular 
 integers and comparing that to the result of a naive convolution, see:
 https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field
-implement arbitrary length FFT using the Bluestein algorithm, test it for complex and modular 
 integer types
-implement a decimation in time FFT
-maybe implement a radix-3 FFT because it may be useful to have FFTs based on different primes
 because when we want to do NTTs, because (i think), the radix is required to be coprime with 
 length of the array (otherwise N-th roots of unity don't exist, or something). Or does it need to
 be coprime to the modulus? or both? -> figure out..or maybe the array length must be coprime
 with the modulus? ...we have 3 numbers: length, modulus, radix
-can we also devise a Bluestein-like algorithm for sequences of arbitrary length for modular 
 integers? 
-what about a tranform based on the outer product of a vector: y = A*x, where A = v * v^T
 -> A is also a Kronecker product, so it should be possible to invert by the same algo but having
    elements of v replaced by their reciprocals
 -> what, if we choose v = x? Can we still invert it simply if v (i.e. x) is unknown?
 -> maybe try C = E{x * x^T} where E is the expectation value, C is the covariance matrix, i think
    this matrix can be formed as outer product of the autocorrelation vector of x
-Implement radix2-DIT algo and also radix-3 algo(s) ...maybe we can do a radix-3 NTT with a modulus
 that is a power of 2?
-Implement sequency based ordering for the Hadamard trafo:
   https://en.wikipedia.org/wiki/Walsh_matrix#Sequency_ordering
   https://en.wikipedia.org/wiki/Gray_code
 "The sequency ordering of the rows of the Walsh matrix can be derived from the ordering of the 
  Hadamard matrix by first applying the bit-reversal permutation and then the Gray-code 
  permutation"

Ideas:
-Wavelet transforms split the signal into lowpass and highpass part at each stage, leading to 
 power-of-2 lengths. Maybe try to split into 3 parts at each stage (low/mid/high), maybe using
 something like: yL = (xL+xM+xR)/3, yM = (xR-xL)/2, yH = xM-(xL+xR)/2. where xL,xL,xR are left
 middle, right inputs and yL,yM,yH are low/mid/high outputs. Idea: yL approximates the value, yM
 the 1st derivative, yH the 2nd derivative all at tM which is the time instant of xM. This would 
 lead to power-of-3 lengths. Maybe generalize the idea to splitting into N parts, where the 
 formulas are based on N-th order derivative estimators. What about the integral? That would seem 
 to require to make use of data from the previous block to init an integrator formula at the left 
 end of the block, for example: yM = (xL+xM+xR)/3, yL = yL[tM-3] + yM, yH = (xR-xL)/2 - but that's
 just a constant amount (a single value) of extra storage per block - maybe call it DIV-trafo 
 (derivative, integral, value)...or VID...maybe it's good for video, indeed :-) ..but i think, a 
 good trafo for video should make use of spatial neighbours, too, not only of temporal ones. Maybe
 try that with rendered videos, like the SIRP model. Maybe make it lossy in the last of the 8 bits
 per pixel, amounting to dithering...for video, dithering should be spatial, not temporal...maybe
 something based on Floyd-Steinberg algorithm
-make a class rsTimeFrequencyRepresentation using various transforms

-can fourierRadix2DIF vectorized by having WN a vector of [W W^2 W^4 W^8]?

*/

