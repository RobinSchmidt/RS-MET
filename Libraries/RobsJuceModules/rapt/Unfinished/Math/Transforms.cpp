template<class T>
void smbFft(T *fftBuffer, long fftFrameSize, long sign)
{
  T wr, wi, arg, *p1, *p2, temp;
  T tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
  long i, bitm, j, le, le2, k, logN;

  logN = (long)(log((double)fftFrameSize)/log(2.)+.5); // pass this value as parameter

  // bit-reversed ordering (?):
  for(i = 2; i < 2*fftFrameSize-2; i += 2)
  {
    for(bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1)
    {
      if(i & bitm)
        j++;
      j <<= 1;
    }
    if(i < j)
    {
      p1      = fftBuffer+i;
      p2      = fftBuffer+j;
      temp    = *p1;
      *(p1++) = *p2;
      *(p2++) = temp;
      temp    = *p1;
      *p1     = *p2;
      *p2     = temp;
    }
  }

  // the actual FFT:
  for(k = 0, le = 2; k < logN; k++)
  {
    le  <<= 1;
    le2   = le>>1;
    ur    = 1.0;
    ui    = 0.0;
    arg   = (T)PI/(le2>>1);
    wr    =      cos(arg);
    wi    = sign*sin(arg);     // use trigonometric recursion (maybe with extended precision)
    for(j = 0; j < le2; j += 2)
    {
      p1r = fftBuffer+j;
      p1i = p1r+1;
      p2r = p1r+le2;
      p2i = p2r+1;
      for(i = j; i < 2*fftFrameSize; i += le)
      {
        tr    = *p2r * ur - *p2i * ui;
        ti    = *p2r * ui + *p2i * ur;
        *p2r  = *p1r - tr;
        *p2i  = *p1i - ti;
        *p1r += tr;
        *p1i += ti;
        p1r  += le;
        p1i  += le;
        p2r  += le;
        p2i  += le;
      }
      tr = ur*wr - ui*wi;
      ui = ur*wi + ui*wr;
      ur = tr;
    }
  }
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
  // perhaps not very advisable from a numeric point of view but very nice mathematically. 

  int np = 1;          // NumOfProblems
  int h  = N/2;        // HalfSize -> distance between butterflied values?
  int jf;              // JFirst
  int jl;              // JLast
  std::complex<T> tmp; // Temp
  std::complex<T> wj;  // W (current twiddle factor), maybe rename to wjk or wkj
  std::complex<T> wm = exp(std::complex<T>(T(0), T(-2.0*PI/N))); // multiplier for twiddle factor, todo: use polar() 

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

/*
template<class T>
void rsBluesteinFFT(std::complex<T> *a, int N)
{
  rsError("Not yet implemented.");
}
*/

// \todo explicit template instantiations for float and double
//....
