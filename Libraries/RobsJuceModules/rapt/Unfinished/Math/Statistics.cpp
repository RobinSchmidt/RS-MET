

template<class T>
void rsCrossCorrelationDirect(T x[], T y[], int N, T r[])
{
  for(int k = 0; k < N; k++)
  {
    r[k] = 0;
    for(int n = 0; n < N-k; n++)
      r[k] += x[n+k] * y[n]; // for complex sequences, y must be conjugated?
  }
}

template<class T>
void rsCrossCorrelationFFT(T x[], T y[], int N, T r[])
{
  // create zero-padded sequences such that no foldover occurs (use a power of 2 for the FFT) and
  // obtain the two FFT spectra:
  int Np = 2*rsNextPowerOfTwo(N);
  std::complex<T> *X = new std::complex<T>[Np];
  std::complex<T> *Y = new std::complex<T>[Np];
  rsArrayTools::convert(x, X, N);
  rsArrayTools::convert(y, Y, N);
  rsFFT(X, Np);
  rsFFT(Y, Np);

  // compute cross power spectrum, store it in X (see text below Eq.12.0.11 in: Numerical
  // Recipies in C, 2nd Ed.):
  int n;
  X[0] = X[0] * Y[0];
  for(n = 1; n < Np; n++)
    X[n] = X[n] * Y[Np-n];

  // obtain cross correlation by IFFT, copy to output and cleanup:
  rsIFFT(X, Np);
  for(n = 0; n < N; n++)
    r[n] = X[n].real();
  delete[] X;
  delete[] Y;
}

template<class T>
void rsCrossCorrelation(T x[], T y[], int N, T r[], bool removeBias)
{
  if(N < 64) // ad hoc - measure at which point FFT becomes more efficient
    rsCrossCorrelationDirect(x, y, N, r);
  else
    rsCrossCorrelationFFT(x, y, N, r);

  if(removeBias == true)
    rsRemoveCorrelationBias(r, N, r);
}

template<class T>
void rsAutoCorrelationFFT(T x[], int N, T r[])
{
  int Np = 2*rsNextPowerOfTwo(N);
  std::complex<T> *X = new std::complex<T>[Np];
  rsArrayTools::convert(x, X, N);
  rsFFT(X, Np);
  int n;
  X[0] = X[0] * X[0];
  for(n = 1; n < Np; n++)
    X[n] = X[n] * conj(X[n]);
  // We can't use X[n]*X[Np-n] because the X-array is getting messed inside the loop itself, so
  // we use Eq.12.0.11 as is.
  rsIFFT(X, Np);
  for(n = 0; n < N; n++)
    r[n] = X[n].real();
  delete[] X;
}

template<class T>
void rsRemoveCorrelationBias(T x[], int N, T r[])
{
  for(int k = 0; k < N; k++)
    r[k] = N*x[k]/(N-k);
}

template<class T>
T rsCrossCorrelation(const T *x, int Nx, const T *y, int Ny)
{
  T xx = rsArrayTools::sumOfSquares(x, Nx);
  T yy = rsArrayTools::sumOfSquares(y, Ny);
  T xy = rsArrayTools::sumOfProducts(x, y, rsMin(Nx, Ny));
  if(xx == 0 || yy == 0)
    return 0;
  return xy / sqrt(xx*yy);
  // Assuming that we would divide each sum by the same value, namely max(Nx,Ny), for computing 
  // a mean, the divisor cancels out in the cross-correlation formula
}

template<class T>
T rsStretchedCrossCorrelation(const T *x, int Nx, const T *y, int Ny)
{
  // x should be the longer than y:
  if(Nx < Ny)
    return rsStretchedCrossCorrelation(y, Ny, x, Nx); // recursive call with swapped arguments
  if(Nx == Ny)
    return rsCrossCorrelation(x, Nx, y, Ny);          // function for same-length arrays

  double a = double(Ny) / double(Nx);  // readout speed for y-array
  int    nx;                           // index into the x-array
  double xn;                           // value of x at nx
  T yn;                                // (interpolated) value of y at a*nx
  T xx = 0;                            // sum of squares of x-values   
  T yy = 0;                            // sum of squares of (interpolated) y-values
  T xy = 0;                            // sum of x*y products
  for(nx = 0; nx < Nx; nx++)
  {
    xn  = x[nx];
    yn  = rsArrayTools::interpolatedValueAt(y, Ny, a*nx);
    xx += xn*xn;
    yy += yn*yn;
    xy += xn*yn;
  }
  if(xx == 0 || yy == 0)
    return 0;
  return xy / sqrt(xx*yy);

  // Is the use of double instead of T for a and xn intentional? If so, document, why. It seems
  // rsArrayTools::interpolatedValueAt also expects a double for the interpolated position - is 
  // this the reason? If so, why does interpolatedValueAt expect double? Maybe we should use TSig
  // for the signal and TPos for the position. I actually think, it makes sense when T is some
  // vector-type, like for stereo signals. What about complex signals? Should we use complex 
  // conjugation here in the code somewhere? Currently, this function should be used only for
  // real signals
}