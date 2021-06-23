#include "OouraFFT/fft4g.cpp"

// construction/destruction:

template<class T>
rsFourierTransformerRadix2<T>::rsFourierTransformerRadix2()
{
  N                   = 0;
  logN                = 0;
  direction           = FORWARD;
  normalizationMode   = NORMALIZE_ON_INVERSE_TRAFO;
  normalizationFactor = 1.0;
  w                   = NULL;
  ip                  = NULL;
  tmpBuffer           = NULL;

  setBlockSize(256);
}

template<class T>
rsFourierTransformerRadix2<T>::~rsFourierTransformerRadix2()
{
  // free dynamically allocated memory:
  if( w != NULL )
    delete[] w;
  if( ip != NULL )
    delete[] ip;
  if( tmpBuffer != NULL )
    delete[] tmpBuffer;
}

// setup:

template<class T>
void rsFourierTransformerRadix2<T>::setBlockSize(int newBlockSize)
{
  // check new blocksize for validity:
  if( newBlockSize >= 2 && rsIsPowerOfTwo(newBlockSize) )
  {
    // check, if the new blocksize is actually different from the old one in order to avoid
    // unnecesarry re-allocations and re-computations:
    if( newBlockSize != N )
    {
      N    = newBlockSize;
      logN = (int) floor( rsLog2((T) N + 0.5 ) );
      updateNormalizationFactor();

      if( w != NULL )
        delete[] w;
      w    = new T[2*N];

      if( ip != NULL )
        delete[] ip;
      ip    = new int[(int) ceil(4.0+rsSqrt((T)N))];
      ip[0] = 0; // indicate that re-initialization is necesarry

      if( tmpBuffer != NULL )
        delete[] tmpBuffer;
      tmpBuffer = new std::complex<T>[N];
    }
  }
  else if( !rsIsPowerOfTwo(newBlockSize) || newBlockSize <= 1 )
    RS_DEBUG_BREAK; // this class can only deal with blocksizes >= 2 that are a power of two
}

template<class T>
void rsFourierTransformerRadix2<T>::setDirection(int newDirection)
{
  if( newDirection >= FORWARD && newDirection <= INVERSE )
  {
    // only when the new direction is actually different form the old one, we have to conjugate
    // all the twiddle-factors, otherwise everything must stay as is:
    if( newDirection != direction )
    {
      direction = newDirection;
      updateNormalizationFactor();
    }
  }
  else
    RS_DEBUG_BREAK; // passed int-parameter does not correspond to any meaningful enum-field
}

template<class T>
void rsFourierTransformerRadix2<T>::setNormalizationMode(int newNormalizationMode)
{
  if( newNormalizationMode >= NORMALIZE_ON_FORWARD_TRAFO &&
      newNormalizationMode <= NEVER_NORMALIZE )
  {
    normalizationMode = newNormalizationMode;
    updateNormalizationFactor();
  }
  else
    RS_DEBUG_BREAK; // passed int-parameter does not correspond to any meaningful enum-field
}

/*
void rsFourierTransformerRadix2::setRealSignalMode(bool willBeUsedForRealSignals)
{
  ip[0] = 0; // retriggers twiddle-factor computation
}
*/

// signal processing:

template<class T>
void rsFourierTransformerRadix2<T>::transformComplexBufferInPlace(std::complex<T> *buffer)
{
  // retrieve the adresses of the real part of the first array entries in order to treat the
  // Complex arrays as arrays of two successive T-numbers:
  //T* d_buffer = &(buffer[0].real());
  T* d_buffer = (T*) &buffer[0];

  // normalize the FFT-input, if required:
  if( normalizationFactor != 1.0 )
  {
    for(int n=0; n<2*N; n++)
      d_buffer[n] *= normalizationFactor;
  }

  // use Ooura's routine:
  int sign;
  if( direction == FORWARD )
    sign = -1;
  else
    sign = +1;
  cdft(2*N, sign, d_buffer, ip, w);
}

template<class T>
void rsFourierTransformerRadix2<T>::transformComplexBuffer(std::complex<T> *inBuffer,
                                                        std::complex<T> *outBuffer)
{
  // retrieve the adresses of the real part of the first array entries in order to treat the
  // Complex arrays as arrays of two successive T-numbers:
  T* d_inBuffer  = (T*) &inBuffer[0];
  T* d_outBuffer = (T*) &outBuffer[0];

  // copy the input into the output for the in-place routine (thereby normalize, if necesarry):
  int n;
  if( normalizationFactor != 1.0 )
  {
    for(n=0; n<2*N; n++)
      d_outBuffer[n] = d_inBuffer[n] * normalizationFactor;
  }
  else
  {
    for(n=0; n<2*N; n++)
      d_outBuffer[n] = d_inBuffer[n];
  }

  // use Ooura's routine:
  int sign;
  if( direction == FORWARD )
    sign = -1;
  else
    sign = +1;
  cdft(2*N, sign, d_outBuffer, ip, w);
}

// convenience functions for real signals:

template<class T>
void rsFourierTransformerRadix2<T>::transformRealSignal(T *inSignal, std::complex<T> *outSpectrum)
{
  setDirection(FORWARD);

  // retrieve the adress of the real part of the first array entry of the output array in order to
  // treat the Complex array as array of two successive T-numbers:
  T* d_outBuffer = (T*) &outSpectrum[0];

  // copy the input into the output for the in-place routine (thereby normalize, if necesarry):
  int n;
  if( normalizationFactor != 1.0 )
  {
    for(n=0; n<N; n++)
      d_outBuffer[n] = inSignal[n] * normalizationFactor;
  }
  else
  {
    for(n=0; n<N; n++)
      d_outBuffer[n] = inSignal[n];
  }

  // use Ooura's routine:
  rdft(N, 1, d_outBuffer, ip, w);

  // for some reason, this routine returns the second half of the spectrum (the complex conjugate
  // values of the desired first half), so we need to take the complex conjugates:
  for(n=3; n<N; n+=2) // start at n=3 (imaginary part of the first bin after DC)
    d_outBuffer[n] = -d_outBuffer[n];
}

template<class T>
void rsFourierTransformerRadix2<T>::transformRealSignal(T *signal, T *reAndIm)
{
  std::complex<T>* c_reAndIm = (std::complex<T>*) &(reAndIm[0]);
  transformRealSignal(signal, c_reAndIm);
}

template<class T>
void rsFourierTransformerRadix2<T>::getRealSignalMagnitudesAndPhases(
  T *signal, T *magnitudes, T *phases)
{
  transformRealSignal(signal, tmpBuffer);

  // store the two purely real transform values at DC and Nyquist-frequency in the first fields of
  // the magnitude- and phase- arrays respectively:
  magnitudes[0] = tmpBuffer[0].real();
  phases[0]     = tmpBuffer[0].real();  // wait - what? this seems wrong should this be tmpBuffer[N/2]
                                        // or something?
  // actually, magnitudes should always be a positive value - and taking the real part does not ensure 
  // this - we should probably set the phase to pi, when tmpBuffer[0] is negative

  // fill the rest of the array with the magnitudes and phases of the regular bins:
  T* dBuffer = (T*) &tmpBuffer[0];
  T  re, im;
  int k;
  for(k=1; k<N/2; k++)
  {
    re = dBuffer[2*k];
    im = dBuffer[2*k+1];
    magnitudes[k] = rsSqrt(re*re + im*im);
    if( re == 0.0 && im == 0.0 )
      phases[k] = 0.0;
    else
      phases[k] = atan2(im, re);
  }
}

template<class T>
void rsFourierTransformerRadix2<T>::getRealSignalMagnitudes(T *signal, T *magnitudes)
{
  transformRealSignal(signal, tmpBuffer);
  magnitudes[0] = tmpBuffer[0].real();

  T* dBuffer = (T*) &tmpBuffer[0];
  T  re, im;
  int     k;
  for(k=1; k<N/2; k++)
  {
    re            = dBuffer[2*k];
    im            = dBuffer[2*k+1];
    magnitudes[k] = rsSqrt(re*re + im*im);
  }
}

template<class T>
void rsFourierTransformerRadix2<T>::transformSymmetricSpectrum(std::complex<T> *inSpectrum,
                                                            T *outSignal)
{
  setDirection(INVERSE);

  // retrieve the adress of the real part of the first array entry of the output array in order to
  // treat the Complex array as array of two successive T-numbers:
  T* d_inBuffer = (T*) &inSpectrum[0];

  // copy the input into the output for the in-place routine (thereby normalize, if necesarry):
  int n;
  if( normalizationFactor != 1.0 )
  {
    for(n=0; n<N; n++)
      outSignal[n] = T(2.0) * d_inBuffer[n] * normalizationFactor;
  }
  else
  {
    for(n=0; n<N; n++)
      outSignal[n] = T(2.0) * d_inBuffer[n];
  }

  // for some reason, the subsequent routine expects the second half of the spectrum (the complex
  // conjugate values of the first half), so we need to take the complex conjugates:
  for(n=3; n<N; n+=2) // start at n=3 (imaginary part of the first bin after DC)
    outSignal[n] = -outSignal[n];

  // use Ooura's routine:
  rdft(N, -1, outSignal, ip, w);
}

template<class T>
void rsFourierTransformerRadix2<T>::transformSymmetricSpectrum(T *reAndIm, T *signal)
{
  std::complex<T>* c_reAndIm = (std::complex<T>*) &(reAndIm[0]);
  transformSymmetricSpectrum(c_reAndIm, signal);
}

template<class T>
void rsFourierTransformerRadix2<T>::getRealSignalFromMagnitudesAndPhases(T *magnitudes,
                                                                      T *phases,
                                                                      T *signal)
{
  tmpBuffer[0].real(magnitudes[0]);
  tmpBuffer[0].real(phases[0]);      // ?BUG? should it be assigned to imag?

  int k;
  T* dBuffer = (T*) &tmpBuffer[0];
  T  s, c;
  for(k=1; k<N/2; k++)
  {
    rsSinCos(phases[k], &s, &c);
    dBuffer[2*k]   = magnitudes[k] * c;
    dBuffer[2*k+1] = magnitudes[k] * s;
  }

  transformSymmetricSpectrum(tmpBuffer, signal);
}

template<class T>
T rsFourierTransformerRadix2<T>::getNormalizationFactor(int N, int dir, int mode)
{
  if( (mode == NORMALIZE_ON_FORWARD_TRAFO && dir == FORWARD) ||
      (mode == NORMALIZE_ON_INVERSE_TRAFO && dir == INVERSE)    )
    return T(1.0) / (T) N;
  else if( mode == ORTHONORMAL_TRAFO )
    return T(1.0) / rsSqrt((T) N);
  else
    return T(1.0);
}

// pre-calculations:

template<class T>
void rsFourierTransformerRadix2<T>::updateNormalizationFactor()
{
  normalizationFactor = getNormalizationFactor(N, direction, normalizationMode);
}

//=================================================================================================

// construction/destruction:

template<class T>
rsFourierTransformerBluestein<T>::rsFourierTransformerBluestein()
{
  h                     = NULL;
  c                     = NULL;
  y                     = NULL;
  N                     = 0;
  M                     = 0;
  direction             = rsFourierTransformerRadix2<T>::FORWARD;
  normalizationMode     = rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO;
  normalizationFactor   = 1.0;
  blockSizeIsPowerOfTwo = false;
}

template<class T>
rsFourierTransformerBluestein<T>::~rsFourierTransformerBluestein()
{
  if( h != NULL )
    delete[] h;
  if( c != NULL )
    delete[] c;
  if( y != NULL )
    delete[] y;
}

// parameter settings:

template<class T>
void rsFourierTransformerBluestein<T>::setBlockSize(int newBlockSize)
{
  // check new blocksize for validity and update our related member variables when it is valid (and
  // different from the old one):
  if( newBlockSize > 0 && newBlockSize != N )
  {
    N = newBlockSize;

    // do the preprocessing for the Bluestein algorithm only when we really need it:
    blockSizeIsPowerOfTwo = rsIsPowerOfTwo(N);
    if( !blockSizeIsPowerOfTwo )
    {
      M = rsNextPowerOfTwo(2*N-1);
      generateChirp();
    }
    else
      M = N;

    updateNormalizationFactor();

    // free old and allocate new memory for the internal buffer of size M:
    if( y != NULL )
      delete[] y;
    y = new std::complex<T>[M];
  }
}

template<class T>
void rsFourierTransformerBluestein<T>::setDirection(int newDirection)
{
  if( newDirection >= rsFourierTransformerRadix2<T>::FORWARD &&
      newDirection <= rsFourierTransformerRadix2<T>::INVERSE )
  {
    // only when the new direction is actually different form the old one, we have to conjugate
    // all the twiddle-factors, otherwise everything must stay as is:
    if( newDirection != direction )
    {
      direction = newDirection;
      updateNormalizationFactor();

      // we need different chirps for forward and inverse transform:
      if( !blockSizeIsPowerOfTwo )
        generateChirp();
    }
  }
  else
    rsError("passed int-parameter does not correspond to any meaningful enum-field");
}

template<class T>
void rsFourierTransformerBluestein<T>::setNormalizationMode(int newNormalizationMode)
{
  if( newNormalizationMode >= rsFourierTransformerRadix2<T>::NORMALIZE_ON_FORWARD_TRAFO &&
      newNormalizationMode <= rsFourierTransformerRadix2<T>::NEVER_NORMALIZE )
  {
    normalizationMode = newNormalizationMode;
    updateNormalizationFactor();
  }
  else
    rsError("passed int-parameter does not correspond to any meaningful enum-field");
}

// signal processing:

template<class T>
void rsFourierTransformerBluestein<T>::transformComplexBufferInPlace(std::complex<T> *buffer)
{
  // use the embedded FourierTransformerRadix2-object directly on the input data for the special
  // case where the blockSize is a power of two:
  if( blockSizeIsPowerOfTwo )
  {
    transformerRadix2.setBlockSize(N);
    transformerRadix2.setDirection(direction);
    transformerRadix2.setNormalizationMode(normalizationMode);
    transformerRadix2.transformComplexBufferInPlace(buffer);
    return;
  }

  // O.K. at this point we know that the block-size is not a power of of two and we must employ the
  // full blown Bluestein algorithm...

  // modulate the input buffer (of size N) by the chirp-factors and store the results in the first
  // N elements of the internal y-buffer (of size M), thereby apply normalization if required:
  int n;
  if( normalizationFactor != 1.0 )
  {
    for(n=0; n<N; n++)
      y[n] = buffer[n] * c[n] *  normalizationFactor;
  }
  else
  {
    for(n=0; n<N; n++)
      y[n] = buffer[n] * c[n];
  }
  // ...and pad the rest of the y-buffer (elements N...M-1) with zeros:
  for(n=N; n<M; n++)
    y[n] = std::complex<T>(0.0, 0.0);

  // compute the radix-2 FFT of size M of the so defined y-buffer:
  transformerRadix2.setDirection(rsFourierTransformerRadix2<T>::FORWARD);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // compute the spectral product between the y-spectrum and the h-spectrum - we can overwrite the
  // y-buffer with these values instead of using a separate buffer:
  for(n=0; n<M; n++)
    y[n] = h[n] * y[n];

  // compute the inverse radix-2 FFT of size M of the so obtained y-buffer:
  transformerRadix2.setDirection(rsFourierTransformerRadix2<T>::INVERSE);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // the top N elements of y now contain the DFT divided by a chirp-factor of the input buffer -
  // multiply by the chirp-fcator and write them back into the input buffer:
  for(n=0; n<N; n++)
    buffer[n] = y[n] * c[n];

  // That's it! Bluestein-FFT done.
  // ToDo: check, that this works also for the inverse transform, i.e. arbitrary size inverse FFT
}

template<class T>
void rsFourierTransformerBluestein<T>::transformComplexBuffer(std::complex<T> *inBuffer,
  std::complex<T> *outBuffer)
{
  // copy the inBuffer into the outBuffer...
  for(int n=0; n<N; n++)
    outBuffer[n] = inBuffer[n];

  // ...and re-use the code in the in place routine on the outBuffer
  transformComplexBufferInPlace(outBuffer);
}

// pre-calculations:

template<class T>
void rsFourierTransformerBluestein<T>::generateChirp()
{
  if( h != NULL )
    delete[] h;
  h = new std::complex<T>[M];

  if( c != NULL )
    delete[] c;
  c = new std::complex<T>[N];

  T  theta = T(2.0*PI) / (T) N;
  std::complex<T> w_N;
  if( direction == rsFourierTransformerRadix2<T>::FORWARD )
    w_N = exp(std::complex<T>(0.0, -theta));
  else
    w_N = exp(std::complex<T>(0.0, theta));

  // compute the first N h-values and the corresponding c-values (their complex conjugates):
  int     n;
  std::complex<T> exponent;
  for(n = 0; n < N; n++)
  {
    exponent = std::complex<T>(T(-0.5*n*n), T(0));
    h[n]     = pow(w_N, exponent);
    c[n]     = conj(h[n]);
  }

  // cyclically wrap the h-sequence to fill up the remaining M-N values:
  for(n = 1; n < N; n++)
    h[M-n] = h[n];

  // obtain the discrete Fourier transform of the h-sequence by means of a radix-2 FFT of length M
  // - this can be done in place, because only the transformed h-values are needed in the rest of
  // the algorithm:
  transformerRadix2.setDirection(rsFourierTransformerRadix2<T>::FORWARD);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2<T>::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(h);
}

template<class T>
void rsFourierTransformerBluestein<T>::updateNormalizationFactor()
{
  normalizationFactor = rsFourierTransformerRadix2<T>::getNormalizationFactor(
    N, direction, normalizationMode);
}

/*


https://dsp.stackexchange.com/questions/24375/fastest-implementation-of-fft-in-c
https://github.com/project-gemmi/benchmarking-fft/

*/