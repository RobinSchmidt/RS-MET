#include "OouraFFT/fft4g.c"
using namespace RSLib;

// construction/destruction:

rsFourierTransformerRadix2::rsFourierTransformerRadix2()
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

rsFourierTransformerRadix2::~rsFourierTransformerRadix2()
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

void rsFourierTransformerRadix2::setBlockSize(int newBlockSize)
{
  // check new blocksize for validity:
  if( newBlockSize >= 2 && rsIsPowerOfTwo(newBlockSize) )
  {
    // check, if the new blocksize is actually different from the old one in order to avoid 
    // unnecesarry re-allocations and re-computations:
    if( newBlockSize != N )
    {
      N    = newBlockSize;
      logN = (int) floor( rsLog2((double) N + 0.5 ) );
      updateNormalizationFactor();

      if( w != NULL )
        delete[] w;
      w    = new double[2*N];

      if( ip != NULL )
        delete[] ip;
      ip    = new int[(int) ceil(4.0+rsSqrt((double)N))];
      ip[0] = 0; // indicate that re-initialization is necesarry

      if( tmpBuffer != NULL )
        delete[] tmpBuffer;
      tmpBuffer = new rsComplexDbl[N];
    }
  }
  else if( !rsIsPowerOfTwo(newBlockSize) || newBlockSize <= 1 )
    RS_DEBUG_BREAK; // this class can only deal with blocksizes >= 2 that are a power of two
}

void rsFourierTransformerRadix2::setDirection(int newDirection)
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

void rsFourierTransformerRadix2::setNormalizationMode(int newNormalizationMode)
{
  if( newNormalizationMode >= NORMALIZE_ON_FORWARD_TRAFO && 
      newNormalizationMode <= ORTHONORMAL_TRAFO )
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

void rsFourierTransformerRadix2::transformComplexBufferInPlace(rsComplexDbl *buffer)
{
  // retrieve the adresses of the real part of the first array entries in order to treat the 
  // Complex arrays as arrays of two successive double-numbers:
  double* d_buffer = &(buffer[0].re);

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

void rsFourierTransformerRadix2::transformComplexBuffer(rsComplexDbl *inBuffer, 
                                                        rsComplexDbl *outBuffer)
{
  // retrieve the adresses of the real part of the first array entries in order to treat the 
  // Complex arrays as arrays of two successive double-numbers:
  double* d_inBuffer  = &(inBuffer[0].re);
  double* d_outBuffer = &(outBuffer[0].re);

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

void rsFourierTransformerRadix2::transformRealSignal(double *inSignal, rsComplexDbl *outSpectrum)
{
  setDirection(FORWARD);

  // retrieve the adress of the real part of the first array entry of the output array in order to
  // treat the Complex array as array of two successive double-numbers:
  double* d_outBuffer = &(outSpectrum[0].re);

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

void rsFourierTransformerRadix2::transformRealSignal(double *signal, double *reAndIm)
{
  rsComplexDbl* c_reAndIm = (rsComplexDbl*) &(reAndIm[0]);
  transformRealSignal(signal, c_reAndIm);
}

void rsFourierTransformerRadix2::getRealSignalMagnitudesAndPhases(double *signal, 
                                                                  double *magnitudes, 
                                                                  double *phases)
{
  transformRealSignal(signal, tmpBuffer);

  // store the two purely real transform values at DC and Nyquist-frequency in the first fields of 
  // the magnitude- and phase- arrays respectively:
  magnitudes[0] = tmpBuffer[0].re;
  phases[0]     = tmpBuffer[0].im;

  // fill the rest of the array with the magnitudes and phases of the regular bins:
  double* dBuffer = &(tmpBuffer[0].re);
  double  re, im;
  int     k;
  for(k=1; k<N/2; k++)
  {
    re            = dBuffer[2*k];
    im            = dBuffer[2*k+1];
    magnitudes[k] = rsSqrt(re*re + im*im);
    if( re == 0.0 && im == 0.0 )
      phases[k] = 0.0;
    else
      phases[k] = atan2(im, re);
  }
}

void rsFourierTransformerRadix2::getRealSignalMagnitudes(double *signal, double *magnitudes)
{
  transformRealSignal(signal, tmpBuffer);
  magnitudes[0] = tmpBuffer[0].re;

  double* dBuffer = &(tmpBuffer[0].re);
  double  re, im;
  int     k;
  for(k=1; k<N/2; k++)
  {
    re            = dBuffer[2*k];
    im            = dBuffer[2*k+1];
    magnitudes[k] = rsSqrt(re*re + im*im);
  }
}

void rsFourierTransformerRadix2::transformSymmetricSpectrum(rsComplexDbl *inSpectrum, 
                                                            double *outSignal)
{
  setDirection(INVERSE);

  // retrieve the adress of the real part of the first array entry of the output array in order to
  // treat the Complex array as array of two successive double-numbers:
  double* d_inBuffer = &(inSpectrum[0].re);

  // copy the input into the output for the in-place routine (thereby normalize, if necesarry):
  int n;
  if( normalizationFactor != 1.0 )
  {
    for(n=0; n<N; n++)
      outSignal[n] = 2.0 * d_inBuffer[n] * normalizationFactor;
  }
  else
  {
    for(n=0; n<N; n++)
      outSignal[n] = 2.0 * d_inBuffer[n];
  }

  // for some reason, the subsequent routine expects the second half of the spectrum (the complex 
  // conjugate values of the first half), so we need to take the complex conjugates:
  for(n=3; n<N; n+=2) // start at n=3 (imaginary part of the first bin after DC)
    outSignal[n] = -outSignal[n];

  // use Ooura's routine:
  rdft(N, -1, outSignal, ip, w);
}

void rsFourierTransformerRadix2::transformSymmetricSpectrum(double *reAndIm, double *signal)
{
  rsComplexDbl* c_reAndIm = (rsComplexDbl*) &(reAndIm[0]);
  transformSymmetricSpectrum(c_reAndIm, signal);
}

void rsFourierTransformerRadix2::getRealSignalFromMagnitudesAndPhases(double *magnitudes, 
                                                                      double *phases, 
                                                                      double *signal)
{
  tmpBuffer[0].re = magnitudes[0];
  tmpBuffer[0].im = phases[0];

  int k;
  double* dBuffer = &(tmpBuffer[0].re);
  double  s, c;
  for(k=1; k<N/2; k++)
  {
    rsSinCos(phases[k], &s, &c);
    dBuffer[2*k]   = magnitudes[k] * c;
    dBuffer[2*k+1] = magnitudes[k] * s;
  }

  transformSymmetricSpectrum(tmpBuffer, signal);
}

// pre-calculations:

void rsFourierTransformerRadix2::updateNormalizationFactor()
{
  if( (normalizationMode == NORMALIZE_ON_FORWARD_TRAFO && direction == FORWARD) ||
      (normalizationMode == NORMALIZE_ON_INVERSE_TRAFO && direction == INVERSE)    )
  {
    normalizationFactor = 1.0 / (double) N;
  }
  else if( normalizationMode == ORTHONORMAL_TRAFO )
  {
    normalizationFactor = 1.0 / rsSqrt((double) N);
  }
  else
    normalizationFactor = 1.0;
}

//=================================================================================================

// construction/destruction:

rsFourierTransformerBluestein::rsFourierTransformerBluestein()
{
  h                     = NULL;
  c                     = NULL;
  y                     = NULL;
  N                     = 0;
  M                     = 0;
  direction             = rsFourierTransformerRadix2::FORWARD;
  normalizationMode     = rsFourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO;
  normalizationFactor   = 1.0;
  blockSizeIsPowerOfTwo = false;
}

rsFourierTransformerBluestein::~rsFourierTransformerBluestein()
{
  if( h != NULL )
    delete[] h;
  if( c != NULL )
    delete[] c;
  if( y != NULL )
    delete[] y;
}

// parameter settings:

void rsFourierTransformerBluestein::setBlockSize(int newBlockSize)
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
    y = new rsComplexDbl[M];
  }
}

void rsFourierTransformerBluestein::setDirection(int newDirection)
{
  if( newDirection >= rsFourierTransformerRadix2::FORWARD && 
      newDirection <= rsFourierTransformerRadix2::INVERSE )
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

void rsFourierTransformerBluestein::setNormalizationMode(int newNormalizationMode)
{
  if( newNormalizationMode >= rsFourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO && 
      newNormalizationMode <= rsFourierTransformerRadix2::ORTHONORMAL_TRAFO )
  {
    normalizationMode = newNormalizationMode;
    updateNormalizationFactor();
  }
  else
    rsError("passed int-parameter does not correspond to any meaningful enum-field");
}

// signal processing:

void rsFourierTransformerBluestein::transformComplexBufferInPlace(rsComplexDbl *buffer)
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
    y[n] = rsComplexDbl(0.0, 0.0);

  // compute the radix-2 FFT of size M of the so defined y-buffer:
  transformerRadix2.setDirection(rsFourierTransformerRadix2::FORWARD);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // compute the spectral product between the y-spectrum and the h-spectrum - we can overwrite the
  // y-buffer with these values instead of using a separate buffer:
  for(n=0; n<M; n++)
    y[n] = h[n] * y[n];

  // compute the inverse radix-2 FFT of size M of the so obtained y-buffer:
  transformerRadix2.setDirection(rsFourierTransformerRadix2::INVERSE);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // the top N elements of y now contain the DFT divided by a chirp-factor of the input buffer - 
  // multiply by the chirp-fcator and write them back into the input buffer:
  for(n=0; n<N; n++)
    buffer[n] = y[n] * c[n];

  // That's it! Bluestein-FFT done.
}

void rsFourierTransformerBluestein::transformComplexBuffer(rsComplexDbl *inBuffer, 
  rsComplexDbl *outBuffer)
{
  // copy the inBuffer into the outBuffer...
  for(int n=0; n<N; n++)
    outBuffer[n] = inBuffer[n];

  // ...and re-use the code in the in place routine on the outBuffer
  transformComplexBufferInPlace(outBuffer);
}

// pre-calculations:

void rsFourierTransformerBluestein::generateChirp()
{
  if( h != NULL )
    delete[] h;
  h = new rsComplexDbl[M];

  if( c != NULL )
    delete[] c;
  c = new rsComplexDbl[N];

  double  theta = 2.0 * PI / (double) N;
  rsComplexDbl w_N;
  if( direction == rsFourierTransformerRadix2::FORWARD )
    w_N = rsExpC(rsComplexDbl(0.0, -theta));
  else
    w_N = rsExpC(rsComplexDbl(0.0, theta));

  // compute the first N h-values and the corresponding c-values (their complex conjugates):
  int     n;  
  rsComplexDbl exponent;
  for(n = 0; n < N; n++) 
  {
    exponent = rsComplexDbl(-0.5*n*n, 0.0);
    h[n]     = rsPowC(w_N, exponent);
    c[n]     = h[n].getConjugate();
  }

  // cyclically wrap the h-sequence to fill up the remaining M-N values:
  for(n = 1; n < N; n++)
    h[M-n] = h[n];

  // obtain the discrete Fourier transform of the h-sequence by means of a radix-2 FFT of length M
  // - this can be done in place, because only the transformed h-values are needed in the rest of 
  // the algorithm:
  transformerRadix2.setDirection(rsFourierTransformerRadix2::FORWARD);
  transformerRadix2.setNormalizationMode(rsFourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(h);
}

void rsFourierTransformerBluestein::updateNormalizationFactor()
{
  if( (normalizationMode == rsFourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO && 
       direction == rsFourierTransformerRadix2::FORWARD) ||
      (normalizationMode == rsFourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO && 
       direction == rsFourierTransformerRadix2::INVERSE)    )
  {
    normalizationFactor = 1.0 / (double) N;
  }
  else if( normalizationMode == rsFourierTransformerRadix2::ORTHONORMAL_TRAFO )
  {
    normalizationFactor = 1.0 / sqrt((double) N);
  }
  else
    normalizationFactor = 1.0;
}

