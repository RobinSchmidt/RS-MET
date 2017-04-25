#include "rosic_FourierTransformerBluestein.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FourierTransformerBluestein::FourierTransformerBluestein()
{
  h                     = NULL;
  c                     = NULL;
  y                     = NULL;
  N                     = 0;
  M                     = 0;
  direction             = FourierTransformerRadix2::FORWARD;
  normalizationMode     = FourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO;
  normalizationFactor   = 1.0;
  blockSizeIsPowerOfTwo = false;
}

FourierTransformerBluestein::~FourierTransformerBluestein()
{
  if( h != NULL )
    delete[] h;
  if( c != NULL )
    delete[] c;
  if( y != NULL )
    delete[] y;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void FourierTransformerBluestein::setBlockSize(int newBlockSize)
{
  // check new blocksize for validity and update our related member variables when it is valid (and 
  // different from the old one):
  if( newBlockSize > 0 && newBlockSize != N )
  {
    N = newBlockSize;

    // do the preprocessing for the Bluestein algorithm only when we really need it:
    blockSizeIsPowerOfTwo = isPowerOfTwo(N);
    if( !blockSizeIsPowerOfTwo )
    {
      M = nextPowerOfTwo(2*N-1);
      generateChirp();
    }
    else
      M = N;

    updateNormalizationFactor();

    // free old and allocate new memory for the internal buffer of size M:
    if( y != NULL )
      delete[] y;
    y = new Complex[M];
  }
}

void FourierTransformerBluestein::setDirection(int newDirection)
{
  if( newDirection >= FourierTransformerRadix2::FORWARD && 
      newDirection <= FourierTransformerRadix2::INVERSE )
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
    DEBUG_BREAK; // passed int-parameter does not correspond to any meaningful enum-field
}

void FourierTransformerBluestein::setNormalizationMode(int newNormalizationMode)
{
  if( newNormalizationMode >= FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO && 
      newNormalizationMode <= FourierTransformerRadix2::ORTHONORMAL_TRAFO )
  {
    normalizationMode = newNormalizationMode;
    updateNormalizationFactor();
  }
  else
    DEBUG_BREAK; // passed int-parameter does not correspond to any meaningful enum-field
}

//-------------------------------------------------------------------------------------------------
// signal processing:

void FourierTransformerBluestein::transformComplexBufferInPlace(Complex *buffer)
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
    y[n] = Complex(0.0, 0.0);

  // compute the radix-2 FFT of size M of the so defined y-buffer:
  transformerRadix2.setDirection(FourierTransformerRadix2::FORWARD);
  transformerRadix2.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // compute the spectral product between the y-spectrum and the h-spectrum - we can overwrite the
  // y-buffer with these values instead of using a separate buffer:
  for(n=0; n<M; n++)
    y[n] = h[n] * y[n];

  // compute the inverse radix-2 FFT of size M of the so obtained y-buffer:
  transformerRadix2.setDirection(FourierTransformerRadix2::INVERSE);
  transformerRadix2.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(y);

  // the top N elements of y now contain the DFT divided by a chirp-factor of the input buffer - 
  // multiply by the chirp-fcator and write them back into the input buffer:
  for(n=0; n<N; n++)
    buffer[n] = y[n] * c[n];

  // That's it! Bluestein-FFT done.
}

void FourierTransformerBluestein::transformComplexBuffer(Complex *inBuffer, Complex *outBuffer)
{
  // copy the inBuffer into the outBuffer...
  for(int n=0; n<N; n++)
    outBuffer[n] = inBuffer[n];

  // ...and re-use the code in the in place routine on the outBuffer
  transformComplexBufferInPlace(outBuffer);
}

//-------------------------------------------------------------------------------------------------
// pre-calculations:

void FourierTransformerBluestein::generateChirp()
{
  if( h != NULL )
    delete[] h;
  h = new Complex[M];

  if( c != NULL )
    delete[] c;
  c = new Complex[N];

  double  theta = 2.0 * PI / (double) N;
  Complex w_N;
  if( direction == FourierTransformerRadix2::FORWARD )
    w_N = expC(Complex(0.0, -theta));
  else
    w_N = expC(Complex(0.0, theta));

  // compute the first N h-values and the corresponding c-values (their complex conjugates):
  int     n;  
  Complex exponent;
  for(n = 0; n < N; n++) 
  {
    exponent = Complex(-0.5*n*n, 0.0);
    h[n]     = powC(w_N, exponent);
    c[n]     = h[n].getConjugate();
  }

  // cyclically wrap the h-sequence to fill up the remaining M-N values:
  for(n = 1; n < N; n++)
    h[M-n] = h[n];

  // obtain the discrete Fourier transform of the h-sequence by means of a radix-2 FFT of length M
  // - this can be done in place, because only the transformed h-values are needed in the rest of 
  // the algorithm:
  transformerRadix2.setDirection(FourierTransformerRadix2::FORWARD);
  transformerRadix2.setNormalizationMode(FourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO);
  transformerRadix2.setBlockSize(M);
  transformerRadix2.transformComplexBufferInPlace(h);
}

void FourierTransformerBluestein::updateNormalizationFactor()
{
  if( (normalizationMode == FourierTransformerRadix2::NORMALIZE_ON_FORWARD_TRAFO && 
       direction == FourierTransformerRadix2::FORWARD) ||
      (normalizationMode == FourierTransformerRadix2::NORMALIZE_ON_INVERSE_TRAFO && 
       direction == FourierTransformerRadix2::INVERSE)    )
  {
    normalizationFactor = 1.0 / (double) N;
  }
  else if( normalizationMode == FourierTransformerRadix2::ORTHONORMAL_TRAFO )
  {
    normalizationFactor = 1.0 / sqrt((double) N);
  }
  else
    normalizationFactor = 1.0;
}



