//#include "rosic_ConvolverFFT.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ConvolverFFT::ConvolverFFT()
{
  x             = NULL;
  y1            = NULL;
  y2            = NULL;
  h             = NULL;
  X             = NULL;
  H             = NULL;
  L             = 0;
  M             = 0;
  writeCounter  = 0;
  readCounter1  = 0;
  readCounter2  = 0;
  useOutBuffer2 = false;

  forwardTransformer.setDirection(FourierTransformerRadix2::FORWARD);
  inverseTransformer.setDirection(FourierTransformerRadix2::INVERSE);
}

ConvolverFFT::~ConvolverFFT()
{
  if( x  != NULL ) delete[] x;
  if( y1 != NULL ) delete[] y1;
  if( y2 != NULL ) delete[] y2;
  if( h  != NULL ) delete[] h;
  if( X  != NULL ) delete[] X;
  if( H  != NULL ) delete[] H;
}
// remove the ifs - it's ok to delete nulltprs!

//-------------------------------------------------------------------------------------------------
// parameter settings:

void ConvolverFFT::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  if( newLength < 0 )
  {
    DEBUG_BREAK;
    return;
  }

  // re-allocate memory (if necessarry):
  allocateBuffers(newLength);   // will also update members M and L

  // write the new impulse response into an internal buffer, thereby zero-pad it:
  int k;
  for(k=0; k<L; k++)
    h[k] = newImpulseResponse[k];
  for(k=L; k<M; k++)
    h[k] = 0.0;

  // transform the impulse response into the FFT domain:
  forwardTransformer.transformRealSignal(h, H);
}

//-------------------------------------------------------------------------------------------------
// others:

void ConvolverFFT::clearImpulseResponse()
{
  h[0] = 1.0;
  for(int k=1; k<M; k++)
    h[k] = 0.0;
  forwardTransformer.transformRealSignal(h, H);
}

void ConvolverFFT::clearInputBuffer()
{
  for(int k=0; k<M; k++)
    x[k] = 0.0;
}

void ConvolverFFT::clearOutputBuffer1()
{
  for(int k=0; k<M; k++)
    y1[k] = 0.0;
}

void ConvolverFFT::clearOutputBuffer2()
{
  for(int k=0; k<M; k++)
    y2[k] = 0.0;
}

void ConvolverFFT::clearBuffers()
{
  clearInputBuffer();
  clearOutputBuffer1();
  clearOutputBuffer2();
}

void ConvolverFFT::allocateBuffers(int newImpulseResponseLength)
{
  if( newImpulseResponseLength != L )
  {
    L = newImpulseResponseLength;
    M = RAPT::rsNextPowerOfTwo(2*L);

    if( x  != NULL ) delete[] x;
    if( y1 != NULL ) delete[] y1;
    if( y2 != NULL ) delete[] y2;
    if( h  != NULL ) delete[] h;
    if( X  != NULL ) delete[] X;
    if( H  != NULL ) delete[] H;
    // remove the ifs - it's ok to delete nulltprs!

    x  = new double[M];
    y1 = new double[M];
    y2 = new double[M];
    h  = new double[M];
    X  = new Complex[M/2];
    H  = new Complex[M/2];

    writeCounter  = 0;
    readCounter1  = M/2;
    readCounter2  = 0;
    useOutBuffer2 = false;

    forwardTransformer.setBlockSize(M);
    inverseTransformer.setBlockSize(M);

    clearBuffers();
  }
}


/*


Ideas:

NTT Convolver:
-Implement a convolver based on a number theorectic transform rsConvolverNTT
-Figure out if modular arithemtic is faster than complex. Certainly, if the modulus is a power of 2 
 such the we may use bitmasking for the modulo operation...but for NTT, we need the modulus to be a 
 prime (i think)
-Maybe we do not need to take the remainder after each multiplcation but only one at the very
 end (or whenever there's an risk of overflow)
-Maybe we can assume that the number is always <= 2*m (m being the modulus), if that's the case, we 
 can replace the mod operation by: if(x >= m) x -= m; That could also be done branchless like:
 x = (1-c)*x + c*(x-m)  where  c = (x < m), complexities:
          mul                           add
 complex: 4 mul, 1 add, 1 sub           1 add
 modular: 3 mul, 1 cmp, 1 add, 1 sub    2 add, 1 cmp, 2 mul, 2 sub
     or:  1 mul, 1 mod                  1 add, 1 mod
-Map the float range -1.0..+1.0 to a range large enough to represent 24 bit integers...and maybe 
 spend 2 or 3 extra bits..Let's say, we use the range of 26 bit integers, then multipying 2 of them 
 needs at most 52 for the resul, so with 64 bit integers, we should have enough headroom to avoid 
 overflow. ..Oh but the overflow occurs already at the modulus, Maybe choose the largest prime 
 below 2^64
-Maybe we can try modular arithmetic based on a power of 2 and do radix-3 NTTs?

FHT Convolver:
-Implement convolution based on the fast discrete Hartley transform -> avoids complex arithmetic:
 https://en.wikipedia.org/wiki/Discrete_Hartley_transform
 RADIX-2 FAST HARTLEY TRANSFORM REVISITED: https://arxiv.org/ftp/arxiv/papers/1503/1503.03794.pdf
 https://www.researchgate.net/publication/291912452_FAST_ALGORITHMS_FOR_THE_DISCRETE_HARTLEY_TRANSFORM

*/

