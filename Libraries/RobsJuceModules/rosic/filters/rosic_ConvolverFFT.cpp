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


