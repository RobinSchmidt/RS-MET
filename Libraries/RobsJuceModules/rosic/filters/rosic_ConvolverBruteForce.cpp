//#include "rosic_ConvolverBruteForce.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ConvolverBruteForce::ConvolverBruteForce()
{
  mutex.lock();
  x = NULL;
  h = NULL;
  M = 0;
  mutex.unlock();
}

ConvolverBruteForce::~ConvolverBruteForce()
{
  mutex.lock();
  if( x != NULL ) delete[] x;
  if( h != NULL ) delete[] h;
  mutex.unlock();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void ConvolverBruteForce::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  mutex.lock();

  if( newLength < 0 )
  {
    DEBUG_BREAK;
    return;
  }

  if( newLength != M )
  {
    M = newLength;
    if( x != NULL ) delete[] x;
    if( h != NULL ) delete[] h;
    x = new double[M];
    h = new double[M];
    clearInputBuffer();
  }
  for(int k=0; k<M; k++)
    h[k] = newImpulseResponse[k];

  mutex.unlock();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void ConvolverBruteForce::clearImpulseResponse()
{
  mutex.lock();
  h[0] = 1.0;
  for(int k=1; k<M; k++)
    h[k] = 0.0;
  mutex.unlock();
}

void ConvolverBruteForce::clearInputBuffer()
{
  mutex.lock();
  for(int k=0; k<M; k++)
    x[k] = 0.0;
  mutex.unlock();
}


/*

ToDo:
-Get rid of the mutex stuff. Thread safety should not be the reponsibility of a DSP class.
-Maybe implement some sort of nonlinear generalization using instead of a sum ove h[k]*x[n-k] some 
 user-defined function f(h[k],x[n-k],k). May try something like f(x,h,k) = h*x + c*h*x^3

*/