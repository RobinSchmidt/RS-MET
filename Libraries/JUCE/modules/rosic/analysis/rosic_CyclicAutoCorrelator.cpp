#include "rosic_CyclicAutoCorrelator.h"
using namespace rosic;

// Construction/Destruction:

CyclicAutoCorrelator::CyclicAutoCorrelator(int newMaxBufferSize)
{
  if( newMaxBufferSize >= 256 )
    maxBufferSize = newMaxBufferSize;
  else
  {   
    maxBufferSize = 256;
    DEBUG_BREAK;  // we need some minimum buffer-size (the value is somewhat arbitrary)
  }

  // allocate memory for the buffer:
  buffer = NULL;
  buffer = new double[maxBufferSize];

  oldCycleLength = maxBufferSize;
  oldMeanSquare  = 0.0;

  // initialize the content of the buffer with zeros and reset the acumulators:
  resetAccumulators(); 
  resetBuffer(); 
}

CyclicAutoCorrelator::~CyclicAutoCorrelator()
{
  if( buffer != NULL )
    delete[] buffer;
}

// Miscellaneous:

void CyclicAutoCorrelator::resetAccumulators()
{
  sumOfProducts   = 0.0;
  sumOfSquares    = 0.0;
  sampleCounter   = 0;
}

void CyclicAutoCorrelator::resetBuffer()
{
  for(int i=0; i<maxBufferSize; i++)
    buffer[i] = 0.0;
}