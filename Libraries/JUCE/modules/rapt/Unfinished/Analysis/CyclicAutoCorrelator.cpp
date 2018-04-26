using namespace RSLib;

// Construction/Destruction:

rsCyclicAutoCorrelator::rsCyclicAutoCorrelator(int newMaxBufferSize)
{
  if( newMaxBufferSize >= 256 )
    maxBufferSize = newMaxBufferSize;
  else
  {   
    maxBufferSize = 256;
    RS_DEBUG_BREAK;  // we need some minimum buffer-size (the value is somewhat arbitrary)
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

rsCyclicAutoCorrelator::~rsCyclicAutoCorrelator()
{
  if( buffer != NULL )
    delete[] buffer;
}

// Misc:

void rsCyclicAutoCorrelator::resetAccumulators()
{
  sumOfProducts   = 0.0;
  sumOfSquares    = 0.0;
  sampleCounter   = 0;
}

void rsCyclicAutoCorrelator::resetBuffer()
{
  for(int i = 0; i < maxBufferSize; i++)
    buffer[i] = 0.0;
}
