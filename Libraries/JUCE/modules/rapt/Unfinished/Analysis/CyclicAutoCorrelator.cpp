// Construction/Destruction:

template<class T>
rsCyclicAutoCorrelator<T>::rsCyclicAutoCorrelator(int newMaxBufferSize)
{
  if( newMaxBufferSize >= 256 )
    maxBufferSize = newMaxBufferSize;
  else
  {   
    maxBufferSize = 256;
    RS_DEBUG_BREAK;  // we need some minimum buffer-size (the value is somewhat arbitrary)
  }

  // allocate memory for the buffer:
  buffer = new T[maxBufferSize];

  oldCycleLength = maxBufferSize;
  oldMeanSquare  = 0.0;

  // initialize the content of the buffer with zeros and reset the acumulators:
  resetAccumulators(); 
  resetBuffer(); 
}

template<class T>
rsCyclicAutoCorrelator<T>::~rsCyclicAutoCorrelator()
{
  delete[] buffer;
}

// Misc:

template<class T>
void rsCyclicAutoCorrelator<T>::resetAccumulators()
{
  sumOfProducts   = 0.0;
  sumOfSquares    = 0.0;
  sampleCounter   = 0;
}

template<class T>
void rsCyclicAutoCorrelator<T>::resetBuffer()
{
  for(int i = 0; i < maxBufferSize; i++)
    buffer[i] = 0.0;
}
