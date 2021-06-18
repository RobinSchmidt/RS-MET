//-------------------------------------------------------------------------------------------------
// construction/destruction:

ConvolverPartitioned::ConvolverPartitioned()
{
  mutex.lock();
  M                = 0;
  fftConvolvers    = NULL;
  numFftConvolvers = 0;
  mutex.unlock();
}

ConvolverPartitioned::~ConvolverPartitioned()
{
  mutex.lock();
  if( fftConvolvers != NULL ) delete[] fftConvolvers;
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// setup:

void ConvolverPartitioned::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  mutex.lock();

  if( newLength < 0 )
  {
    DEBUG_BREAK;
    mutex.unlock();
    return;
  }

  if( newLength != M )  // this is wrong!
  {
    M = newLength;

    directConvolver.setImpulseResponse(newImpulseResponse, RAPT::rsMin(newLength, 
      directConvolutionLength));

    int accu          = directConvolutionLength;
    int currentLength = directConvolutionLength;
    numFftConvolvers  = 0;
    while( newLength > accu )
    {
      numFftConvolvers += 1;
      accu             += currentLength;
      currentLength    *= 2;
    }

    if( fftConvolvers != NULL )
      delete[] fftConvolvers;
    fftConvolvers = new ConvolverFFT[numFftConvolvers];

    int currentStart = directConvolutionLength;
    currentLength    = directConvolutionLength;
    for(int c=0; c<numFftConvolvers; c++)
    {
      if( c == numFftConvolvers-1 )
      {
        // Last block might be shorter than currentLength, so we pass a zero-padded version:
        double *finalBlock  = new double[currentLength];
        int     finalLength = newLength-currentStart;    // length of non-zero part
        int k;
        for(k=0; k<finalLength; k++)
          finalBlock[k] = newImpulseResponse[currentStart+k];
        for(k=finalLength; k<currentLength; k++)
          finalBlock[k] = 0.0;
        fftConvolvers[c].setImpulseResponse(finalBlock, currentLength);
        delete[] finalBlock;
        // ToDo: try to avoid the memory allocation, maybe setImpulseResponse should take an 
        // optional parameter "zeroPadding"
      }
      else
        fftConvolvers[c].setImpulseResponse(&(newImpulseResponse[currentStart]), currentLength);

      currentStart  += currentLength;
      currentLength *= 2;
    }
  }

  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// others:

void ConvolverPartitioned::clearImpulseResponse()
{
  mutex.lock();
  directConvolver.clearImpulseResponse();
  for(int c=0; c<numFftConvolvers; c++)
    fftConvolvers[c].clearImpulseResponse();
  mutex.unlock();
}

void ConvolverPartitioned::clearInputBuffers()
{
  mutex.lock();
  directConvolver.clearInputBuffer();
  for(int c=0; c<numFftConvolvers; c++)
    fftConvolvers[c].clearInputBuffer();
  mutex.unlock();
}