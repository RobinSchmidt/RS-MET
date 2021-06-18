//-------------------------------------------------------------------------------------------------
// construction/destruction:

//ConvolverPartitioned::ConvolverPartitioned()
//{
//  M = 0;
//}
//
//ConvolverPartitioned::~ConvolverPartitioned()
//{
//
//}

//-------------------------------------------------------------------------------------------------
// setup:

void ConvolverPartitioned::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  if( newLength < 0 ) { RAPT::rsError("Length must be >= 0"); return; }

  M = newLength;
  directConvolver.setImpulseResponse(newImpulseResponse, RAPT::rsMin(newLength,
    directConvolutionLength));

  int accu             = directConvolutionLength;
  int currentLength    = directConvolutionLength;
  int numFftConvolvers = 0;
  while(newLength > accu)
  {
    numFftConvolvers += 1;
    accu             += currentLength;
    currentLength    *= 2;
  }
  fftConvolvers.resize(numFftConvolvers);

  int currentStart = directConvolutionLength;
  currentLength    = directConvolutionLength;
  for(int c=0; c<numFftConvolvers; c++)
  {
    if(c == numFftConvolvers-1)
    {
      // Last block might be shorter than currentLength, so we pass a zero-padded version:
      double* finalBlock  = new double[currentLength];
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

//-------------------------------------------------------------------------------------------------
// others:

void ConvolverPartitioned::clearImpulseResponse()
{
  directConvolver.clearImpulseResponse();
  for(size_t c = 0; c < fftConvolvers.size(); c++)
    fftConvolvers[c].clearImpulseResponse();
}

void ConvolverPartitioned::clearInputBuffers()
{
  directConvolver.clearInputBuffer();
  for(size_t c = 0; c < fftConvolvers.size(); c++)
    fftConvolvers[c].clearInputBuffer();
}