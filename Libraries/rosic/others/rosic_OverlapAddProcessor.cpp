#include "rosic_OverlapAddProcessor.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OverlapAddProcessor::OverlapAddProcessor(int maxBlockSize, int maxOverlapFactor, 
                                         int maxPaddingFactor)
{
  this->maxBlockSize     = maxBlockSize;
  this->maxOverlapFactor = maxOverlapFactor;
  this->maxPaddingFactor = maxPaddingFactor;
  blockSize              = rmin(1024, maxBlockSize);
  overlapFactor          = rmin(2, maxOverlapFactor);
  paddingFactor          = rmin(2, maxPaddingFactor);
  hopSize                = blockSize/overlapFactor;
  windowPower            = 2;
  useInputWindow         = true;
  useOutputWindow        = false;
  x                      = new double[maxBlockSize];
  w                      = new double[maxBlockSize];
  tmp                    = new double[maxBlockSize*maxPaddingFactor];
  readPositions          = new int[maxPaddingFactor*maxOverlapFactor];
  y                      = new double*[maxPaddingFactor*maxOverlapFactor];
  for(int i=0; i<maxPaddingFactor*maxOverlapFactor; i++)
    y[i] = new double[maxPaddingFactor*maxBlockSize];
  initInternalState();
}

OverlapAddProcessor::~OverlapAddProcessor()
{
  for(int i=0; i<maxPaddingFactor*maxOverlapFactor; i++)
    delete[] y[i];
  delete[] y;
  delete[] readPositions;
  delete[] tmp;
  delete[] w;
  delete[] x;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void OverlapAddProcessor::setInputBlockSize(int newSize)
{
  if( newSize != blockSize && newSize <= maxBlockSize )
  {
    blockSize = newSize;
    hopSize   = blockSize/overlapFactor;
    initInternalState();
  }
}

void OverlapAddProcessor::setOverlapFactor(int newFactor)
{
  if( newFactor != overlapFactor && newFactor <= maxOverlapFactor )
  {
    overlapFactor = newFactor;
    hopSize       = blockSize/overlapFactor;
    initInternalState();
  }
}

void OverlapAddProcessor::setPaddingFactor(int newFactor)
{
  if( newFactor != paddingFactor && newFactor <= maxPaddingFactor )
  {
    paddingFactor = newFactor;
    initInternalState();
  }
}

void OverlapAddProcessor::setWindowPower(int newPower)
{
  if( newPower >= 1 )
    windowPower = newPower;
  generateWindowFunction();
}

//-------------------------------------------------------------------------------------------------
// others:

void OverlapAddProcessor::clearBuffers()
{
  for(int n=0; n<maxBlockSize; n++)
    x[n] = 0.0;

  for(int n=0; n<maxBlockSize*maxPaddingFactor; n++)
    tmp[n] = 0.0;

  for(int i=0; i<maxPaddingFactor*maxOverlapFactor; i++)
  {
    for(int n=0; n<maxBlockSize*maxPaddingFactor; n++)
      y[i][n] = 0.0;
  }
}

//-------------------------------------------------------------------------------------------------
// internal fucntions:

void OverlapAddProcessor::prepareBlockForProcessing()
{
  // copy (a part of) the circular input buffer into the temporary processing buffer:
  int rp = writePosition - blockSize;  // read-position in circular buffer
  int wp = 0;                          // write-position in tmp-buffer
  if( rp < 0 )
    rp += maxBlockSize;
  for(wp=0; wp<blockSize; wp++)   
  {
    tmp[wp] = x[rp];
    rp++;
    if( rp >= maxBlockSize ) 
      rp = 0;
  }

  // apply window function, if desired:
  if( useInputWindow == true )
  {
    for(wp=0; wp<blockSize; wp++)   
      tmp[wp] *= w[wp];
  }

  // zero-pad tmp-buffer, if desired:
  for(wp=blockSize; wp<paddingFactor*blockSize; wp++)   
    tmp[wp] = 0.0;
}

void OverlapAddProcessor::postProcessBlock()
{
  // copy content of the temporary buffer into the appropriate output buffer:
  int n;
  for(n=0; n<paddingFactor*blockSize; n++)
    y[nextOutBuffer][n] = tmp[n];

  // apply output window, if desired:
  if( useOutputWindow == true )
  {
    for(n=0; n<blockSize; n++)   
      y[nextOutBuffer][n] *= w[n];
    for(n=blockSize; n<paddingFactor*blockSize; n++)   
      y[nextOutBuffer][n] = 0.0;
  }

  // reset the read-potsition in the output-buffer:
  readPositions[nextOutBuffer] = 0;

  // update the variable for the next output-buffer to write into:
  nextOutBuffer++;
  if( nextOutBuffer >= overlapFactor*paddingFactor )
    nextOutBuffer = 0;
}

void OverlapAddProcessor::initInternalState()
{
  nextOutBuffer = 0;
  initWritePointer();
  initReadPointers();
  clearBuffers();
  generateWindowFunction();
  calculateCompensationGain();
}

void OverlapAddProcessor::initWritePointer()
{
  writePosition = 0;
}

void OverlapAddProcessor::initReadPointers()
{
  for(int i=0; i<maxOverlapFactor*maxPaddingFactor; i++)
    readPositions[i] = i*hopSize;  // is this correct?
}

void OverlapAddProcessor::generateWindowFunction()
{
  makeCosinePowerWindow(w, blockSize, (double) windowPower);
}

void OverlapAddProcessor::calculateCompensationGain()
{
  if( useInputWindow == false && useOutputWindow == false )
    gain = 1.0 / overlapFactor;
  else
  {
    bool   squareWindow = useInputWindow == true && useOutputWindow == true;
    double accu         = 0.0;
    for(int i=0; i<overlapFactor; i++)
    {
      if( squareWindow )
        accu += w[i*hopSize] * w[i*hopSize];
      else
        accu += w[i*hopSize];
    }
    gain = 1.0 / accu;
  }
}