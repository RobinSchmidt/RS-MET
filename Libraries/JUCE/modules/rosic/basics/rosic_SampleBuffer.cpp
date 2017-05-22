//#include "rosic_SampleBuffer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SampleBuffer::SampleBuffer()
{
  mutex.lock();

  numChannels     = 0;
  numSamples      = 0;
  sampleRate      = 44100.0;
  sampleData      = NULL;
  channelPointers = NULL;

  mutex.unlock();
}

SampleBuffer::~SampleBuffer()
{
  freeMemory();
}

//-------------------------------------------------------------------------------------------------
// setup:

bool SampleBuffer::setSampleData(float **newSampleData, int newNumSamples, int newNumChannels, 
                                 float newSampleRate)
{
  if( newSampleData == NULL || newNumSamples == 0 || newNumChannels == 0 )
  {
    DEBUG_BREAK;     // nah - if NULL is passed, this should invalidate the data
    return false;
  }
  sampleRate = newSampleRate;

    
  mutex.lock();
  bool success = reAllocateMemoryIfNecessary(newNumSamples, newNumChannels);
  int  c, n;
  if( success == true )
  {
    // pointers have been sucesfully assigned to memory and members numSamples and numChannels have 
    // been updated - we now do the actual copying work, thereby repeat the first few samples at 
    // the end of each channel buffer for the interpolator:
    for(c=0; c<numChannels; c++)
    {
      for(n=0; n<numSamples; n++)
        channelPointers[c][n] = newSampleData[c][n];
      if( numSamples >= interpolatorMargin )
      {
        for(n=numSamples; n<numSamples+interpolatorMargin; n++)
          channelPointers[c][n] = newSampleData[c][n-numSamples];
      }
    }
  }
  mutex.unlock();
  return success;
}

void SampleBuffer::copyDataFrom(const SampleBuffer &source)
{
  /*
  // create deep copies of the member arrays (WARNING: this will create memory leaks when this 
  // object already has it's pointers assigned - i.e. the copyDataFrom was not called immediately
  // after construction):
  if( source.sampleName != NULL )
  {
    sampleNameLength     = source.sampleNameLength;
    sampleName           = new char[sampleNameLength+1];
    for(int c=0; c<=sampleNameLength; c++) // the <= is valid here, because we have one more cell allocated
      sampleName[c] = source.sampleName[c];
  }

  if( source.sampleData != NULL )
  {
    mutex.lock();
    sampleData = new double[numSamples+1];
    for(int n=0; n<numSamples; n++)
      sampleData[n] = source.sampleData[n];
    sampleData[numSamples] = sampleData[0];
    mutex.unlock();
  }
  */
}

//-------------------------------------------------------------------------------------------------
// inquiry and analysis:

double SampleBuffer::findNearestUpwardZeroCrossing(double position)
{
  mutex.lock();
  if( sampleData != NULL && numSamples > 0 && numChannels > 0)
    position = rosic::findNearestUpwardZeroCrossing(channelPointers[0], numSamples, position);
  mutex.unlock();
  return position;
}

/*
float* SampleBuffer::getChannelSampleData(int channel)
{
  if( channel > numChannels-1 )
    return NULL;
  else
    return channelPointers[channel];
}
*/

//-------------------------------------------------------------------------------------------------
// audio data manipulation:

void SampleBuffer::fillWithZeros()
{
  mutex.lock();
  for(int c=0; c<numChannels; c++)
  {
    for(int n=0; n<numSamples+interpolatorMargin; n++)
    {
      channelPointers[c][n] = 0.0;
    }
  }
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// others:

bool SampleBuffer::reAllocateMemoryIfNecessary(int newNumSamples, int newNumChannels)
{
  mutex.lock();

  if( newNumSamples == numSamples && newNumChannels == numChannels )
  {
    // no memory allocation to do - we can re-use the old memory area:
    fillWithZeros();
    mutex.unlock();
    return true;
  }

  // free old and allocate new memory:
  freeMemory();
  sampleData = new float[newNumChannels*(newNumSamples+interpolatorMargin)]; 
  if( sampleData == NULL )
  {
    mutex.unlock();
    return false;
  }
  channelPointers = new float*[newNumChannels];
  if( channelPointers == NULL )
  {
    mutex.unlock();
    return false;
  }

  // update members and initialize memory with zeros:
  numSamples  = newNumSamples;
  numChannels = newNumChannels;
  for(int c=0; c<numChannels; c++)
    channelPointers[c] = &( sampleData[c*(numSamples+interpolatorMargin)] );
  fillWithZeros();

  // memory allocation successfully done:   
  mutex.unlock();
  return true;
}

void SampleBuffer::freeMemory()
{
  mutex.lock();
  if( sampleData != NULL )
  {
    delete[] sampleData;
    sampleData = NULL;
    numSamples = 0;
  }
  if( channelPointers != NULL )
  {
    delete[] channelPointers;
    channelPointers = NULL;
    numChannels     = 0;
  }
  mutex.unlock();
}