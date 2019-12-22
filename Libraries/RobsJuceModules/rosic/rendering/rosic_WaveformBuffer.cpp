//#include "rosic_WaveformBuffer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveformBuffer::WaveformBuffer()
{
  buffer      = NULL; 
  numSamples  = 0;
  numChannels = 0;   
}

WaveformBuffer::~WaveformBuffer()
{
  delete[] buffer;
}

//-------------------------------------------------------------------------------------------------
// setup:

void WaveformBuffer::setWaveform(double* newWaveForm, int newLength, rsString newName)
{
  allocateMemory(newLength);
  RAPT::rsArrayTools::copy(newWaveForm, buffer, numSamples);
  name = newName;
}

void WaveformBuffer::setWaveform(float* newWaveForm, int newLength, rsString newName)
{
  double* tmpBuffer = new double[newLength];
  for(int i=0; i<newLength; i++)
    tmpBuffer[i] = (double) newWaveForm[i];
  setWaveform(tmpBuffer, newLength, newName);
  delete[] tmpBuffer;
}
    
void WaveformBuffer::initWaveform(int newLength, rsString newName)
{
  allocateMemory(newLength);
  RAPT::rsArrayTools::fillWithZeros(buffer, numSamples);
  name = newName;
}

//-------------------------------------------------------------------------------------------------
// waveform retrieval:
    
void WaveformBuffer::getWaveform(double *targetBuffer, int length)
{
  if( buffer == NULL )
  {
    RAPT::rsArrayTools::fillWithZeros(targetBuffer, length);
    return;
  }
  if( length == numSamples )
    RAPT::rsArrayTools::copy(buffer, targetBuffer, numSamples);
  else
    RAPT::rsArrayTools::copyBufferWithLinearInterpolation(buffer, numSamples, targetBuffer, length);
}

//-------------------------------------------------------------------------------------------------
// others:
    
void WaveformBuffer::allocateMemory(int newNumSamples)
{
  if( newNumSamples != numSamples )
  {
    delete[] buffer;
    buffer = new double[newNumSamples];
    numSamples = newNumSamples;
  }
}