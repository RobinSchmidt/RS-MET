//#include "romos_ProcessingStatus.h"
//using namespace romos;

ProcessingStatus processingStatus;  // definition of the global object

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ProcessingStatus::ProcessingStatus()
{
  setSystemSampleRate(44100.0);
  setTempo(120.0);
  setBufferSize(66);  // 66 yields best performance with IdentityChain performance test
}

ProcessingStatus::~ProcessingStatus()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void ProcessingStatus::setSystemSampleRate(double newSampleRate)
{
  systemSampleRate   = newSampleRate;
  systemSamplePeriod = 1.0      / systemSampleRate;
  freqToOmegaFactor  = (2.0*PI) / systemSampleRate;
}

void ProcessingStatus::setTempo(double newTempo)
{
  tempo = newTempo;
}

void ProcessingStatus::setBufferSize(int newBufferSize)
{
  if( newBufferSize <= maxBufferSize )
    bufferSize = newBufferSize;
}

