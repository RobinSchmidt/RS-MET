#include "SampleLauncherVoice.h"

//----------------------------------------------------------------------------
// construction/destruction:

SampleLauncherVoice::SampleLauncherVoice()
{
 // init the pointers to the table with a NULL-pointer:
 sampleL = NULL;
 sampleR = NULL;    

 // init parameters:
 sampleLength = 0;
 position     = 0;
}

SampleLauncherVoice::~SampleLauncherVoice()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void SampleLauncherVoice::setSampleAddresses(float* newSampleAddressL, 
                                             float* newSampleAddressR)
{
 sampleL = newSampleAddressL;
 if( newSampleAddressR == NULL )
  sampleR = sampleL;
 else
  sampleR = newSampleAddressR;
}

void SampleLauncherVoice::setSampleLength(int newSampleLength)
{
 if( newSampleLength <= MAX_SAMPLE_LENGTH )
  sampleLength = newSampleLength;

 // reset the phase-pointer:
 reset();
}



//----------------------------------------------------------------------------
// others:

void SampleLauncherVoice::reset()
{
 position = 0;
}

