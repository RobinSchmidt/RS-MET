#include "DemoVersionNoiseEmitter.h"

//----------------------------------------------------------------------------
// construction/destruction:

DemoVersionNoiseEmitter::DemoVersionNoiseEmitter()
{
 burstIntervalInSeconds = 120.0;
 sampleRate             = 44100.0;
 burstIntervalInSamples = (unsigned int) (burstIntervalInSeconds * sampleRate);
 burstScaler            = sampleRate/44100.0;
 sampleCounter          = 0;

 // for scaling and offsetting the range in order to have outputs bewteen 
 // -1.0...+1.0:
 scale  =  2.0 / (double) RAND_MAX;
 offset = -1.0;

 ampEnv.setSampleRate(sampleRate);

 // set up the amplitude-envelope for the burst:
 ampEnv.setStart(-300.0);  // value is in dB
 ampEnv.setAttack(1.0);
 ampEnv.setPeak(-24.0);      // also in dB
 ampEnv.setHold(1.0);      
 ampEnv.setDecay(1.0);  
 ampEnv.setSustain(-300.0);  
 ampEnv.setRelease(1.0);  
 ampEnv.setEnd(-300.0);  
}

DemoVersionNoiseEmitter::~DemoVersionNoiseEmitter()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void DemoVersionNoiseEmitter::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

 burstIntervalInSamples = (unsigned int) (burstIntervalInSeconds * sampleRate);
 burstScaler            = sampleRate/44100.0;
 ampEnv.setSampleRate(sampleRate);
}

void DemoVersionNoiseEmitter::setBurstInterval(double newBurstInterval)
{
 if( newBurstInterval >= 10.0 )
  burstIntervalInSeconds = newBurstInterval;

 burstIntervalInSamples = (unsigned int) (burstIntervalInSeconds * sampleRate);
}
