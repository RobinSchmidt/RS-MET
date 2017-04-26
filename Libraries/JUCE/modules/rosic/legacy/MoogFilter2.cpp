#include "MoogFilter2.hpp"

//construction/destruction:
MoogFilter2::MoogFilter2()
{
 sampleRate      = 44100.0;
 sampleRateRec   = 1.0/sampleRate;
 overSampling    = 10;
 overSamplingRec = 1.0/(sample)overSampling;

 outputTap       = 3;
	cutoff          = 20000.0;
	feedback        = 0.0;
	drive           = 1.0;

 V_t             = 1.0;
 beta            = 1.0/(2*V_t);

 //set up the anti-alias filter
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
}

MoogFilter2::~MoogFilter2()
{}
//------------------------------------------------------------------------------------------------------------
//parameter settings:
void MoogFilter2::setSampleRate(sample SampleRate)
{
 if(SampleRate>0.0)
  sampleRate = SampleRate;
 sampleRateRec = 1/sampleRate;
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
}
void MoogFilter2::setOverSamp(long OverSamp)
{
 if(OverSamp>=1 && OverSamp<=maxOverSampling) 
  overSampling = OverSamp;
 overSamplingRec = 1.0/(sample)overSampling;
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
}

void MoogFilter2::setOutputTap(long OutputTap)
{
 if(OutputTap>=1 && OutputTap<=4)
  outputTap = OutputTap-1;    //taps are indexed from 0-3 instead of 1-4
}

void MoogFilter2::setDrive(sample Drive)
{
 drive = DB2AMP(Drive);
}

//------------------------------------------------------------------------------------------------------------
//audio processing:


//------------------------------------------------------------------------------------------------------------
//others:
void MoogFilter2::resetBuffers()
{
 uint32 i;
 for(i=0; i<4; i++)
  y_1[i] = 0;
 for(i=0; i<3; i++)
  W_1[i] = 0;
 prevInSamp  = 0.0;
 prevOutSamp = 0.0;
}