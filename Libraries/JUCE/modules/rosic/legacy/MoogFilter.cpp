#include "MoogFilter.h"

//construction/destruction:
MoogFilter::MoogFilter()
{
 sampleRate      = 44100.0;
 sampleRateRec   = 1.0 / sampleRate;
 overSampling    = 8;
 overSamplingRec = 1.0 / (double)overSampling;

 outputTap       = 3;
	cutoff          = 20000.0;
	feedback        = 0.0;
	drive           = 1.0;

	setMode(1);     //mode 1 means moog lowpass

 teeBeeFactor  = 33.0/18.0;

 resetBuffers();

 //set up the anti-alias filter:
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
}

MoogFilter::~MoogFilter()
{

}

//----------------------------------------------------------------------------
//parameter settings:
void MoogFilter::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 sampleRateRec = 1.0 / sampleRate;
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
}


void MoogFilter::setOverSamp(int newOverSamp)
{
 /*
 if( (newOverSamp>=1) && (newOverSamp<=maxOverSampling) ) 
  overSampling = newOverSamp;

 overSamplingRec = 1.0 / (double)overSampling;
 antiAliasFlt.setSampleRate(overSampling*sampleRate);
 */
}

void MoogFilter::setMode(int newMode)
{
 mode = newMode;
}

void MoogFilter::setOutputTap(int newOutputTap)
{
 if( (newOutputTap>=1) && (newOutputTap<=4) )
  outputTap = newOutputTap-1;    //taps are indexed from 0-3 instead of 1-4
}

void MoogFilter::setDrive(double newDrive)
{
 drive = DB2AMP(newDrive);
}

//----------------------------------------------------------------------------
//others:
void MoogFilter::resetBuffers()
{
 for(long i=0; i<4; i++)
 {
  x_1[i] = 0;
  y_1[i] = 0;
 }

 previousIn  = 0.0;
 previousOut = 0.0;
}