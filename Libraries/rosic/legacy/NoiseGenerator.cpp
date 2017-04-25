#include "NoiseGenerator.h"

//----------------------------------------------------------------------------
// construction/destruction:

NoiseGenerator::NoiseGenerator()
{
	setSampleRate(44100.0);  //sampleRate = 44100 Hz by default
 lowpassFilter.setMode(IirDesigner::LOWPASS);
 highpassFilter.setMode(IirDesigner::HIGHPASS);

 setLowpassOrder(2);
 setLowpassCutoff(5.0);

 setHighpassOrder(2);
 setHighpassCutoff(0.01);
}

NoiseGenerator::~NoiseGenerator()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void NoiseGenerator::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 lowpassFilter.setSampleRate(sampleRate);
 highpassFilter.setSampleRate(sampleRate);

 // for scaling and offsetting the range in order to have outputs bewteen 
 // -1.0...+1.0:
 scale  =  2.0 / (double) RAND_MAX;
 offset = -1.0;
}

void NoiseGenerator::setLowpassOrder(int newLowpassOrder)
{
 lowpassFilter.setSlope(newLowpassOrder);
}

void NoiseGenerator::setLowpassCutoff(double newLowpassCutoff)
{
 lowpassFilter.setFreq1(newLowpassCutoff);
 lowpassFilter.setGlobalGain(sqrt(sampleRate/newLowpassCutoff));
}

void NoiseGenerator::setHighpassOrder(int newHighpassOrder)
{
 highpassFilter.setSlope(newHighpassOrder);
}

void NoiseGenerator::setHighpassCutoff(double newHighpassCutoff)
{
 highpassFilter.setFreq1(newHighpassCutoff);
 highpassFilter.setGlobalGain(sqrt(sampleRate/(sampleRate-newHighpassCutoff)));
}

//----------------------------------------------------------------------------
//others:
