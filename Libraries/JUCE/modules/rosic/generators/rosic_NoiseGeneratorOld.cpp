//#include "rosic_NoiseGeneratorOld.h"
//using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

NoiseGeneratorOld::NoiseGeneratorOld()
{
  setSampleRate(44100.0);  //sampleRate = 44100 Hz by default
  lowpassFilter.setMode(FourPoleFilterParameters::LOWPASS_12);
  highpassFilter.setMode(FourPoleFilterParameters::HIGHPASS_12);

  setLowpassCutoff(5.0);
  setHighpassCutoff(0.01);
}

NoiseGeneratorOld::~NoiseGeneratorOld()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void NoiseGeneratorOld::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  lowpassFilter.setSampleRate(sampleRate);
  highpassFilter.setSampleRate(sampleRate);

  // for scaling and offsetting the range in order to have outputs bewteen 
  // -1.0...+1.0:
  scale  =  2.0 / (double) RAND_MAX;
  offset = -1.0;
}

/*
void NoiseGeneratorOld::setLowpassOrder(int newLowpassOrder)
{
  lowpassFilter.setSlope(newLowpassOrder);
}
*/

void NoiseGeneratorOld::setLowpassCutoff(double newLowpassCutoff)
{
  lowpassFilter.setFrequency(newLowpassCutoff);
  //lowpassFilter.setGlobalGain(sqrt(sampleRate/newLowpassCutoff));
}

/*
void NoiseGeneratorOld::setHighpassOrder(int newHighpassOrder)
{
  highpassFilter.setSlope(newHighpassOrder);
}
*/

void NoiseGeneratorOld::setHighpassCutoff(double newHighpassCutoff)
{
  highpassFilter.setFrequency(newHighpassCutoff);
  //highpassFilter.setGlobalGain(sqrt(sampleRate/(sampleRate-newHighpassCutoff)));
}

//----------------------------------------------------------------------------
//others:
