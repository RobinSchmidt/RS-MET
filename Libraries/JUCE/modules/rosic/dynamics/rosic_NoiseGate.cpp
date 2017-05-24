//#include "rosic_NoiseGate.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

NoiseGate::NoiseGate(int newLookAheadBufferSize) : DynamicsProcessorBase(newLookAheadBufferSize)
{
  threshold      =  -60.0;
  hysteresis     =    6.0;
  holdTime       =  10.0;
  gateClosed     = true;
  numHoldSamples = roundToInt(0.001*holdTime*sampleRate);
  holdCounter    = 0;
  updateThresholds();
  levelDetector.setMode(LevelDetector::MEAN_ABS);
  levelDetector.setAttackTime(0.0);
  setAttackTime(0.0);    // affects only the AR/follower 
  setReleaseTime(10.0);  // affects level-detector and AR/follower
}

NoiseGate::~NoiseGate()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void NoiseGate::setSampleRate(double newSampleRate)    
{ 
  DynamicsProcessorBase::setSampleRate(newSampleRate);
  numHoldSamples = roundToInt(0.001*holdTime*sampleRate);
}

void NoiseGate::setThreshold(double newThreshold)
{
  threshold = newThreshold;
  updateThresholds();
}

void NoiseGate::setHysteresis(double newHysteresis)
{
  hysteresis = newHysteresis;
  updateThresholds();
}

void NoiseGate::setHoldTime(double newHoldTime)
{
  holdTime       = newHoldTime;
  numHoldSamples = roundToInt(0.001*holdTime*sampleRate);
}

//-------------------------------------------------------------------------------------------------
// interanal functions:

void NoiseGate::updateThresholds()
{
  openingThreshold = dB2amp(threshold);
  closingThreshold = dB2amp(threshold-hysteresis);
}

