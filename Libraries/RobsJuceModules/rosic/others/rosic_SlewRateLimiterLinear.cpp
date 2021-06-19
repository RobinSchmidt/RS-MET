//#include "rosic_SlewRateLimiterLinear.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SlewRateLimiterLinear::SlewRateLimiterLinear()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateUpwardLimit();
  calculateDownwardLimit();
}

SlewRateLimiterLinear::~SlewRateLimiterLinear()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void SlewRateLimiterLinear::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateUpwardLimit();
    calculateDownwardLimit();
  }
}

void SlewRateLimiterLinear::setAttackTime(double newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateUpwardLimit();
  }
}

void SlewRateLimiterLinear::setReleaseTime(double newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateDownwardLimit();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SlewRateLimiterLinear::reset()
{
  y1 = 0.0;
}

void SlewRateLimiterLinear::calculateUpwardLimit()
{
  upwardLimit = calculateStepLimit(attackTime);
}

void SlewRateLimiterLinear::calculateDownwardLimit()
{
  downwardLimit = calculateStepLimit(releaseTime);
}

double SlewRateLimiterLinear::calculateStepLimit(double unitStepTime)
{
  // return 1.0 / (0.001*unitStepTime * sampleRate);   // reaches target one sample too early
  return 1.0 / ((0.001*unitStepTime * sampleRate)+1);  // yep, that seems to work
  // ToDo: let user set the time in samples, the formual is then:
  // 1.0 / (stepTimeInSamples+1)
}

double SlewRateLimiterLinear::getSample(double in) // inline this
{
  y1 += RAPT::rsClip(in-y1, -downwardLimit, upwardLimit);
  return y1;
  // This will be branchless, if the rsClip function is branchless (which it is not yet, but it 
  // can be made so using min/max)
}