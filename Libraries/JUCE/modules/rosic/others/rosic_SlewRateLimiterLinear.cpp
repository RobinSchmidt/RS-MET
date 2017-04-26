#include "rosic_SlewRateLimiterLinear.h"
using namespace rosic;

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
  return 1.0 / (0.001*unitStepTime * sampleRate);
}

double SlewRateLimiterLinear::getSample(double in) // inline this
{
  y1 += clip(in-y1, -downwardLimit, upwardLimit);
  return y1;
}