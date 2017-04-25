#include "rosic_SlewRateLimiter.h"
using namespace rosic;

// Construction/Destruction:

SlewRateLimiter::SlewRateLimiter()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateAttackCoefficient();
  calculateReleaseCoefficient();
}

SlewRateLimiter::~SlewRateLimiter()
{

}

// Setup:

void SlewRateLimiter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateAttackCoefficient();
    calculateReleaseCoefficient();
  }
}

void SlewRateLimiter::setAttackTime(double newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateAttackCoefficient();
  }
}

void SlewRateLimiter::setReleaseTime(double newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateReleaseCoefficient();
  }
}

// Miscellaneous:

void SlewRateLimiter::reset()
{
  y1 = 0.0;
}

void SlewRateLimiter::calculateAttackCoefficient()
{
  double tau = 0.001*attackTime; // in seconds
  if( tau > 0.0 )
    coeffAttack = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffAttack = 0.0;
}

void SlewRateLimiter::calculateReleaseCoefficient()
{
  double tau = 0.001*releaseTime;
  if( tau > 0.0 )
    coeffRelease = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffRelease = 0.0;
}
