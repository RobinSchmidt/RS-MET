using namespace RSLib;

// Construction/Destruction:

rsSlewRateLimiter::rsSlewRateLimiter()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateAttackCoefficient();
  calculateReleaseCoefficient();
}

rsSlewRateLimiter::~rsSlewRateLimiter()
{

}

// Setup:

void rsSlewRateLimiter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateAttackCoefficient();
    calculateReleaseCoefficient();
  }
}

void rsSlewRateLimiter::setAttackTime(double newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateAttackCoefficient();
  }
}

void rsSlewRateLimiter::setReleaseTime(double newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateReleaseCoefficient();
  }
}

// Misc:

void rsSlewRateLimiter::reset()
{
  y1 = 0.0;
}

void rsSlewRateLimiter::calculateAttackCoefficient()
{
  double tau = 0.001*attackTime; // in seconds
  if( tau > 0.0 )
    coeffAttack = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffAttack = 0.0;
}

void rsSlewRateLimiter::calculateReleaseCoefficient()
{
  double tau = 0.001*releaseTime;
  if( tau > 0.0 )
    coeffRelease = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffRelease = 0.0;
}
