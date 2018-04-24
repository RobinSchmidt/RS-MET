namespace RSLib
{

// Construction/Destruction:

rsSlewRateLimiterLinear::rsSlewRateLimiterLinear()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateUpwardLimit();
  calculateDownwardLimit();
}

rsSlewRateLimiterLinear::~rsSlewRateLimiterLinear()
{

}

// Setup:

void rsSlewRateLimiterLinear::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateUpwardLimit();
    calculateDownwardLimit();
  }
}

void rsSlewRateLimiterLinear::setAttackTime(double newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateUpwardLimit();
  }
}

void rsSlewRateLimiterLinear::setReleaseTime(double newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateDownwardLimit();
  }
}

// Misc:

void rsSlewRateLimiterLinear::reset()
{
  y1 = 0.0;
}

void rsSlewRateLimiterLinear::calculateUpwardLimit()
{
  upwardLimit = calculateStepLimit(attackTime);
}

void rsSlewRateLimiterLinear::calculateDownwardLimit()
{
  downwardLimit = calculateStepLimit(releaseTime);
}

double rsSlewRateLimiterLinear::calculateStepLimit(double unitStepTime)
{
  return 1.0 / (0.001*unitStepTime * sampleRate);
}

}