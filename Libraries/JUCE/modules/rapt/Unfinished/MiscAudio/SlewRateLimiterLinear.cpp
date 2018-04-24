// Construction/Destruction:

template<class TSig, class TPar>
rsSlewRateLimiterLinear<TSig, TPar>::rsSlewRateLimiterLinear()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateUpwardLimit();
  calculateDownwardLimit();
}

// Setup:

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateUpwardLimit();
    calculateDownwardLimit();
  }
}

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::setAttackTime(TPar newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateUpwardLimit();
  }
}

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::setReleaseTime(TPar newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateDownwardLimit();
  }
}

// Misc:

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::reset()
{
  y1 = 0.0;
}

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::calculateUpwardLimit()
{
  upwardLimit = calculateStepLimit(attackTime);
}

template<class TSig, class TPar>
void rsSlewRateLimiterLinear<TSig, TPar>::calculateDownwardLimit()
{
  downwardLimit = calculateStepLimit(releaseTime);
}

template<class TSig, class TPar>
TPar rsSlewRateLimiterLinear<TSig, TPar>::calculateStepLimit(TPar unitStepTime)
{
  return 1.0 / (0.001*unitStepTime * sampleRate);
}