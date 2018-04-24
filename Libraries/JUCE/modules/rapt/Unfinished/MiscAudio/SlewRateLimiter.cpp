// Construction/Destruction:

template<class TSig, class TPar>
rsSlewRateLimiter<TSig, TPar>::rsSlewRateLimiter()
{
  sampleRate  = 44100.0f; 
  attackTime  = 10.0f;    
  releaseTime = 100.0f;   
  y1          = 0.0;

  calculateAttackCoefficient();
  calculateReleaseCoefficient();
}

// Setup:

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    sampleRate = newSampleRate;
    calculateAttackCoefficient();
    calculateReleaseCoefficient();
  }
}

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::setAttackTime(TPar newAttackTime)
{
  if( newAttackTime >= 0.0 && newAttackTime != attackTime )
  {
    attackTime = newAttackTime; 
    calculateAttackCoefficient();
  }
}

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::setReleaseTime(TPar newReleaseTime)
{
  if( newReleaseTime >= 0.0 && newReleaseTime != releaseTime )
  {
    releaseTime = newReleaseTime; 
    calculateReleaseCoefficient();
  }
}

// Misc:

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::reset()
{
  y1 = 0.0;
}

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::calculateAttackCoefficient()
{
  double tau = 0.001*attackTime; // in seconds
  if( tau > 0.0 )
    coeffAttack = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffAttack = 0.0;
}

template<class TSig, class TPar>
void rsSlewRateLimiter<TSig, TPar>::calculateReleaseCoefficient()
{
  double tau = 0.001*releaseTime;
  if( tau > 0.0 )
    coeffRelease = exp( -1.0 / (sampleRate*tau)  );
  else
    coeffRelease = 0.0;
}
