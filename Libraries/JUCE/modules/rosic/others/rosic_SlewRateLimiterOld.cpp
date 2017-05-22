//#include "rosic_SlewRateLimiter.h"
//using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

SlewRateLimiter::SlewRateLimiter()
{
	// parameter initializition:
 sampleRate  = 44100.0; // default sample rate
	attackTime  = 0.0;    // 10 ms attack time constant by default
	releaseTime = 0.0;    // 50 ms release time constant by default

	// calculation of internal variables:
	coeffAttack  = 0.0;
	coeffRelease = 0.0;

	// init previous output sample:
	y_1 = 0.0;
}

SlewRateLimiter::~SlewRateLimiter()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void SlewRateLimiter::setSampleRate(double newSampleRate)
{
	if(newSampleRate > 0.01)
		sampleRate = newSampleRate;

 setAttackTime(attackTime);
 setReleaseTime(releaseTime);

	//coeffAttack  = exp( -1.0 / (sampleRate*attackTime)  );
	//coeffRelease = exp( -1.0 / (sampleRate*releaseTime) );
}

void SlewRateLimiter::setAttackTime(double newAttackTime)
{
	if(newAttackTime > 0.0)   
 {
		attackTime  = newAttackTime;
	 coeffAttack = exp( -1.0 / (sampleRate*attackTime)  );
 }
 else
 {
  attackTime  = 0.0;
	 coeffAttack = 0.0;;
 }
}

void SlewRateLimiter::setReleaseTime(double newReleaseTime)
{
	if(newReleaseTime > 0.0)  
 {
		releaseTime  = newReleaseTime; 
	 coeffRelease = exp( -1.0 / (sampleRate*releaseTime) );
 }
 else
 {
		releaseTime  = 0.0; 
	 coeffRelease = 0.0;
 }
}

//----------------------------------------------------------------------------
// others:

void SlewRateLimiter::reset()
{
 y_1 = 0;
}