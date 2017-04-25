#include "SlewRateLimiter.h"

//----------------------------------------------------------------------------
// construction/destruction:

SlewRateLimiter::SlewRateLimiter()
{
	// parameter initializition:
 sampleRate  = 44100.0; // default sample rate
	attackTime  = 0.01;    // 10 ms attack time constant by default
	releaseTime = 0.05;    // 50 ms release time constant by default

	// calculation of internal variables:
	coeffAttack  = exp( -1.0 / (sampleRate*attackTime)  );
	coeffRelease = exp( -1.0 / (sampleRate*releaseTime) );

	// init previous output sample:
	y_1 = 0.0;
}

SlewRateLimiter::~SlewRateLimiter()
{

}

//----------------------------------------------------------------------------
//parameter settings:

void SlewRateLimiter::setSampleRate(double newSampleRate)
{
	if(newSampleRate > 0.01)
		sampleRate = newSampleRate;

	coeffAttack  = exp( -1.0 / (sampleRate*attackTime)  );
	coeffRelease = exp( -1.0 / (sampleRate*releaseTime) );
}

void SlewRateLimiter::setAttackTime(double newAttackTime)
{
	if(newAttackTime > 0.0001)    // 1/10 of a millisecond minimum
		attackTime = newAttackTime;

	coeffAttack  = exp( -1.0 / (sampleRate*attackTime)  );
}

void SlewRateLimiter::setReleaseTime(double newReleaseTime)
{
	if(newReleaseTime > 0.0001)      // 1/10 of a millisecond minimum
		releaseTime = newReleaseTime; 

	coeffRelease = exp( -1.0 / (sampleRate*releaseTime) );
}

//----------------------------------------------------------------------------
// others:

void SlewRateLimiter::reset()
{
 y_1 = 0;
}