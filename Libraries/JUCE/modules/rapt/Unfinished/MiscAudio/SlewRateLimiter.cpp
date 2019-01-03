// Construction/Destruction:

rsSlewRateLimiter::rsSlewRateLimiter()
{
	sampleRate = 44100.0f;
	attackTime = 10.0f;
	releaseTime = 100.0f;
	y1 = 0.0;

	calculateAttackCoefficient();
	calculateReleaseCoefficient();
}

// Setup:

void rsSlewRateLimiter::setSampleRate(double newSampleRate)
{
	if (newSampleRate > 0.0)
	{
		sampleRate = newSampleRate;
		calculateAttackCoefficient();
		calculateReleaseCoefficient();
	}
}

void rsSlewRateLimiter::setAttackTime(double newAttackTime)
{
	if (newAttackTime >= 0.0 && newAttackTime != attackTime)
	{
		attackTime = newAttackTime;
		calculateAttackCoefficient();
	}
}

void rsSlewRateLimiter::setReleaseTime(double newReleaseTime)
{
	if (newReleaseTime >= 0.0 && newReleaseTime != releaseTime)
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
	if (tau > 0.0)
		coeffAttack = exp(-1.0 / (sampleRate*tau));
	else
		coeffAttack = 0.0;
}

void rsSlewRateLimiter::calculateReleaseCoefficient()
{
	double tau = 0.001*releaseTime;
	if (tau > 0.0)
		coeffRelease = exp(-1.0 / (sampleRate*tau));
	else
		coeffRelease = 0.0;
}

double rsSlewRateLimiter::getSample(double in)
{
	if (y1 < in)
		y1 = in + coeffAttack * (y1 - in);
	else
		y1 = in + coeffRelease * (y1 - in);
	return y1;
}

/** Returns the attack-time (in milliseconds) - this time which it takes to rise to
 1-1/e = 0.63 when the input signal makes an upward step from 0 to 1. */
double rsSlewRateLimiter::getAttackTime() const { return attackTime; }

/** Returns the release-time (in milliseconds) - this time which it takes to fall to 1/e = 0.37
 when the input signal makes a downward step from 1 to 0. */
double rsSlewRateLimiter::getReleaseTime() const { return releaseTime; }


double rsSlewRateLimiterWithHold::getSample(double in)
{
	if (y1 > in) {
		if (sampleCounter >= holdSamples)
			y1 = in + coeffRelease * (y1 - in);
		else
			sampleCounter++;
	}
	else {
		y1 = in + coeffAttack * (y1 - in);
		sampleCounter = 0;
	}
	return y1;
}

void rsSlewRateLimiterWithHold::setSampleRate(double newSampleRate) //override ..is only for runtime polymorphism
{
	rsSlewRateLimiter::setSampleRate(newSampleRate);
	holdSamples = (int)round(holdTime * sampleRate);
}

void rsSlewRateLimiterWithHold::setHoldTime(double newHoldTime)
{
	holdTime = newHoldTime;
	holdSamples = (int)round(holdTime * sampleRate);
}

void rsSlewRateLimiterWithHold::reset() // override
{
	y1 = 0.0;
	sampleCounter = 0;
}