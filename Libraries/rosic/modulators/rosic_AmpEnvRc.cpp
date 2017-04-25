#include "rosic_AmpEnvRc.h"
using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:
AmpEnvRc::AmpEnvRc()
{
	sampleRate     = 44100.0;
	sampleRateRec  = 1.0/sampleRate;
	startAmp       = 0.0;
	attackTime     = 0.0;
	peakAmp        = 1.0;
 holdTime       = 0.0;
	decayTime      = 0.1;
	sustainAmp     = 0.5;
	releaseTime    = 0.01;
	endAmp         = 0.0+DB2AMP(-300.0); // is not exactly zero to 
                                      // avoid denorm-problems
	time           = 0.0;
	timeScale      = 1.0;
	increment      = timeScale/sampleRate;
	tauScale       = 1.0;
	peakScale      = 1.0;
	noteIsOn       = false;
 outputIsZero   = true;

	previousOutput = 0.0;

 // call some internal functions to trigge the coefficient calculations:
 setAttack(attackTime);
 setDecay(decayTime);
 setRelease(releaseTime);
}

AmpEnvRc::~AmpEnvRc()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void AmpEnvRc::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

 // adjust time increment:
 increment = timeScale/sampleRate;

	//re-calculate coefficients for the 3 filters:
	setAttack (attackTime);
	setDecay  (decayTime);
	setRelease(releaseTime);
}

void AmpEnvRc::setStart(double newStartValue)
{
 startAmp = DB2AMP(newStartValue);
}

void AmpEnvRc::setAttack(double newAttackTime)
{
	double lTauAttack;

	if( newAttackTime > 0 )
	{
		attackTime = newAttackTime;

		// calculate the filters time constant tau (which is the time of the impulse
  // response to decay to 36.8 % - in the case of a step-function entering the
  // filter (as we have here) it is also the the time, which the filter needs
  // to reach 63.2 % of the step height:
		lTauAttack = (sampleRate*attackTime);

		// we don't want the envelope to reach only 63.2 % of the peak value - we
  // want it to be quite close to the peak value - so we may want make the
  // time constant smaller in order to reach the peak faster:
		lTauAttack *= tauScale;

		// we have also take into account, that there is a timeScale variable in this
  // class, which should affect the time constant too:
		lTauAttack /= timeScale;

		// calculate the feedforward and feedback coefficient:
		recursionCoeffAttack  = exp( -1 / lTauAttack ); 
		inputCoeffAttack      = 1 - recursionCoeffAttack;
	}
	else // newAttackTime <= 0
	{
		attackTime            = 0.0;
		recursionCoeffAttack  = 0.0; 
		inputCoeffAttack      = 1.0; 
	}

 // calculate the accumulated time-values:
 attPlusHld               = attackTime + holdTime;
 attPlusHldPlusDec        = attPlusHld + decayTime;
 attPlusHldPlusDecPlusRel = attPlusHldPlusDec + releaseTime;
}

void AmpEnvRc::setPeak(double newPeakValue)
{
 peakAmp = DB2AMP(newPeakValue);
}

void AmpEnvRc::setHold(double newHoldTime)
{
 if( newHoldTime >= 0 ) 
  holdTime = newHoldTime;

 // calculate the accumulated time-values:
 attPlusHld               = attackTime + holdTime;
 attPlusHldPlusDec        = attPlusHld + decayTime;
 attPlusHldPlusDecPlusRel = attPlusHldPlusDec + releaseTime;
}

void AmpEnvRc::setDecay(double newDecayTime)
{
	double lTauDecay;

	if( newDecayTime > 0 )
	{
		decayTime = newDecayTime;

		// calculate the filters time constant tau (see calculation of the attack
  // constant for information about the scale factors):
		lTauDecay  = (sampleRate*decayTime);
		lTauDecay *= tauScale;
		lTauDecay /= timeScale;

		//calculate the feedforward and feedback coefficient:
		recursionCoeffDecay  = exp( -1 / lTauDecay  );
		inputCoeffDecay      = 1 - recursionCoeffDecay;
	}
	else // newDecayTime <= 0
	{
		decayTime            = 0.0;
		recursionCoeffDecay  = 0.0;
		inputCoeffDecay      = 1.0; 
	}

 // calculate the accumulated time-values:
 attPlusHld               = attackTime + holdTime;
 attPlusHldPlusDec        = attPlusHld + decayTime;
 attPlusHldPlusDecPlusRel = attPlusHldPlusDec + releaseTime;
}

void AmpEnvRc::setSustain(double newSustainValue)
{
 sustainAmp = DB2AMP(newSustainValue);
}

void AmpEnvRc::setRelease(double newReleaseTime)
{
	double lTauRelease;

	if( newReleaseTime > 0 )
	{
		releaseTime = newReleaseTime;

		// calculate the filters time constant tau (see calculation of the attack
  // constant for information about the scale factors):
		lTauRelease  = (sampleRate*releaseTime);
		lTauRelease *= tauScale;
		lTauRelease /= timeScale;

  // calculate the feedforward and feedback coefficient:
		recursionCoeffRelease  = exp( -1 / lTauRelease  );
		inputCoeffRelease      = 1 - recursionCoeffRelease;
	}
	else // newReleaseTime <= 0
	{
		releaseTime            = 0.0;
		recursionCoeffRelease  = 0.0;
		inputCoeffRelease      = 1.0; 
	}

 //calculate the accumulated time-values:
 attPlusHld               = attackTime + holdTime;
 attPlusHldPlusDec        = attPlusHld + decayTime;
 attPlusHldPlusDecPlusRel = attPlusHldPlusDec + releaseTime;
}

void AmpEnvRc::setEnd(double newEndValue)
{
 endAmp = DB2AMP(newEndValue);
}

void AmpEnvRc::setTimeScale(double newTimeScale)
{
	if( newTimeScale > 0 )
  timeScale = newTimeScale;

	increment  = timeScale/sampleRate;

	//re-calculate coefficients for the 3 filters:
	setAttack (attackTime);
	setDecay  (decayTime);
	setRelease(releaseTime);
}

void AmpEnvRc::setTauScale(double newTauScale)
{
 if( newTauScale > 0 )
  tauScale = newTauScale;

 setAttack(attackTime);
 setDecay(decayTime);
 setRelease(releaseTime);
}

void AmpEnvRc::setPeakScale(double newPeakScale)
{
 if( newPeakScale > 0 )
  peakScale = newPeakScale;
}
//----------------------------------------------------------------------------
// others:

void AmpEnvRc::reset()
{
 time = 0.0;
}

void AmpEnvRc::trigger(bool startFromCurrentValue)
{
	if( !startFromCurrentValue )
  previousOutput = startAmp;  // may lead to clicks

 //reset time for the new note
 time         = 0.0;   
 noteIsOn     = true;
 outputIsZero = false;
}

void AmpEnvRc::noteOff()
{
 noteIsOn = false; //begin the release phase

 //advance time to the beginnig of the release phase:
 time = (attackTime + holdTime + decayTime + increment);
}

bool AmpEnvRc::endIsReached()
{
 /*
 if( time > attPlusHldPlusDec + 7.0 * releaseTime )
  return true;
 else 
  return false;
 */

 if( noteIsOn == false && previousOutput < 0.001 )
  return true;
 else
  return false;
}