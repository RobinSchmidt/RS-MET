// ExponentialEnvelope.cpp: implementation of the ExponentialEnvelope class.
//
//////////////////////////////////////////////////////////////////////

#include "ExponentialEnvelope.h"

//-------------------------------------------------------------------------------------------------------
//construction/destruction:
ExponentialEnvelope::ExponentialEnvelope()
{
	sampleRate   = 44100.0;
	startValue   = 0.0;
	attackTime   = 0.0;
	peakValue    = 1.0;
 holdTime     = 0.0;
	decayTime    = 0.1;
	sustainValue = 0.5;
	releaseTime  = 0.01;
	endValue     = 0.0;

	time         = 0.0;
	timeScale    = 1.0;
 increment    = 1.0/sampleRate;
	tauScale     = 1.0;
	noteIsOn     = false;
 endIsReached = false;

	prevOutput  = 0.0;

 //test:
 setAttack(0.1);
 setRelease(0.1);
}

ExponentialEnvelope::~ExponentialEnvelope()
{

}

//-------------------------------------------------------------------------------------------------------
//parameter settings:
void ExponentialEnvelope::setSampleRate(sample SampleRate)
{
	if(SampleRate>0)
		sampleRate = sampleRate;

 //adjust time increment:
 increment = 1.0/SampleRate;

	//re-calculate coefficients for the 3 filters:
	setAttack (attackTime);
	setDecay  (decayTime);
	setRelease(releaseTime);
}

void ExponentialEnvelope::setStart(sample Start)
{
 startValue = Start;
}

void ExponentialEnvelope::setAttack(sample Attack)
{
	sample tauAttack;

	if(timeScale*tauScale*sampleRate*Attack>0)
	{
		attackTime        = Attack;
  scaledAttackTime  = timeScale*attackTime;
  scaledOverallTime = scaledAttackTime+scaledHoldTime+scaledDecayTime+scaledReleaseTime;

		//calculate the filters time constant tau (which is the time of the impulse response to decay to
		//36.8 % - in the case of a step-function entering the filter (as we have here) it is also the
		//the time, which the filter needs to reach 63.2 % of the step height:
		tauAttack = (sampleRate*scaledAttackTime);

		//we don't want the envelope to reach only 63.2 % of the peak value - we want it to be quite close
		//to the peak value - so we make the time constant smaller in order to reach the peak faster:
		tauAttack = tauScale*tauAttack;

		//calculate the feedforward and feedback coefficient:
		recursionCoeffAttack  = exp( -1 / tauAttack ); 
		inputCoeffAttack      = 1 - recursionCoeffAttack;
	}
	else
	{
		attackTime            = 0.0;
		scaledAttackTime      = 0.0;
		recursionCoeffAttack  = 0.0; //this pair of coefficients leads to no filtering at all
		inputCoeffAttack      = 1.0; //the signal is passed through as it is
	}

}

void ExponentialEnvelope::setPeak(sample Peak)
{
 peakValue = Peak;
}

void ExponentialEnvelope::setHold(sample Hold)
{
 if(Hold>=0) 
  holdTime = Hold;
 scaledHoldTime    = timeScale*holdTime;
 scaledOverallTime = scaledAttackTime+scaledHoldTime+scaledDecayTime+scaledReleaseTime;

}

void ExponentialEnvelope::setDecay(sample Decay)
{
	sample tauDecay;

	if(timeScale*tauScale*sampleRate*Decay>0)
	{
		decayTime         = Decay;
  scaledDecayTime   = timeScale*decayTime;
  scaledOverallTime = scaledAttackTime+scaledHoldTime+scaledDecayTime+scaledReleaseTime;

		//calculate the filters time constant tau (see calculation of the attack constant for information
		//about the scale factors):
		tauDecay = (sampleRate*scaledDecayTime);
		tauDecay = tauScale*tauDecay;

		//calculate the feedforward and feedback coefficient:
		recursionCoeffDecay  = exp( -1 / tauDecay  );
		inputCoeffDecay      = 1 - recursionCoeffDecay;
	}
	else
	{
		decayTime            = 0.0;
		scaledDecayTime      = 0.0;
		recursionCoeffDecay  = 0.0; //this pair of coefficients leads to no filtering at all
		inputCoeffDecay      = 1.0; //the signal is passed through as it is
	}

}

void ExponentialEnvelope::setSustain(sample Sustain)
{
 sustainValue = Sustain;
}

void ExponentialEnvelope::setRelease(sample Release)
{
	sample tauRelease;

	if(timeScale*tauScale*sampleRate*Release>0)
	{
		releaseTime       = Release;
  scaledReleaseTime = timeScale*releaseTime;
  scaledOverallTime = scaledAttackTime+scaledHoldTime+scaledDecayTime+scaledReleaseTime;

		//calculate the filters time constant tau (see calculation of the attack constant for information
		//about the scale factors):
		tauRelease = (sampleRate*scaledReleaseTime);
		tauRelease = tauScale*tauRelease;

  //calculate the feedforward and feedback coefficient:
		recursionCoeffRelease  = exp( -1 / tauRelease  );
		inputCoeffRelease      = 1 - recursionCoeffRelease;
	}
	else
	{
		releaseTime            = 0.0;
		scaledReleaseTime      = 0.0;
		recursionCoeffRelease  = 0.0; //this pair of coefficients leads to no filtering at all
		inputCoeffRelease      = 1.0; //the signal is passed through as it is
	}

}

void ExponentialEnvelope::setEnd(sample End)
{
 endValue = End;
}

void ExponentialEnvelope::setTimeScale(sample TimeScale)
{
	if(TimeScale>0.001)
  timeScale = TimeScale;

	//re-calculate coefficients for the 3 filters:
	setAttack (attackTime);
 setHold   (holdTime);
	setDecay  (decayTime);
	setRelease(releaseTime);
}

void ExponentialEnvelope::setTauScale(sample TauScale)
{
 if(TauScale>0.001)
  tauScale = TauScale;
 setAttack(attackTime);
 setDecay(decayTime);
 setRelease(releaseTime);
}

//-------------------------------------------------------------------------------------------------------
//audio processing:
sample ExponentialEnvelope::getSample()
{
 static sample outSamp;

 //envelope is in attack or hold phase:
 if(time <= scaledAttackTime+scaledHoldTime)  //noteIsOn has not to be checked, because, time is advanced to the beginning
                          //of the release phase in noteOff()
 {
  outSamp = inputCoeffAttack*peakValue + recursionCoeffAttack*prevOutput;
  time   += increment;
 }

 //envelope is in decay phase:
 else if(time <= (scaledAttackTime+scaledHoldTime+scaledDecayTime)) //noteIsOn has not to be checked, see above
 {
  outSamp = inputCoeffDecay*sustainValue + recursionCoeffDecay*prevOutput;
  time   += increment;
 }

 //envelope is in sustain phase:
 //else if( (time <= (attackTime+decayTime+releaseTime)) && noteIsOn)
	else if(noteIsOn)
  outSamp = inputCoeffDecay*sustainValue + recursionCoeffDecay*prevOutput;
  //time is not incremented in sustain

 //envelope is in release phase:
 else if( (time <= (scaledAttackTime+scaledHoldTime+scaledDecayTime+scaledReleaseTime)) && !noteIsOn)
 {
  outSamp = inputCoeffRelease*endValue + recursionCoeffRelease*prevOutput;
  time   += increment;
 }
 else if(time <= scaledOverallTime+2*releaseTime)
 {
  outSamp = inputCoeffRelease*endValue + recursionCoeffRelease*prevOutput;
  time   += increment;
 }
 //envelope has reached its end point:
 else
 {
  outSamp = inputCoeffRelease*endValue + recursionCoeffRelease*prevOutput;
  time   += increment;
  endIsReached = true;
 }

	//store output sample for next call:
	prevOutput = outSamp;

 return outSamp;
}

//------------------------------------------------------------------------------------------------------------
//others:
void ExponentialEnvelope::reset()
{
 time = 0.0;
}

void ExponentialEnvelope::trigger()
{
	if( endIsReached )
  //start envelope from start point only when the end was reached, otherwise start whereever
		//it is:
  prevOutput = startValue;

 //reset time for the new note
 time         = 0.0;   
 noteIsOn     = true;
 endIsReached = false;
}

void ExponentialEnvelope::noteOff()
{
 noteIsOn = false; //begin the release phase

 //advance time to the beginnig of the release phase:
 time = (scaledAttackTime+scaledHoldTime+scaledDecayTime+increment);
}

/*
bool ExponentialEnvelope::endIsReached()
{
 if(time >= scaledOverallTime+2*releaseTime)
  return true;
 else
  return false;
}
*/