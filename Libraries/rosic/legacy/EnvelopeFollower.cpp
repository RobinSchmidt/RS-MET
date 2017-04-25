// EnvelopeFollower.cpp: implementation of the EnvelopeFollower class.
//
//////////////////////////////////////////////////////////////////////

#include "EnvelopeFollower.h"


//construction/destruction:
EnvelopeFollower::EnvelopeFollower()
{
	//parameter initializition:
 sampleRate  = 44100.0; //default sample rate
	attackTime  = 0.01;    //10 ms attack time constant by default
	releaseTime = 0.05;    //50 ms release time constant by default

	//calculation of internal variables:
	coeffAttack  = exp( -1 / (sampleRate*attackTime)  );
	coeffRelease = exp( -1 / (sampleRate*releaseTime) );

	//init previous output sample:
	y_1 = 0;

	//init mode to averaged absolute value:
	mode = 1;
}

EnvelopeFollower::~EnvelopeFollower()
{

}

//-----------------------------------------------------------------------------------------------------------
//parameter settings:
void EnvelopeFollower::setSampleRate(sample SampleRate)
{
	if(SampleRate>0)
		sampleRate = SampleRate;

 if(attackTime>0)
	 coeffAttack  = exp( -1 / (sampleRate*attackTime)  );
 else
  coeffAttack = 0.0;
 if(releaseTime>0)
	 coeffRelease  = exp( -1 / (sampleRate*releaseTime)  );
 else
  coeffRelease = 0.0;
}

void EnvelopeFollower::setMode(short Mode)
{
 if( (Mode>=1) && (Mode<=3) )
		mode = Mode;
}

void EnvelopeFollower::setAttackTime(sample AttackTime)
{
 /*
	if(0.001*AttackTime>0)
 {
		attackTime   = 0.001*AttackTime; //convert milliseconds to seconds
	 coeffAttack  = exp( -1 / (sampleRate*attackTime)  );
 }
	else
 {
		attackTime  = 0.0;
		coeffAttack = 0.0;
 }
 */
 if(AttackTime>0)
  attackTime = 0.001*AttackTime; //convert milliseconds to seconds
 else
  attackTime = 1/sampleRate;     //only one sample attack-time (almost instantaneous attack)

 coeffAttack  = exp( -1 / (sampleRate*attackTime)  );
 coeffRelease = exp( -1 / (sampleRate*releaseTime) );
}

void EnvelopeFollower::setReleaseTime(sample ReleaseTime)
{
 /*
	if(0.001*ReleaseTime>0)
 {
		releaseTime  = 0.001*ReleaseTime; //convert milliseconds to seconds
	 coeffRelease = exp( -1 / (sampleRate*releaseTime) );
 }
	else
 {
		releaseTime  = 0.0;
	 coeffRelease = 0.0;
 }
 */
 if(ReleaseTime>0)
  releaseTime = 0.001*ReleaseTime; //convert milliseconds to seconds
 else
  releaseTime = 1/sampleRate;      //only one sample release-time (almost instantaneous release)

 coeffAttack  = exp( -1 / (sampleRate*attackTime)  );
 coeffRelease = exp( -1 / (sampleRate*releaseTime) );
}

//-----------------------------------------------------------------------------------------------------------
//audio processing:
sample EnvelopeFollower::getSample(sample InputSample)
{
	static sample temp, coeff;

	if(mode==1)
		//take absolute value in mode 1:
		temp = absDbl(InputSample);
	else
 	//square the signal in mode 2 and 3:
	 temp = InputSample*InputSample; 

	//decide, if signal is rising or falling and choose appropriate filter coefficient:
	if( y_1 < temp )
		coeff = coeffAttack;
	else
		coeff = coeffRelease;

	//perform the filter operation:
	temp = (1-coeff)*temp + coeff*y_1;

	//update the memorized output for next call:
	y_1 = temp;

	//extract the square-root in mode 3:
	if(mode==3)
		temp = sqrt(temp);

	return temp;
}

//-----------------------------------------------------------------------------------------------------------
//others:
void EnvelopeFollower::reset()
{
 y_1 = 0;
}