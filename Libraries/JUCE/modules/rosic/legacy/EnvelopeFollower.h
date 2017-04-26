#ifndef EnvelopeFollower_h
#define EnvelopeFollower_h

/**

This is an envelope follower with user adjustable attack and release time
constants. It squares the incoming signal, smoothes this signal via an
attack-release (AR) averager (see Undo Zoelzer: DAFX, page 84) and extracts
the square-root. So the output value represents the instantaneous RMS-value of
the incoming signal. The averaging time is determined by the attack and
release time constants.

© Braindoc 2004 (www.braindoc.de)

*/

/*
#if !defined(AFX_ENVELOPEFOLLOWER_H__E44408BB_3158_42A1_8015_02D8C22206F9__INCLUDED_)
#define AFX_ENVELOPEFOLLOWER_H__E44408BB_3158_42A1_8015_02D8C22206F9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
*/

//#include "math.h"
//#include "VSTools.h"

class EnvelopeFollower  
{
public:
	//construction/destruction:
	EnvelopeFollower();
	virtual ~EnvelopeFollower();

	//public member functions:
	//parameter settings:
	void setSampleRate (sample SampleRate);
	void setMode       (short  mode);       //3 modes: 1: AR-averaged absolute value, 
	                                        //         2: AR-averaged squared value,
	                                        //         3: square root of AR-averaged squared value (RMS)
	void setAttackTime (sample AttackTime); //exspects the value in milliseconds
	void setReleaseTime(sample ReleaseTime);

	//audio processing:
	sample getSample(sample InputSample);

	//others:
	void reset();

protected:
	//protected member variables:
	sample sampleRate;
	short  mode;                    //switch between averaged absolute, squared or RMS value
	sample attackTime, releaseTime; //time constants in seconds

	//the recursion coefficient for the recursive RC-Filter - the smoothing is performed via the recursive
	//filter equation: y[n] = (1-coeff)*x[n] + coeff*y[n-1]. Which coeff of the two will be used is 
	//determined by the signal - if it is rising ( y[n-1] < x^2[n] ), the attack-coefficient is used,
 //if it is falling, the release coeeficient will be used
	sample coeffAttack, coeffRelease; 

	//previous ouput sample of the unit:
 sample y_1;

};

#endif // !defined(AFX_ENVELOPEFOLLOWER_H__E44408BB_3158_42A1_8015_02D8C22206F9__INCLUDED_)
