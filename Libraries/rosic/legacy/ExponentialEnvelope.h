// ExponentialEnvelope.h: interface for the ExponentialEnvelope class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(ExponentialEnvelope_h_Included)
#define ExponentialEnvelope_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "VSTools.h"

class ExponentialEnvelope  
{
public:
	//construction/destruction:
	ExponentialEnvelope();
	virtual ~ExponentialEnvelope();

	//parameter settings:
	void setSampleRate(sample SampleRate);  
 void setStart     (sample Start);     //this is the point where the envelope starts (in dB)
	void setAttack    (sample Attack);    //length of attack phase (in seconds)
 void setPeak      (sample Peak);      //highest point of the envelope (in dB)  
 void setHold      (sample Hold);      //hold time (in seconds)
	void setDecay     (sample Decay);     //length of decay phase (in seconds)
	void setSustain   (sample Sustain);   //sustain level
	void setRelease   (sample Release);   //length of release phase (in seconds)
	void setEnd       (sample End);       //end point of the envelope
	void setTimeScale (sample TimeScale); //scales the A,D and R times by adjusting the
  //increment (1 if not used), a timescale of 2 means the envelope is twice as fast, 0.5 means half 
  //as fast -> useful for implementing a key/velocity-tracking feature for the overall length of the envelope
 void setTauScale  (sample TauScale);

	//triggering:
 void trigger();  //causes the envelope to start with its attack-phase
	void noteOff();  //causes the envelope to start with its release-phase
	void reset();    //reset the time variable

	//info for the user:
	//bool endIsReached(); //tells the user if envelope has reached its end
 bool endIsReached;

	//audio processing:
 sample getSample();           

protected:
 //protected member variables:
	sample sampleRate, sampleRateRec; //sampleRateRec is the reciprocal of sampleRate

	//amplitude values for the stairstep input signal (as amplitude values - not in dB):
 sample startValue, peakValue, sustainValue, endValue;   
	
	//time values in seconds:
	sample attackTime, holdTime, decayTime, releaseTime;
 sample scaledAttackTime, scaledHoldTime,   //these times are scaled by timeScale
        scaledDecayTime, scaledReleaseTime;
 sample scaledOverallTime;  //sum of the 4 (scaled) times

	//time reference variables:
	sample time;      //this variable represents the time (the unit is seconds)
	sample timeScale; //scale the time constants
	sample increment; //increment for the time variable per sample (usually 1/sampleRate but can be
                   //modified for changing the overall duration of the envelope), also influences
	                  //the filter coefficients

	sample tauScale;  //scale factor for the time constants of the filters

	bool   noteIsOn;  //indicates if note is on, if not, the envelope starts with its release phase

	//filter coefficients:
	sample inputCoeffAttack,  recursionCoeffAttack, //the inputCoeffs scales the current input sample
		      inputCoeffDecay,   recursionCoeffDecay,  //whereas the recursionCoeff scales the previous 
		      inputCoeffRelease, recursionCoeffRelease;//output sample - which pair of coefficients will be
	       //used is determined by the state of the envelope (which is given by the time variable and
	       //note on/off status)

	//buffering:
	sample prevOutput;//holds the previous output sample

};

#endif // !defined(ExponentialEnvelope_h_Included)
