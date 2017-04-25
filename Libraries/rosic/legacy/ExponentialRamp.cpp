#include "ExponentialRamp.h"

//construction/destruction:
ExponentialRamp::ExponentialRamp()
{
	//set up parameters to neutral values (output in always 1.0):
	sampleRate = 44100.0;
	start      = 1.0;
	time       = 0.0;
	end        = 1.0;
	tauScale   = 1.0;

	// there was no previous output:
	out        = 0.0;

	//filter coefficients are calculated by the member functions:
	setStart(start);
	setTime (time);
	setEnd  (end);
}

ExponentialRamp::~ExponentialRamp()
{

}

//parameter settings:
void ExponentialRamp::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

	//re-calculate filter coefficients:
	setTime(time);
}

void ExponentialRamp::setStart(double newStart)
{
	start = newStart;
}

void ExponentialRamp::setTime(double newTime)
{
	double tau;

	if( newTime > 0.0 )
	{
		time  = newTime;                  //time parameter is expected in ms
		//tau   = 0.001*sampleRate*time;    //calculate the time constant
  tau   = sampleRate*time;    //calculate the time constant
  tau  *= tauScale;                 //scale the time constant

		//calculate filter coefficients:
  recursionCoeff = exp( -1.0 / tau ); 
  inputCoeff     = 1.0 - recursionCoeff;
	}
	else
	{
		time = 0.0;
  recursionCoeff = 0.0; //this pair of coefficients leads to no filtering at all
  inputCoeff     = 1.0; //the signal is passed through as it is
	}	
}

void ExponentialRamp::setEnd(double newEnd)
{
	end = newEnd;
}

void ExponentialRamp::setTauScale(double newTauScale)
{
	if( newTauScale > 0.0 )
		tauScale = newTauScale;

	// re-calculate filter coefficients:
	setTime(time); //setTime exspects ms
}

//event processing:
void ExponentialRamp::trigger()
{
 out = start;
}
