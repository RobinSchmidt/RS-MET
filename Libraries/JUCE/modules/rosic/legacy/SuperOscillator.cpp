#include "SuperOscillator.h"

SuperOscillator::SuperOscillator()
{
	numVoices     = 1;
	detune        = 0.0;
 freqSpacing   = 0;
 freqsAreDirty = true;
 phaseSpread   = 2.0/3.0;
	ampScale      = 1.0 / sqrt((double)numVoices);

	for(int i=0; i<maxVoices; i++)
	{
		increments1[i]  = 0.0;
		increments2[i]  = 0.0;
		phaseIndices[i] = 0.0;
	}

 detuneRatio = 0.5*(sqrt(5.0)-1.0);
 detuneDamp  = 1.0;

 sampleCount = 0;
}

SuperOscillator::~SuperOscillator()
{

}

//----------------------------------------------------------------------------
//parameter settings:

void SuperOscillator::setSampleRate(double newSampleRate)
{
	if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

	sampleRateRec   = 1.0 / sampleRate;

 basicIncrement  = tableLengthDbl*freq*sampleRateRec;

	increments1[0] = basicIncrement*pulseFactor1;
	increments2[0] = basicIncrement*pulseFactor2;

 freqsAreDirty = true;
}

void SuperOscillator::setStartPhase(double newStartPhase)
{
 if( (newStartPhase>=0.0) && (newStartPhase<=360.0) )
  startIndex = (newStartPhase/360.0)*tableLengthDbl;
}

void SuperOscillator::setNumVoices(int newNumVoices)
{
	if( isOdd(newNumVoices) ) //check, if it is an odd number - only those are allowed
		numVoices = newNumVoices;
	ampScale   = 1.0 / sqrt((double)numVoices);
 freqsAreDirty = true;
}

void SuperOscillator::setFreqSpacing(int newFreqSpacing)
{
	if( newFreqSpacing >= 0 && newFreqSpacing <= 127)
  freqSpacing = newFreqSpacing;
 freqsAreDirty = true;
}

void SuperOscillator::setDetuneRatio(double newDetuneRatio)
{
 if( (newDetuneRatio>=0.0) && (newDetuneRatio<=1.0) )
  detuneRatio = newDetuneRatio;
 freqsAreDirty = true;
}

void SuperOscillator::setDetunePhaseSpread(double newDetunePhaseSpread)
{
 phaseSpread = newDetunePhaseSpread;
}

//----------------------------------------------------------------------------
//event processing:
void SuperOscillator::resetPhase()
{
	for(long i=0; i<numVoices; i++)
  phaseIndices[i] = startIndex + phaseSpread*(i*tableLengthDbl/numVoices);
 sampleCount = 0;
}


//-----------------------------------------------------------------------------
//others:
