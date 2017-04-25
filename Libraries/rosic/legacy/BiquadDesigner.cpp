#include "BiquadDesigner.h"

//----------------------------------------------------------------------------
//construction/destruction:

BiquadDesigner::BiquadDesigner()
{
	setSampleRate(44100.0); // samplerate = 44100 Hz by default
	setMode      (0);       // bypass by default
	setFreq      (1000.0);  // cutoff/center frequency = 1 kHz by default 
	setQ         (1.0);     // Q = 1 by default 
	setGain      (0.0);     // Gain = 0 dB by default 
}

BiquadDesigner::~BiquadDesigner()
{}

//----------------------------------------------------------------------------
//parameter settings:

void BiquadDesigner::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

 sampleRateRec = 1.0 / sampleRate;
}

void BiquadDesigner::setMode(int newMode)
{
 mode = newMode;
}




