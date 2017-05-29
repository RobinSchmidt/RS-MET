//#include "rosic_AudioToMidi.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AudioToMidi::AudioToMidi()
{
  // initialize parameters:
  sampleRate     = 44100.0;
  minFreq        = 20.0;
  maxFreq        = 4000.0;

}

AudioToMidi::~AudioToMidi()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AudioToMidi::setSampleRate(double newSampleRate)
{
  pitchDetector.setSampleRate(sampleRate);

}

void AudioToMidi::setMinFrequency(double newMinFrequency)
{
  if( newMinFrequency > 10.0 && newMinFrequency < 5000.0 && newMinFrequency < maxFreq)
    minFreq = newMinFrequency;
}

void AudioToMidi::setMaxFrequency(double newMaxFrequency)
{
  if( newMaxFrequency > 100.0 && newMaxFrequency < 20000.0 && newMaxFrequency > minFreq)
    maxFreq = newMaxFrequency;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double AudioToMidi::getMinFrequency()
{
  return minFreq;
}

double AudioToMidi::getMaxFrequency()
{
  return maxFreq;
}



//-------------------------------------------------------------------------------------------------
// others:

void AudioToMidi::reset()
{

}

//-------------------------------------------------------------------------------------------------
// internal functions:



