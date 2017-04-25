#include "AudioModule.h"

AudioModule::AudioModule()
{
 sampleRate = 44100.0;
}

AudioModule::~AudioModule()
{
}

void AudioModule::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;
}
