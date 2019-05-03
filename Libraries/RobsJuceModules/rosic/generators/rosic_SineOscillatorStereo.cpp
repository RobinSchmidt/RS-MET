//#include "rosic_SineOscillatorStereo.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SineOscillatorStereo::SineOscillatorStereo()
{
  sampleRate = 44100.0; 
  frequency  = 1000.0;
  amplitude  = 1.0;
  startPhase = 0.0;
  omega      = 2.0*PI*frequency/sampleRate;
  a1         = 2.0*cos(omega);
  s1         = 0.0;
  s2         = 0.0;
  c1         = 0.0;
  c2         = 0.0;
  trigger();
}

SineOscillatorStereo::~SineOscillatorStereo()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void SineOscillatorStereo::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.0)
    sampleRate = newSampleRate;

  omega = 2.0*PI*frequency/sampleRate;
  a1    = 2.0*cos(omega);
}

void SineOscillatorStereo::setFrequency(double newFrequency)
{
  if(newFrequency > 0.0)
    frequency = newFrequency;

  omega = 2.0*PI*frequency/sampleRate;
  a1    = 2.0*cos(omega);
}

void SineOscillatorStereo::setAmplitude(double newAmplitude)
{
  amplitude = newAmplitude;
}

//-------------------------------------------------------------------------------------------------
// event handling:

void SineOscillatorStereo::trigger()
{
  s1 = sin(startPhase -     omega);
  s2 = sin(startPhase - 2.0*omega);
  c1 = cos(startPhase -     omega);
  c2 = cos(startPhase - 2.0*omega);
}


