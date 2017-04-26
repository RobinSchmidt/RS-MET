#include "rosic_SineOscillator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SineOscillator::SineOscillator()
{
  sampleRate = 44100.0; 
  frequency  = 1000.0;
  startPhase = 0.0;
  omega      = 2.0*PI*frequency/sampleRate;
  a1         = 2.0*cos(omega);
  trigger();
}

SineOscillator::~SineOscillator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void SineOscillator::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.0)
    sampleRate = newSampleRate;
  setOmega(2.0*PI*frequency/sampleRate);
}

void SineOscillator::setFrequency(double newFrequency)
{
  if(newFrequency > 0.0)
    frequency = newFrequency;
  setOmega(2.0*PI*frequency/sampleRate);
}

void SineOscillator::setOmega(double newOmega)
{
  if( newOmega != omega )
  {
    mutex.lock();
    double currentPhase = asin(getSample()); 
    omega     = newOmega;
    frequency = omega*sampleRate / (2.0*PI);
    a1        = 2.0*cos(omega);
    triggerWithPhase(currentPhase);
    mutex.unlock();
  }
}

void SineOscillator::setStartPhase(double newStartPhase)
{
  startPhase = newStartPhase;
}

//-------------------------------------------------------------------------------------------------
// event handling:

void SineOscillator::trigger()
{   
  mutex.lock();
  s1 = sin(startPhase -     omega);
  s2 = sin(startPhase - 2.0*omega);
  mutex.unlock();
}

void SineOscillator::triggerWithPhase(double phase)
{  
  mutex.lock();
  s1 = sin(phase -     omega);
  s2 = sin(phase - 2.0*omega);
  mutex.unlock();
}


