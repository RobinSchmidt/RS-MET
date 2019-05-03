//#include "rosic_RingModulator.h"
//using namespace rosic;

//=================================================================================================
// class RingModulator:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RingModulator::RingModulator()
{
  yOld               = 0.0;
  modulatorFrequency = 0.0;
  feedbackFactor     = 0.0;
  antiAlias          = false;
  setSampleRate(44100.0);
  upsamplingFilter.setSubDivision(2);
  downsamplingFilter.setSubDivision(2);
  sineOscillator.setStartPhase(0.0);
}

RingModulator::~RingModulator()
{

}

//-------------------------------------------------------------------------------------------------
// setup:
    
void RingModulator::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    mutex.lock();
    sampleRate = newSampleRate;
    if( antiAlias == true )
      sineOscillator.setSampleRate(2*sampleRate);
    else
      sineOscillator.setSampleRate(sampleRate);
    mutex.unlock();
  }
}

void RingModulator::setModulatorFrequency(double newFrequency)
{
  mutex.lock();
  sineOscillator.setFrequency(newFrequency); 
  mutex.unlock();
}

void RingModulator::setAntiAliasing(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
  mutex.lock();
  if( antiAlias == true )
    sineOscillator.setSampleRate(2*sampleRate);
  else
    sineOscillator.setSampleRate(sampleRate);
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// others:
    
void RingModulator::reset()
{
  yOld = 0.0;
  upsamplingFilter.reset();
  downsamplingFilter.reset();
  // nyquistBlocker.reset();
}


//=================================================================================================
// class RingModulatorStereo:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RingModulatorStereo::RingModulatorStereo()
{
  dry                 = 0.0;
  wet                 = 1.0;
  modulatorFrequency  = 0.0;
  stereoOffset        = 0.0;
}

RingModulatorStereo::~RingModulatorStereo()
{

}




