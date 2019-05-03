//#include "rosic_SingleSidebandModulator.h"
//using namespace rosic;

//=================================================================================================
// class SingleSidebandModulator:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SingleSidebandModulator::SingleSidebandModulator()
{
  yOld               = 0.0;
  modulatorFrequency = 0.0;
  feedbackFactor     = 0.0;
  usbFactor          = 1.0;
  lsbFactor          = 1.0;
  antiAlias          = false;
  setSampleRate(44100.0);
  upsamplingFilter.setSubDivision(2);
  downsamplingFilter.setSubDivision(2);
    
  mutex.lock();
  sineOscillator.trigger();
  mutex.unlock();
}

SingleSidebandModulator::~SingleSidebandModulator()
{

}

//-------------------------------------------------------------------------------------------------
// setup:
    
void SingleSidebandModulator::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    mutex.lock();
    sampleRate = newSampleRate;
    if( antiAlias == true )
      sineOscillator.setSampleRate(2*sampleRate);
    else
      sineOscillator.setSampleRate(sampleRate);
    sineOscillator.trigger();
    mutex.unlock();
  }
}

void SingleSidebandModulator::setModulatorFrequency(double newFrequency)
{
  mutex.lock();
  sineOscillator.setFrequency(newFrequency); 
  sineOscillator.trigger();
  mutex.unlock();
}

void SingleSidebandModulator::setAntiAliasing(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
  mutex.lock();
  if( antiAlias == true )
    sineOscillator.setSampleRate(2*sampleRate);
  else
    sineOscillator.setSampleRate(sampleRate);
  sineOscillator.trigger();
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// others:
    
void SingleSidebandModulator::reset()
{
  yOld = 0.0;
  upsamplingFilter.reset();
  downsamplingFilter.reset();
  quadratureNetwork.reset();
  mutex.lock();
  sineOscillator.trigger();
  mutex.unlock();
}


//=================================================================================================
// class SingleSidebandModulatorStereo:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SingleSidebandModulatorStereo::SingleSidebandModulatorStereo()
{
  dry                 = 0.0;
  wet                 = 1.0;
  modulatorFrequency  = 0.0;
  stereoOffset        = 0.0;
}

SingleSidebandModulatorStereo::~SingleSidebandModulatorStereo()
{

}

//-------------------------------------------------------------------------------------------------
// setup:
    


