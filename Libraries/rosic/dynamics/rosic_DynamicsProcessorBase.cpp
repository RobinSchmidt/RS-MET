#include "rosic_DynamicsProcessorBase.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DynamicsProcessorBase::DynamicsProcessorBase(int newLookAheadBufferSize)
{
  lookAheadBufferLength = rmax(newLookAheadBufferSize, 1);
  lookAheadBufferL      = new double[lookAheadBufferLength];
  lookAheadBufferR      = new double[lookAheadBufferLength];
  tapIn                 = 0;
  tapOut                = 0;
  lookAheadTime         = 0.0;
  sampleRate            = 44100.0;
  dry                   = 0.0;
  wet                   = 1.0;
  inputGainFactor       = 1.0;
  outputGainFactor      = 1.0;

  levelDetector.setAttackTime(1.0);
  levelDetector.setReleaseTime(1.0);
  reset();
  setLookAheadTime(0.0); 
}

DynamicsProcessorBase::~DynamicsProcessorBase()
{
  delete[] lookAheadBufferL;
  delete[] lookAheadBufferR;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void DynamicsProcessorBase::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupLookAhead();
    levelDetector.setSampleRate(sampleRate);
    attackReleaseEnveloper.setSampleRate(sampleRate);
  }
}

void DynamicsProcessorBase::setAttackTime(double newAttackTime)
{ 
  attackReleaseEnveloper.setAttackTime(newAttackTime); 
  levelDetector.setAttackTime(newAttackTime); 
}
  
void DynamicsProcessorBase::setReleaseTime(double newReleaseTime)
{ 
  attackReleaseEnveloper.setReleaseTime(newReleaseTime); 
  levelDetector.setReleaseTime(newReleaseTime); 
}

void DynamicsProcessorBase::setLookAheadTime(double newLookAheadTime)
{
  lookAheadTime = newLookAheadTime;  
  setupLookAhead();
}

//-------------------------------------------------------------------------------------------------
// others:

void DynamicsProcessorBase::reset()
{
  for(int i=0; i<lookAheadBufferLength; i++)
  {
    lookAheadBufferL[i] = 0.0;
    lookAheadBufferR[i] = 0.0;
  }
  attackReleaseEnveloper.reset();
  levelDetector.reset();
}

void DynamicsProcessorBase::setupLookAhead()
{
  double lookAheadInSeconds = 0.001*lookAheadTime;
  double lookAheadInSamples = sampleRate*lookAheadInSeconds;
  lookAheadInSamples        = clip(lookAheadInSamples, 0.0, (double)lookAheadBufferLength);
  int dInt                  = (int) round(lookAheadInSamples);
  tapOut                    = tapIn - dInt;
  tapOut                    = wrapAround(tapOut, lookAheadBufferLength);
}
