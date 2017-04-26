#include "rosic_LowFrequencyOscillator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

LowFrequencyOscillator::LowFrequencyOscillator()
{
  positionInt   = 0;
  positionFrac  = 0.0;
  incrementInt  = 1;
  incrementFrac = 0.0;
  waveTable     = new WaveTable;
  parameters    = new LowFrequencyOscillatorParameters;
  tableIsOwned  = true;
  setSampleRate(44100.0);
  slewRateLimiterL.setAttackTime(0.0);
  slewRateLimiterL.setReleaseTime(0.0);
  slewRateLimiterR.setAttackTime(0.0);
  slewRateLimiterR.setReleaseTime(0.0);
}

LowFrequencyOscillator::~LowFrequencyOscillator()
{  
  delete parameters;
  if( tableIsOwned )
    delete waveTable;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void LowFrequencyOscillator::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
  {
    parameters->sampleRate = newSampleRate;
    slewRateLimiterL.setSampleRate(newSampleRate);  
    slewRateLimiterR.setSampleRate(newSampleRate);  
    calculateIncrement();
  }
}

void LowFrequencyOscillator::setWaveTableToUse(WaveTable *newTableToUse)
{
  if( newTableToUse != waveTable )
  {
    if( tableIsOwned )
      delete waveTable;
    waveTable    = newTableToUse;
    tableIsOwned = false;
    reset();
    calculateIncrement();
  }
}
  
/*
void LowFrequencyOscillator::setWaveform(int newWaveform)
{ 
  if(waveTable != NULL) 
    waveTable->setWaveform(newWaveform); 
}

void LowFrequencyOscillator::setWaveform(double *newWaveformBuffer, int newLength, char* name)
{
  if(waveTable != NULL) 
    waveTable->setWaveform(newWaveformBuffer, newLength, name); 
}

void LowFrequencyOscillator::setWaveform(float *newWaveformBuffer, int newLength, char* name)
{
  if(waveTable != NULL) 
    waveTable->setWaveform(newWaveformBuffer, newLength, name); 
}
*/

void LowFrequencyOscillator::setCycleLength(double newCycleLength)
{
  if( newCycleLength != parameters->cycleLength )
  {
    parameters->cycleLength = newCycleLength;
    calculateIncrement();
  }
}

void LowFrequencyOscillator::setBeatsPerMinute(double newBPM)
{
  if( newBPM > 0.0 && newBPM != parameters->bpm )
  {
    parameters->bpm = newBPM;
    calculateIncrement();
  }
}

void LowFrequencyOscillator::setTempoSync(bool shouldSync)
{
  parameters->sync = shouldSync;
  calculateIncrement();
}

void LowFrequencyOscillator::setStartPhase(double newStartPhase)
{
  if( newStartPhase >= 0.0 && newStartPhase <= 360 )
    parameters->startPhase = newStartPhase;
}

void LowFrequencyOscillator::setStereoPhase(double newStereoPhase)
{
  if( newStereoPhase >= 0.0 && newStereoPhase <= 360 )
    parameters->stereoPhase = newStereoPhase;
}

void LowFrequencyOscillator::setMute(bool shouldBeMuted)
{
  parameters->mute = shouldBeMuted;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double LowFrequencyOscillator::getFrequency() const 
{
  double cycleLengthInSeconds;
  if( parameters->sync )
    cycleLengthInSeconds = beatsToSeconds(parameters->cycleLength, parameters->bpm);
  else
    cycleLengthInSeconds = parameters->cycleLength;
  return 1.0/cycleLengthInSeconds;
}

/*
bool LowFrequencyOscillator::isMuted()
{
  return parameters->mute;
}

double LowFrequencyOscillator::getStartPhase()
{
  return parameters->startPhase;
}
*/

/*
void LowFrequencyOscillator::getWaveformForDisplay(double *targetBuffer, int numSamplesToShow)
{
  for(int n=0; n<numSamplesToShow; n++)
    targetBuffer[n] = 0.0;

  // preliminary - use min/max extraction later:
  if( waveTable != NULL )
    copyBufferWithLinearInterpolation(waveTable->prototypeBufferL, waveTable->tableLength, 
    targetBuffer, numSamplesToShow);
}

char* LowFrequencyOscillator::getWaveformName() const
{
  if( waveTable != NULL )
    return waveTable->getWaveformName();
  else
    return "Error";
}
*/

//-------------------------------------------------------------------------------------------------
// others:

void LowFrequencyOscillator::trigger()
{
  triggerWithPhase(parameters->startPhase);
}

void LowFrequencyOscillator::triggerWithPhase(double phase)
{
  double startPosition = waveTable->tableLength * (phase/360.0);
  positionInt          = floorInt(startPosition);
  positionFrac         = startPosition - positionInt;
}

void LowFrequencyOscillator::reset()
{
  slewRateLimiterL.reset();
  slewRateLimiterR.reset();
  trigger();
}