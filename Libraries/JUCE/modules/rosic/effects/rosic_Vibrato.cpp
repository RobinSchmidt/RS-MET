#include "rosic_Vibrato.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Vibrato::Vibrato(int bufferLengthToAllocate) : DelayLineStereo(bufferLengthToAllocate)
{
  dryWetRatio  = 1.0;
  setAverageDelayTime(10.0);
}

Vibrato::~Vibrato()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void Vibrato::setSampleRate(double newSampleRate)
{
  ModulationEffect::setSampleRate(newSampleRate);
  calculateDepthInSamples();
}

void Vibrato::setCycleLength(double newCycleLength)
{
  ModulationEffect::setCycleLength(newCycleLength);
  calculateDepthInSamples();
}

void Vibrato::setDepth(double newDepth)
{
  ModulationEffect::setDepth(newDepth);
  calculateDepthInSamples();
}

void Vibrato::setTempoSync(bool shouldTempoSync)
{
  ModulationEffect::setTempoSync(shouldTempoSync);
  calculateDepthInSamples();
}

void Vibrato::setTempoInBPM(double newTempoInBPM)
{
  ModulationEffect::setTempoInBPM(newTempoInBPM);
  calculateDepthInSamples();
}

void Vibrato::setAverageDelayTime(double newAverageDelayTime)
{
  if( newAverageDelayTime >= 0.0 )
    averageDelay = newAverageDelayTime;
  dA = roundToInt(0.001*averageDelay*sampleRate);           // average delay in samples
  calculateDepthInSamples();
}

//-------------------------------------------------------------------------------------------------
// others:

void Vibrato::reset()
{
  clearBuffers();
  resetOscillatorPhases();
}

void Vibrato::calculateDepthInSamples()
{
  double ratio = pitchOffsetToFreqFactor(0.5*depth);        // resampling ratio
  d  = (ratio-1.0)*sampleRate / (2*PI*lfo.getFrequency()); // delay modulation depth in samples
  if( d > dA-3.0 )
  {
    d = dA-3.0;
    // DEBUG_BREAK; // delay excursion larger than average delay (minus interpolator margin)
  }
}