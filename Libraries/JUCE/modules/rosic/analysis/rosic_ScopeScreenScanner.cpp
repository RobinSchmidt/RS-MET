#include "rosic_ScopeScreenScanner.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ScopeScreenScanner::ScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; 
  numCyclesShown = 5;
  sync = false;
  updateSawIncrement();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void ScopeScreenScanner::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  pitchDetector.setSampleRate(sampleRate);
  updateSawIncrement();
}

void ScopeScreenScanner::setMinFrequency(double newMinFreq)
{
  pitchDetector.setMinFundamental(newMinFreq);
}

void ScopeScreenScanner::setMaxFrequency(double newMaxFreq)
{
  pitchDetector.setMaxFundamental(newMaxFreq);
}

void ScopeScreenScanner::setSync(bool shouldSync)
{
  sync = shouldSync;
  updateSawIncrement();
}

//-------------------------------------------------------------------------------------------------
// others:

void ScopeScreenScanner::updateSawIncrement(double in)
{
  if(sync)
    pitchDetector.estimateFundamentalFrequency(in) / (sampleRate * numCyclesShown);
  else
    sawInc = scanFreq / (sampleRate * numCyclesShown);
}

void ScopeScreenScanner::reset()
{
  sawPhase = 0.0;
}