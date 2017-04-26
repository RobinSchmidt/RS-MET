#include "rosic_WahWah.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WahWah::WahWah()
{
  pitch  = freqToPitch(1000.0);
  dryWet = 1.0;
  filterL.setMode(TwoPoleFilter::BANDPASS);
  filterR.setMode(TwoPoleFilter::BANDPASS);
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void WahWah::setSampleRate(double newSampleRate)
{
  ModulationEffect::setSampleRate(newSampleRate);
  filterL.setSampleRate(newSampleRate);
  filterR.setSampleRate(newSampleRate);
}
   
void WahWah::setFilterMode(int newMode)
{
  filterL.setMode(newMode);
  filterR.setMode(newMode);
}
   
void WahWah::setFrequency(double newFrequency)
{
  pitch = freqToPitch(newFrequency);
}

void WahWah::setGain(double newGain)
{
  filterL.setGain(newGain);
  filterR.setGain(newGain);
}
   
void WahWah::setBandwidth(double newBandwidth)
{
  filterL.setBandwidth(newBandwidth);
  filterR.setBandwidth(newBandwidth);
}

//-------------------------------------------------------------------------------------------------
// others:

void WahWah::reset()
{
  filterL.reset();
  filterR.reset();
}

