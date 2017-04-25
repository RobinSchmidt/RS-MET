#include "rosic_PhaseStereoizer.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PhaseStereoizer::PhaseStereoizer()
{
  setDryWetRatio(0.5);
  channelSwitch = false;
  sinFactor     = 1.0;
  cosFactor     = 0.0;
  setDryWetRatio(1.0);
  setMidSideRatio(0.5);
}

PhaseStereoizer::~PhaseStereoizer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void PhaseStereoizer::setSampleRate(double newSampleRate)
{
  filterM.setSampleRate(newSampleRate);
  filterS.setSampleRate(newSampleRate);
}

void PhaseStereoizer::setPhaseOffset(double newPhaseOffset)
{
  double phi = degreeToRadiant(newPhaseOffset);
  sinFactor  = sin(phi);
  cosFactor  = cos(phi);
}

//-------------------------------------------------------------------------------------------------
// others:

void PhaseStereoizer::reset()
{
  quadratureNetwork.reset();
  filterM.reset();
  filterS.reset();
}
