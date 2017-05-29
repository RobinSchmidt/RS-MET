//#include "rosic_CombStereoizer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombStereoizer::CombStereoizer()
{
  setDryWetRatio(0.5);
  channelSwitch = false;
  delayLine.setDelayInMilliseconds(20.0);
}

CombStereoizer::~CombStereoizer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void CombStereoizer::setSampleRate(double newSampleRate)
{
  delayLine.setSampleRate(newSampleRate);
  wetFilter.setSampleRate(newSampleRate);
}

/*
void CombStereoizer::setDryWet(double newDryWet)
{
  if( newDryWet >= 0.0 && newDryWet <= 100.0 )
  {
    wetFactor = 0.01*newDryWet;
    dryFactor = 1.0 - wetFactor;
  }
}
*/

void CombStereoizer::setSwitchWetLeftForRight(bool shouldSwitch)
{
  channelSwitch = shouldSwitch;
}

//-------------------------------------------------------------------------------------------------
// inquiry:
/*
double CombStereoizer::getDryWet()
{
  return 100.0*wetFactor;
}
*/

bool CombStereoizer::isWetSignalChannelSwitched()
{
  return channelSwitch;
}

//-------------------------------------------------------------------------------------------------
// others:

void CombStereoizer::reset()
{
  delayLine.clearDelayBuffer();
  wetFilter.resetBuffers();
}