#include "MidSideMixer.h"

//----------------------------------------------------------------------------
// construction/destruction:

MidSideMixer::MidSideMixer()
{
 // initialize parameters:
 setMidGain(0.0);
 setSideGain(0.0);
 setGlobalGain(0.0);
}

MidSideMixer::~MidSideMixer()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void MidSideMixer::setMidGain(double newMidGain)
{
 midGain = dB2amp(newMidGain);
}

void MidSideMixer::setSideGain(double newSideGain)
{
 sideGain = dB2amp(newSideGain);
}

void MidSideMixer::setGlobalGain(double newGlobalGain)
{
 globalGain = dB2amp(newGlobalGain);
}
