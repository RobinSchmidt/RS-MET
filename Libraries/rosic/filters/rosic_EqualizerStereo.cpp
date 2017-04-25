#include "rosic_EqualizerStereo.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

EqualizerStereo::EqualizerStereo()
{
  stereoMode = STEREO_LINKED;
  bypass     = false;
}

EqualizerStereo::~EqualizerStereo()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

void EqualizerStereo::getMagnitudeResponse(int channel, double *frequencies, double *magnitudes, int numBins)
{
  equalizers[channel].getMagnitudeResponse(frequencies, magnitudes, numBins);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void EqualizerStereo::reset()
{
  equalizers[0].reset();
  equalizers[1].reset();
}
