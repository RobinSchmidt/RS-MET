#include "Upsampler.h"

//----------------------------------------------------------------------------
// construction/destruction:

Upsampler::Upsampler()
{
 oversamplingFactor    = 1;
 oversamplingFactorRec = 1.0;
 previousInput         = 0.0;
}

Upsampler::~Upsampler()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void Upsampler::setOversamplingFactor(int newOversamplingFactor)
{
 if( newOversamplingFactor >= 1 )
  oversamplingFactor = newOversamplingFactor;

 oversamplingFactorRec = 1.0 / oversamplingFactor;
}





//----------------------------------------------------------------------------
//others:



void Upsampler::reset()
{
 previousInput = 0.0;
}
