#include "rosic_DualTwoPoleFilter.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DualTwoPoleFilter::DualTwoPoleFilter()
{
  freq1     = 500.f;
  freq2     = 2000.f;
  freqScale = 1.f;

  gain1 = gain2 = 0.f;
  gainScale = 1.f;

  bandwidth1 = bandwidth2 = 1.f;
  bandwidthScale = 1.f;

  gp = 0.0;
}

DualTwoPoleFilter::~DualTwoPoleFilter()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void DualTwoPoleFilter::setSerialParallelBlend(double newBlend)
{
  if( newBlend >= 0.0 && newBlend <= 1.0 )
    gp = newBlend;
}








//-------------------------------------------------------------------------------------------------
// inquiry:





//-------------------------------------------------------------------------------------------------
// others:





