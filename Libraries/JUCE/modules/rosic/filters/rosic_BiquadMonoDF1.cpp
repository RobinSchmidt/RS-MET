#include "rosic_BiquadMonoDF1.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BiquadMonoDF1::BiquadMonoDF1()
{
  reset();                
}

BiquadMonoDF1::~BiquadMonoDF1()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void BiquadMonoDF1::reset()
{
  x1 = x2 = y1 = y2 = 0.0;
}

