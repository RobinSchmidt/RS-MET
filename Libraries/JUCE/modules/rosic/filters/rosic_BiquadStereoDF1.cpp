#include "rosic_BiquadStereoDF1.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BiquadStereoDF1::BiquadStereoDF1()
{
  reset();                
}

BiquadStereoDF1::~BiquadStereoDF1()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void BiquadStereoDF1::reset()
{
  lx1 = lx2 = ly1 = ly2 = 0.0;
  rx1 = rx2 = ry1 = ry2 = 0.0;
}

