#include "rosic_BiquadStereoDF2.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

BiquadStereoDF2::BiquadStereoDF2()
{
  reset();                
}

BiquadStereoDF2::~BiquadStereoDF2()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void BiquadStereoDF2::reset()
{
  lw1 = lw2 = 0.0;
  rw1 = rw2 = 0.0;
}

