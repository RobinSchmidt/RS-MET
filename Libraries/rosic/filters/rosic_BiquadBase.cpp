#include "rosic_BiquadBase.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BiquadBase::BiquadBase()
{
  initializeCoefficients();             
}

BiquadBase::~BiquadBase()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void BiquadBase::setCoefficients(double newB0, double newB1, double newB2, double newA1, double newA2)
{
  b0 = newB0;
  b1 = newB1;
  b2 = newB2;
  a1 = newA1;
  a2 = newA2;
}

//-------------------------------------------------------------------------------------------------
// others:

void BiquadBase::initializeCoefficients()
{
  setCoefficients(1.0, 0.0, 0.0, 0.0, 0.0);
}


