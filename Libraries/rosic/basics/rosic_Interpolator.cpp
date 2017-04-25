#include "rosic_Interpolator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Interpolator::Interpolator()
{
  previousOutput      = 0;
  interpolationMethod = LINEAR;
}

Interpolator::~Interpolator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Interpolator::setInterpolationMethod(int newMethod)
{
  interpolationMethod = newMethod;
}

//-------------------------------------------------------------------------------------------------
// others:

void Interpolator::reset()
{
  previousOutput = 0;
}