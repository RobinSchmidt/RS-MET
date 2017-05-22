//#include "rosic_WarpedAllpassInterpolator.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WarpedAllpassInterpolator::WarpedAllpassInterpolator()
{
  coeff          = 0.0;
  previousOutput = 0.0;
}

WarpedAllpassInterpolator::~WarpedAllpassInterpolator()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void WarpedAllpassInterpolator::reset()
{
  previousOutput = 0.0;
}