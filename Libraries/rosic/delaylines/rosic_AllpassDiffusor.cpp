#include "rosic_AllpassDiffusor.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AllpassDiffusor::AllpassDiffusor(int maximumDelayInSamples)
: IntegerDelayLine(maximumDelayInSamples)
{
  g = 0.0;
}

AllpassDiffusor::~AllpassDiffusor()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void AllpassDiffusor::setDiffusionAmount(double newAmount)
{
  g = 0.01 * newAmount;
  g = clip(g, -1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
// inquiry (get-, is-, etc. functions):

double AllpassDiffusor::getDiffusionAmount()
{
  return 100.0 * g;
}
