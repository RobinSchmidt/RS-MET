#include "rosic_EnvelopeFollower.h"
using namespace rosic;

// Construction/Destruction:

EnvelopeFollower::EnvelopeFollower()
{
  mode = 1; // mean absolute value
}

EnvelopeFollower::~EnvelopeFollower()
{

}

// Setup:

void EnvelopeFollower::setMode(int Mode)
{
  if( Mode >= MEAN_ABS && Mode < NUM_MODES )
    mode = Mode;
}