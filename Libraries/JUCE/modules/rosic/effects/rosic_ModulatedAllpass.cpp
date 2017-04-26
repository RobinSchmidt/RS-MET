#include "rosic_ModulatedAllpass.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModulatedAllpass::ModulatedAllpass() 
{
  factor = 0.0;
  offset = 0.0;
  reset();
}

ModulatedAllpass::~ModulatedAllpass()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void ModulatedAllpass::setSampleRate(double newSampleRate)
{

}


//-------------------------------------------------------------------------------------------------
// others:

void ModulatedAllpass::reset()
{
  wL = 0.0;
  wR = 0.0;
}

