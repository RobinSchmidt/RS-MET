#include "rosic_Limiter.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Limiter::Limiter(int newLookAheadBufferSize) : DynamicsProcessorBase(newLookAheadBufferSize)
{
  limit = 1.0;
  setAttackTime(0.0);    
  setReleaseTime(10.0); 
}

Limiter::~Limiter()
{

}
