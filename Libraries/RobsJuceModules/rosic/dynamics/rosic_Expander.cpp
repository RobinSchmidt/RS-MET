//#include "rosic_Expander.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Expander::Expander(int newLookAheadBufferSize) : DynamicsProcessorBase(newLookAheadBufferSize)
{
  threshold = 0.0;
  ratio     = 1.0;
  gateMode  = false;
}

Expander::~Expander()
{

}
