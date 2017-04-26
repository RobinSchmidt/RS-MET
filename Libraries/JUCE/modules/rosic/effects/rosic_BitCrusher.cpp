#include "rosic_BitCrusher.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BitCrusher::BitCrusher()
{
  decimationFactor     = 1;
  quantizationInterval = 0.00001f;
  amount               = 1.0;
  reset();
}

BitCrusher::~BitCrusher()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void BitCrusher::reset()
{
  sampleCounter        = 0;
  outSampleL           = 0.0;
  outSampleR           = 0.0;
}