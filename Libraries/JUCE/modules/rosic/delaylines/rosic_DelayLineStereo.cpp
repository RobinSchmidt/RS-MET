#include "rosic_DelayLineStereo.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

DelayLineStereo::DelayLineStereo(int bufferLengthToAllocate)
{
  if( bufferLengthToAllocate < 512 )
    DEBUG_BREAK;  // you need allocate some reasonable amount of memory

  length  = bufferLengthToAllocate;
  bufferL = new(std::nothrow) double[length+interpolatorMargin+1];
  bufferR = new(std::nothrow) double[length+interpolatorMargin+1]; 
  tapIn   = 0;
  clearBuffers();
}

DelayLineStereo::~DelayLineStereo()
{
  delete[] bufferL;
  delete[] bufferR;
}

//-------------------------------------------------------------------------------------------------
// others:

void DelayLineStereo::clearBuffers()
{
  for(int i=0; i<length+interpolatorMargin; i++)
  {
    bufferL[i] = 0.0;
    bufferR[i] = 0.0;
  }
}