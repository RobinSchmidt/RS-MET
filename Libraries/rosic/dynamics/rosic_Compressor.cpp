#include "rosic_Compressor.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Compressor::Compressor(int newLookAheadBufferSize) : DynamicsProcessorBase(newLookAheadBufferSize)
{
  threshold      = 0.0;
  ratio          = 1.0;
  oneOverRatio   = 1.0;
  autoGainFactor = 1.0;
  autoGainActive = false;
  limiterMode    = false;
}

Compressor::~Compressor()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void Compressor::updateAutoGainFactor()
{
  if( autoGainActive == true )
    autoGainFactor = dB2amp(-transferCurveAt(0.0));
  else
    autoGainFactor = 1.0;
}
