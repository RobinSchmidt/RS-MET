#include "rosic_StereoWidth.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StereoWidth::StereoWidth()
{
  // initialize parameters:
  setMidGain(0.0);
  setSideGain(0.0);
  setGlobalGain(0.0);
  mixToMono = false;
}

StereoWidth::~StereoWidth()
{

}


