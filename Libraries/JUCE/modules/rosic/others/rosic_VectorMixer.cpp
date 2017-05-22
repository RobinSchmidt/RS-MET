//#include "rosic_VectorMixer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

VectorMixer::VectorMixer()
{
  x = 0.0;
  y = 0.0;
}

VectorMixer::~VectorMixer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void VectorMixer::setX(double newX)
{
  x = newX;
}

void VectorMixer::setY(double newY)
{
  y = newY;
}