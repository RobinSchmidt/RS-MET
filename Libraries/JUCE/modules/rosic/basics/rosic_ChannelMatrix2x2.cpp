#include "rosic_ChannelMatrix2x2.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ChannelMatrix2x2::ChannelMatrix2x2()
{
  gLL = 1.0;
  gRL = 0.0;
  gLR = 0.0;
  gRR = 1.0;
}

ChannelMatrix2x2::~ChannelMatrix2x2()
{

}
