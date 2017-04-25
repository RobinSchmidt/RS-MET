#include "rosic_WhiteToPinkFilter.h"
using namespace rosic;

WhiteToPinkFilter::WhiteToPinkFilter()
{
 reset();
}

WhiteToPinkFilter::~WhiteToPinkFilter()
{

}

void WhiteToPinkFilter::reset()
{
 for(int i=0; i<7; i++)
  b[i] = 0.0;
}