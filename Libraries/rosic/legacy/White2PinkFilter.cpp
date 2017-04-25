#include "White2PinkFilter.h"

White2PinkFilter::White2PinkFilter()
{
 reset();
}

White2PinkFilter::~White2PinkFilter()
{

}

void White2PinkFilter::reset()
{
 for(int i=0; i<7; i++)
  b[i] = 0.0;
}