//#include "rosic_HelperFunctions.h"
//using namespace rosic;

// there's a similar (or even same?) function in rapt - so get rid of this:
double rosic::findNearestUpwardZeroCrossing(float* buffer, int length, double searchStart)
{
  // initialize the search pointers for above and below the search start
  int hi1 = (int) ceil(searchStart);
  int hi2 = hi1 + 1;
  int lo1 = (int) floor(searchStart);
  int lo2 = lo1 - 1;

  // catch some special conditions:
  if( lo2 < 0 || hi2 >= length )
    return searchStart;

  // check whether the serch start itself is in between two values which surround a zero 
  // crossing:
  if( buffer[lo1] <= 0.0 && buffer[hi1] > 0.0 )
  {
    return lo1 + buffer[lo1] / (buffer[lo1]-buffer[hi1]);
  }

  // search to the right (ascending direction):
  while( hi2 < length  )
  {
    if( buffer[hi1] <= 0.0 && buffer[hi2] > 0.0 )
      break;
    hi1++;
    hi2++;
  }

  // search to the left (descending direction):
  while( lo2 >= 0  )
  {
    if( buffer[lo2] <= 0.0 && buffer[lo1] > 0.0 )
      break;
    lo1--;
    lo2--;
  }

  // if the search was without success, we return the stating point of the search:
  if( lo2 == 0 && hi2 == length-1 )
    return searchStart;

  // select the closer from both and return the zero crossing of the connecting line:
  if( (searchStart-lo2) <= -(searchStart-hi1) )
    return lo2 + buffer[lo2] / (buffer[lo2]-buffer[lo1]);
  else
    return hi1 + buffer[hi1] / (buffer[hi1]-buffer[hi2]);
}

void rosic::error(const char *message)
{
  DEBUG_BREAK;
  printf("%s %s", message, "\n");
}



