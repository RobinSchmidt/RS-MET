#include "Utilities.h"

int min(int x, int y)
{
  if( x < y )
    return x;
  else
    return y;
}

int max(int x, int y)
{
  if( x > y )
    return x;
  else
    return y;
}

double pitchOffsetToFreqFactor(double pitchOffset)  
{
  return pow(2.0, pitchOffset/12.0);
}

void clearBuffer(float *buffer, int length)
{
  memset(buffer, 0, length*sizeof(float));
}

void copyBuffer(float *source, float *destination, int length)
{
  memcpy(destination, source, length*sizeof(float));
}

double random(double min, double max, int seed)
{
  static unsigned int state = 0;
  if(seed >= 0)
    state = seed;
  state = 1664525 * state + 1013904223;
  return min + (max - min) * ((1.0 / 4294967296.0) * state);
}

double cubicFadeIn(double x)
{
  return x*(x*((PI/2-2)*x+(3-PI))+(PI/2));
}

double cubicFadeOut(double x)
{
  return x*x*((2-PI/2)*x+(PI/2-3))+1;
}
