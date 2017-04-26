#include "rosic_SoftKneeCompressor.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SoftKneeCompressor::SoftKneeCompressor(int newLookAheadBufferSize) 
: Compressor(newLookAheadBufferSize)
{
  kneeWidth = 0.0;
  a0        = 0.0;
  a1        = 1.0;  // coeff for x^1
  a2        = 0.0;
  a3        = 0.0;
  a4        = 0.0;
  kneeShape = QUARTIC;
  antiAlias = true;
}

SoftKneeCompressor::~SoftKneeCompressor()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void SoftKneeCompressor::setKneeWidth(double newKneeWidth)
{
  if( newKneeWidth >= 0.0 )
  {
    kneeWidth = newKneeWidth;
    calculateCoefficients();
  }
}

void SoftKneeCompressor::setKneeShape(int newKneeShape)
{
  if( newKneeShape >= CUBIC && newKneeShape <= QUARTIC )
  {
    kneeShape = newKneeShape;
    calculateCoefficients();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SoftKneeCompressor::calculateCoefficients()
{
  if( kneeShape == CUBIC )
    calculateCubicCoefficients();
  else
    calculateQuarticCoefficients();
  updateAutoGainFactor();
}

void SoftKneeCompressor::calculateCubicCoefficients()
{
  double x1  = threshold - 0.5*kneeWidth;
  double y1  = threshold - 0.5*kneeWidth;
  double yd1 = 1.0;
  double x2  = threshold + 0.5*kneeWidth;
  double y2  = threshold + (1.0/ratio)*(0.5*kneeWidth);
  double yd2 = 1.0/ratio;
  if( limiterMode )
  {
    y2  = threshold;
    yd2 = 0.0;
  }
  a4 = 0.0;
  fitCubicWithDerivative(x1, x2, y1, y2, yd1, yd2, &a3, &a2, &a1, &a0);
}

void SoftKneeCompressor::calculateQuarticCoefficients()
{   
  double t = threshold;
  double r = ratio;
  double c = 0.5*kneeWidth;

  double t2 = t*t;
  double t3 = t*t2;
  double t4 = t2*t2;
  double c2 = c*c;
  double c3 = c*c2;
  double c4 = c2*c2;

  if( r == INF || limiterMode == true )
  {
    a0 = (t4 - 6*c2*t2 + 8*c3*t - 3*c4) / (16*c3);
    a1 = -(t3 - 3*c2*t - 2*c3) / (4*c3);
    a2 = (3*t2 - 3*c2) / (8*c3);
    a3 = -t / (4*c3);
    a4 = 1  / (16*c3);
  }
  else
  {
    a0 = (r-1)*t4 + c2*(6-6*r)*t2 + c3*(8*r-8)*t + c4*(3-3*r);
    a0 = a0 / (16*c3*r);
    a1 = (r-1)*t3 + c2*(3-3*r)*t + c3*(-2*r-2);
    a1 = -a1 / (4*c3*r);
    a2 = (3*r-3)*t2 + c2*(3-3*r);
    a2 = a2 / (8*c3*r);
    a3 = -(r-1)*t / (4*c3*r);
    a4 = (r-1) / (16*c3*r);
  }
}
