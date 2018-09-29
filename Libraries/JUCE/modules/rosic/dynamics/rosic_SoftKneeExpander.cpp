//#include "rosic_SoftKneeExpander.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

SoftKneeExpander::SoftKneeExpander(int newLookAheadBufferSize) 
: Expander(newLookAheadBufferSize)
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

SoftKneeExpander::~SoftKneeExpander()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void SoftKneeExpander::setKneeWidth(double newKneeWidth)
{
  if( newKneeWidth >= 0.0 )
  {
    kneeWidth = newKneeWidth;
    calculateCoefficients();
  }
}

void SoftKneeExpander::setKneeShape(int newKneeShape)
{
  if( newKneeShape >= CUBIC && newKneeShape <= QUARTIC )
  {
    kneeShape = newKneeShape;
    calculateCoefficients();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void SoftKneeExpander::calculateCoefficients()
{
  if( kneeShape == CUBIC )
    calculateCubicCoefficients();
  else
    calculateQuarticCoefficients();
}

void SoftKneeExpander::calculateCubicCoefficients()
{
  double x1  = threshold - 0.5*kneeWidth;
  double y1  = threshold - 0.5*ratio*kneeWidth;
  double yd1 = ratio;
  double x2  = threshold + 0.5*kneeWidth;
  double y2  = threshold + 0.5*kneeWidth;
  double yd2 = 1.0;
  a4 = 0.0;
  RAPT::fitCubicWithDerivative(x1, x2, y1, y2, yd1, yd2, &a3, &a2, &a1, &a0);
}

void SoftKneeExpander::calculateQuarticCoefficients()
{   
  double t  = threshold;
  double r  = ratio;
  double c  = 0.5*kneeWidth;
  double t2 = t*t;
  double t3 = t*t2;
  double t4 = t2*t2;
  double c2 = c*c;
  double c3 = c*c2;
  double c4 = c2*c2;

  a0  = (r-1)*t4+(6*c2-6*c2*r)*t2+(8*c3-8*c3*r)*t-3*c4*r+3*c4;
  a0 /= 16*c3;
  a1  = (r-1)*t3+(3*c2-3*c2*r)*t-2*c3*r-2*c3;
  a1 /= -4*c3;
  a2  = (3*r-3)*t2-3*c2*r+3*c2;
  a2 /= 8*c3;
  a3  = (r-1)*t;
  a3 /= -4*c3;
  a4  = (r-1)/(16*c3);
}
