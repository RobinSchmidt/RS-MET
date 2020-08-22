//#include "romos_Interpolation.h"
//#include <math.h>
//using namespace romos;

double interpolateExponentially(double x, double xL, double yL, double xR, double yR, double shape)
{
  double p = (x-xL) / (xR-xL);  // proportion of passedLength/fullLength (Eq.3-30)
  if( shape == 0.0 )
    return yL + (yR-yL) * p;    // needs to treated separately to avoid division-by-zero
  else
    return yL + (yR-yL) * (1.0 - exp(p*shape)) / (1.0 - exp(shape));
}
// this formula is also implemented in RAPT::rsNodeBasedFunction. factor it out there and use
// it here - also with the mapped shape parameter