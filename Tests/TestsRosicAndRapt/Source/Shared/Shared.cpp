#include "Shared.h"

#include "RaptLibraryCode/RaptInstantiations.cpp"

//#include "Plotting/GNUPlotter.cpp"

//#include "Plotting/DSPPlotters.cpp"
//#include "Plotting/Plotting.cpp"
#include "Plotting/rosic_Plotter.cpp" // obsolete - get rid

#include "Prototypes/Prototypes.cpp"
#include "Prototypes/FiniteAutomaton.cpp"
#include "Prototypes/ParticleBouncer.cpp"
#include "Prototypes/ParticleSystem.cpp"
#include "Prototypes/Tensor.cpp"

#include "Utilities/FileWriting.cpp"
#include "Utilities/PerformanceTestTools.cpp"
#include "Utilities/TestInputCreation.cpp"
#include "Utilities/TestUtilities.cpp"

using namespace std;

#undef min

bool areNumbersEqual(double x, double y, double relativeTolerance)
{
  // nan == nan:
  double tmp = RS_NAN(double);
  if( memcmp(&x, &tmp, sizeof(double)) == 0 && memcmp(&y, &tmp, sizeof(double)) == 0 )
    return true;

  // inf == inf:
  tmp = RS_INF(double);
  if( x == tmp && y == tmp )
    return true;

  // -inf == -inf:
  tmp = -tmp;
  if( x == tmp && y == tmp )
    return true;

  // catch case where one or bothe of x, y are denormals - in this case, the 
  // relativeTolerance*rsMax(fabs(x), fabs(y)) yields a zero absolute tolerance whereas the 
  // difference on the left hand side is a very small (denormal) nonzero number
  //double denormThresh = RS_MIN(double);
  double denormThresh = std::numeric_limits<double>::min();
  if( fabs(x) <= denormThresh && fabs(y) <= denormThresh )
    return true;

  // x == y, if the absolute difference is below a relative tolerance:
  return fabs(x-y) <= relativeTolerance * max(fabs(x), fabs(y));
}
