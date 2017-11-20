#include "UnitTestUtilities.h"

void appendTestHeaderToReport(std::string &reportString, const std::string &nameOfTest)
{
  reportString += "Testing " + nameOfTest + "\n";
}

void appendTestResultToReport(std::string &reportString, const std::string &nameOfTest, 
                              bool result)
{
  if( result == true )
    reportString += nameOfTest + ": OK \n";
  else
    reportString += nameOfTest + ": !!! FAILED !!!\n";
}

// a comparison function that compares with a given error tolerance and also returns true also when
// the involved numbers are NaNs:
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
  double denormThresh = RS_MIN(double);
  if( fabs(x) <= denormThresh && fabs(y) <= denormThresh )
    return true;

  // x == y, if the absolute difference is below a relative tolerance:
  return fabs(x-y) <= relativeTolerance*rsMax(fabs(x), fabs(y));
}

