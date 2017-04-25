#ifndef rosic_RealFunctionEvaluationAlgorithms_h
#define rosic_RealFunctionEvaluationAlgorithms_h

// standard includes:
#include <math.h>

// rosic includes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

  /** Performs M iterations of the Landen transformation of an elliptic modulus k and returns the 
  results in the array v which must be of length M. */
  void landen(double k, int M, double* v);


} // end namespace rosic

#endif // #ifndef rosic_RealFunctionEvaluationAlgorithms_h