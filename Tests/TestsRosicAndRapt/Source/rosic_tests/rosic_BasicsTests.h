#ifndef rosic_BasicsTests_h
#define rosic_BasicsTests_h

//#include "../datastructures/rosic_StringTests.h"
#include "rosic/rosic.h"

namespace rotes
{

  void testBinomialCoefficients();

  void testMathFunctions();
  void testWindowFunctions();




  void testInterpolation();

  void testHermiteTwoPoint1();
  void testHermiteTwoPoint2();
  void testHermiteTwoPoint3();
  void testHermiteTwoPointM();


  rosic::Matrix createHermiteInterpolatorImpulseResponses(int inLength, int oversampling, const int M[5], double shape);
  void plotOneSidedInterpolatorContinuousResponses(int M[5], double shape);
  void plotOneSidedInterpolatorPolyphaseResponses(int M, double shape, double d[5]);


  void testAsymmetricPolynomialInterpolatorsOld();



}

#endif 