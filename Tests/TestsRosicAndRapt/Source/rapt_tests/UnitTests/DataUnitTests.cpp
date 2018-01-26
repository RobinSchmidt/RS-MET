#include "DataUnitTests.h"
using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

//double sum(double* a, size_t N)
//{
//  double accu = 0;
//  for(size_t i = 0; i < N; i++)
//    accu += a[i];
//  return accu;
//}

bool float64x2UnitTest()
{
  bool r = true;      // test result

  // test constructors and getters:
  //rsFloat64x2 x00;              // default constructor, two zeros
  //r &= x00.get0() == 0.0;       // ...or actually no, we leave them uninitialized
  //r &= x00.get1() == 0.0;

  // construct from a double:
  rsFloat64x2 x11(1.0);  r &= x11.get0() == 1.0; r &= x11.get1() == 1.0;

  // construct from two doubles:
  rsFloat64x2 x12(1.0, 2.0); r &= x12.get0() == 1.0; r &= x12.get1() == 2.0;

  // construct from array of doubles:
  double arr[2] = { 3.0, 4.0 };
  rsFloat64x2 x34(arr); r &= x34.get0() == 3.0; r &= x34.get1() == 4.0;

  // construct from another instance:
  rsFloat64x2 y(x34); r &= y.get0() == 3.0; r &= y.get1() == 4.0;

  // test setters:
  y.set0(5.0);     r &= y.get0() == 5.0; r &= y.get1() == 4.0;
  y.set1(6.0);     r &= y.get0() == 5.0; r &= y.get1() == 6.0;
  y.set(1.0, 2.0); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y.set(3.0);      r &= y.get0() == 3.0; r &= y.get1() == 3.0;

  // assignment and equality:
  y = x12; r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  //r &= y == x12;

  // binary arithmetic operators:
  y = x12 + x34; r &= y.get0() == 4.0; r &= y.get1() == 6.0;
  y = x34 - x12; r &= y.get0() == 2.0; r &= y.get1() == 2.0;
  y = x12 * x34; r &= y.get0() == 3.0; r &= y.get1() == 8.0;
  y = x34 / x12; r &= y.get0() == 3.0; r &= y.get1() == 2.0;

  // binary arithmetic operators with scalar lhs:
  y  = x12;
  y  = 1.0  + y; r &= y.get0() == 2.0; r &= y.get1() == 3.0;
  y  = 5.0  - y; r &= y.get0() == 3.0; r &= y.get1() == 2.0;
  y  = 2.0  * y; r &= y.get0() == 6.0; r &= y.get1() == 4.0;
  y  = 12.0 / y; r &= y.get0() == 2.0; r &= y.get1() == 3.0;

  // binary arithmetic operators with scalar rhs:
  y  = x34;
  y  = y + 2.0; r &= y.get0() == 5.0; r &= y.get1() == 6.0;
  y  = y - 2.0; r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y  = y * 2.0; r &= y.get0() == 6.0; r &= y.get1() == 8.0;
  y  = y / 2.0; r &= y.get0() == 3.0; r &= y.get1() == 4.0;

  // unary arithmetic operators -,+=,... :
  y  = x12;
  y += x34; r &= y.get0() ==  4.0; r &= y.get1() ==   6.0;
  y -= x12; r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;
  y *= x34; r &= y.get0() ==  9.0; r &= y.get1() ==  16.0;
  y /= x34; r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;
  y  = -y;  r &= y.get0() == -3.0; r &= y.get1() ==  -4.0;
  y  = -y;  r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;

  // functions: sqrt, min, max, clip, abs, sign:
  y = x34 * x34;
  y = rsSqrt(y);       r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y = rsMin(x12, x34); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y = rsMin(x34, x12); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y = rsMax(x12, x34); r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y = rsMax(x34, x12); r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y.set(1.0, 9.0);
  y = rsClip(y, 2.0, 7.0); r &= y.get0() == 2.0; r &= y.get1() == 7.0; 
  y.set(-2.0, -3.0);
  y = rsAbs(y); r &= y.get0() == 2.0; r &= y.get1() == 3.0;
  y.set(-2.0, 3.0);
  y = rsSign(y); r &= y.get0() == -1.0; r &= y.get1() == +1.0;

  // reductions to scalar:
  double s;
  y.set(3.0, 5.0);
  s = y.getSum(); r &= s == 8.0;
  s = y.getMin(); r &= s == 3.0;
  s = y.getMax(); r &= s == 5.0;

  return r;
}

bool complexFloat64x2UnitTest()
{  
  bool r = true;      // test result

  // we have 4 complex numbers z1[0] = 1 + 3i, z1[1] = 2 + 4i, z2[0] = 5 + 7i, z2[1] = 6 + 8i:
  std::complex<double> z10(1, 3), z11(2, 4), z20(5, 7), z21(6, 8), w0, w1;
  rsFloat64x2 re1(1, 2), im1(3, 4), re2(5, 6), im2(7, 8);
  std::complex<rsFloat64x2> z1(re1, im1), z2(re2, im2), w;
  //std::complex<rsFloat64x2> w; // for outputs
  //std::complex<double> w0, w1;  

  // addition:
  w = z1 + z2;
  r &= w.real().get0() ==  6;
  r &= w.imag().get0() == 10;
  r &= w.real().get1() ==  8;
  r &= w.imag().get1() == 12;

  // division:
  w = z2 / z2;
  r &= w.real().get0() == 1;
  r &= w.imag().get0() == 0;
  r &= w.real().get1() == 1;
  r &= w.imag().get1() == 0;
  w   = z2;
  w  /= z2;
  r &= w.real().get0() == 1;
  r &= w.imag().get0() == 0;
  r &= w.real().get1() == 1;
  r &= w.imag().get1() == 0;

  //z = std::exp(z1); // this doesn't work - it doesn't try to invoke the exp for rsFloat64x2
  // i think, we need to implement explicit specializations for the math functions for
  // complex<rsFloat64x2>

  // exponential function:
  w0 = std::exp(z10);
  w1 = std::exp(z11);
  w  = rsExp(z1);
  r &= w0 == get0(w);
  r &= w1 == get1(w);

  return r;
}