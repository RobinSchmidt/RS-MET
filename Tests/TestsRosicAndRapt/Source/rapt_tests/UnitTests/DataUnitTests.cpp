#include "DataUnitTests.h"
using namespace RAPT;

bool float64x2UnitTest()
{
  bool r = true;      // test result

  // test constructors and getters:
  rsFloat64x2 x00;              // default constructor, two zeros
  r &= x00.get0() == 0.0;
  r &= x00.get1() == 0.0;

  rsFloat64x2 x11(1.0);         // construct from a double
  r &= x11.get0() == 1.0;
  r &= x11.get1() == 1.0;

  rsFloat64x2 x12(1.0, 2.0);    // construct from two doubles
  r &= x12.get0() == 1.0;
  r &= x12.get1() == 2.0;

  double arr[2] = { 3.0, 4.0 }; // construct from array of doubles
  rsFloat64x2 x34(arr);
  r &= x34.get0() == 3.0;
  r &= x34.get1() == 4.0;

  rsFloat64x2 y(x34);           // construct from another instance
  r &= y.get0() == 3.0;
  r &= y.get1() == 4.0;

  // test setters:
  //y.set0(5.0);
  //r &= y.get0() == 5.0;
  //r &= y.get1() == 4.0;
  //y.set1(6.0);
  //r &= y.get0() == 5.0;
  //r &= y.get1() == 6.0;
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



	  
  return r;
}