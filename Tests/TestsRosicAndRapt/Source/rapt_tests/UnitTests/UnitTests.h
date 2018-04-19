#ifndef UNITTESTS_INCLUDED
#define UNITTESTS_INCLUDED

#include "DataUnitTests.h"
#include "MathUnitTests.h"
#include "ImageUnitTests.h"
#include "FilterUnitTests.h"

/*
// just to test some compiler error:
template<class T>
class Foo
{
public:
  template<class S>
  static void bar(S x);
};
*/

bool runAllUnitTests();

#endif