#ifndef RAPT_MATHUNITTESTS_H
#define RAPT_MATHUNITTESTS_H

//#include "../../Shared/Shared.h"
#include "Math/LinearAlgebraUnitTests.h"
#include "Math/PolynomialUnitTests.h"

#include "Math/VectorUnitTests.h"
//#include "Math/MatrixUnitTests.h"
//#include "Math/MiscMathUnitTests.h"

bool coordinateMapperUnitTest();
bool interpolatingFunctionUnitTest();
bool rootFinderUnitTest();

bool polynomialRootsUnitTest(); // the new explicit formulas - move to PolynomialUnitTests

#endif