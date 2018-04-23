#ifndef RS_TRANSFORMSTESTS_H
#define RS_TRANSFORMSTESTS_H

#include "../../../Shared/Shared.h"

bool testTransforms();

bool testSmbFFT(std::string &reportString);
bool testRsFFT( std::string &reportString);
bool testFourierTransformerRadix2(std::string &reportString);



bool testCorrelation(std::string &reportString);
 // maybe move to testStatistics or something

bool testFitQuadratic(std::string &reportString);
bool testFitQuarticWithDerivatives(std::string &reportString);
 // maybe move to testCurveFitting


//bool testLinearSystem3x3(std::string &reportString);

#endif
