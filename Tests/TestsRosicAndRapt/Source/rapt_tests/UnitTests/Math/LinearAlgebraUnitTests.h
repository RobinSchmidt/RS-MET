#ifndef RAPT_LINEARALGEBRATESTS_H
#define RAPT_LINEARALGEBRATESTS_H

#include "../../../Shared/Shared.h"

//bool testLinearAlgebra(std::string &reportString);
bool testLinearAlgebra();


bool testBandDiagonalSolver(std::string &reportString);

bool testLinearSystem2x2(std::string &reportString);
bool testLinearSystem3x3(std::string &reportString);
bool testLinearSystemViaGauss(std::string &reportString);
bool testGaussJordanInversion(std::string &reportString);
bool testTridiagonalSystem(std::string &reportString);

// \todo move these two into the test module for RSCore
bool testSquareMatrixTranspose(std::string &reportString);
bool testMatrixVectorMultiply(std::string &reportString);
bool testMatrixMultiply(std::string &reportString);

bool testChangeOfBasis(std::string &reportString);



#endif
