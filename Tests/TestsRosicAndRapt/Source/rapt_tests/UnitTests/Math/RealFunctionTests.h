#ifndef RS_REALFUNCTIONTESTS_H
#define RS_REALFUNCTIONTESTS_H

#include "../../../Shared/Shared.h"

bool testRealFunctions();

bool testAbsAndSign(std::string &reportString);
bool testHyperbolicFunctions(std::string &reportString);
bool testSinc(std::string &reportString);
bool testFunctionIterators(std::string &reportString);


#endif
