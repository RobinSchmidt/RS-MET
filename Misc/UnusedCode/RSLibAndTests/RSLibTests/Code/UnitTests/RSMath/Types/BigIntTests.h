#ifndef RS_BIGINTTESTS_H
#define RS_BIGINTTESTS_H
 // rename to BigNumberTests

#include "../../UnitTestUtilities.h"

bool testArbitraryPrecision(std::string &reportString); // rename to testBigNumbers

bool testDigitArithmetic(std::string &reportString);
bool testBaseChange(std::string &reportString);


bool testBigInt(std::string &reportString);

bool testBigFloat(std::string &reportString);


#endif
