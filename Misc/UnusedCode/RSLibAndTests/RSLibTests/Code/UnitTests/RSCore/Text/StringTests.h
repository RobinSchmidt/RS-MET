#ifndef RS_STRINGTESTS_H
#define RS_STRINGTESTS_H

#include "../../UnitTestUtilities.h"


bool testString(std::string &reportString);

bool testCharacterComparisons(std::string &reportString);
bool testStringComparisons(std::string &reportString);
bool testStringBufferCopying(std::string &reportString);
bool testStringCharacterRemoval(std::string &reportString);
bool testStringIntConversions(std::string &reportString);
bool testStringDoubleConversions(std::string &reportString);
bool testSubStringFunctions(std::string &reportString);

bool testStringDoubleConversionsRandom(std::string &reportString);
bool testStringDoubleConversionsSpecialValues(std::string &reportString);
bool testStringDoubleConversionsDenormals(std::string &reportString);
bool testStringDoubleConversionsLarge(std::string &reportString);
bool testStringDoubleConversionsGeometricProgression(std::string &reportString, double start, double factor);

//void testStringConcatenation();
//void testStringComparison();...

rsString createStringWithAllCharacters();
rsString createStringWithAllPrintableCharacters();


#endif
