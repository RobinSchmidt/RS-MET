#ifndef RS_MISCAUDIOTESTS_H
#define RS_MISCAUDIOTESTS_H

#include "../../UnitTestUtilities.h"

bool testBandwidthConversions(std::string &reportString);
bool testSincInterpolation(   std::string &reportString);
bool testSineParameters(      std::string &reportString);
bool testZeroCrossingFinder(  std::string &reportString);

#endif
