//#ifndef RS_MISCAUDIOTESTS_H
//#define RS_MISCAUDIOTESTS_H

#pragma once

// check, if and why we still need this - and get rid

bool testBandwidthConversions(std::string &reportString);
bool testSincInterpolation(   std::string &reportString);
bool testSineParameters(      std::string &reportString);
bool testZeroCrossingFinder(  std::string &reportString);

//#endif
