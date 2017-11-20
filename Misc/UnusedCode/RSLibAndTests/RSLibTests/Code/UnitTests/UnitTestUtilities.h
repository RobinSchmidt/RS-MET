#ifndef RS_UNITTESTUTILITIES_H
#define RS_UNITTESTUTILITIES_H

#include "../Common/TestUtilities.h"


/** Comparison function that compares with a given error tolerance and also returns true when the 
involved numbers are NaNs or infinities. */
bool areNumbersEqual(double x, double y, double relativeTolerance);

void appendTestHeaderToReport(std::string &reportString, const std::string &nameOfTest);
void appendTestResultToReport(std::string &reportString, const std::string &nameOfTest, 
                              const bool result);

#endif
