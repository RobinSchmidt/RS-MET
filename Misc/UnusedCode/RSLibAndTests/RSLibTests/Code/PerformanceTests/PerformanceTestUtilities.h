#ifndef RS_PERFORMANCETESTUTILITIES_H
#define RS_PERFORMANCETESTUTILITIES_H

#include "../Common/TestUtilities.h"
#include "ProcessorCycleCounter.h"

std::string createPerformanceTestHeader();


void appendResultToReport(std::string &reportString, const std::string &nameOfTest, 
                          const double result);


/*
void appendTestHeaderToReport(std::string &reportString, const std::string &nameOfTest);

*/

#endif
