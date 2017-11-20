#ifndef RSCORE_CALLBACKTESTS_H
#define RSCORE_CALLBACKTESTS_H

#include "../../UnitTestUtilities.h"

bool testCallbacks(std::string &reportString);

bool testCallbackByValueSemantics(std::string &reportString);
bool testCallbackByReferenceSemantics(std::string &reportString);


bool testParameter(std::string &reportString);


//bool testCallback2(std::string &reportString);

// \todo add tests with classes that use various forms of inheritance (single, multi, 
// single-virtual, multi-virtual)

#endif
