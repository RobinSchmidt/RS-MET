#ifndef RS_BUFFERFUNCTIONTESTS_H
#define RS_BUFFERFUNCTIONTESTS_H

//#include "../../UnitTestUtilities.h"
#include "../../../Shared/Shared.h"  // new

bool testBufferFunctions(std::string &reportString);


bool testCopySection(   std::string &reportString);

bool testMoveElements(  std::string &reportString);
bool testRemoveElements(std::string &reportString);

#endif
