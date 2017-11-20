#ifndef RS_ARRAYTESTS_H
#define RS_ARRAYTESTS_H

#include "../../UnitTestUtilities.h"

bool testArray(std::string &reportString);

bool testArrayAppend(std::string &reportString);
bool testArrayInsert(std::string &reportString);
bool testArrayRemove(std::string &reportString);
bool testArrayGrowAndShrink(std::string &reportString);
bool testArrayMisc(std::string &reportString);

bool testFlagArray(std::string &reportString);

bool testInfiniteDataStream(std::string &reportString);
 // this should either be in a separate file or (better) we consolidate all container tests into
 // a single file

#endif
