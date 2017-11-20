#ifndef RS_KEYVALUEMAPTESTS_H
#define RS_KEYVALUEMAPTESTS_H

#include "../../UnitTestUtilities.h"

// we need to have these functions in the RSLib namespace in order to make the friend declaration 
// in the to-be-tested classes work (try to get rid of this):
namespace RSLib
{

  bool testKeyValueMap(std::string &reportString);

  bool testKeyValueMapInsert(std::string &reportString);
  //bool testKeyValueMapRemove(std::string &reportString);
  bool testKeyValueMapFind(std::string &reportString);

}

#endif
