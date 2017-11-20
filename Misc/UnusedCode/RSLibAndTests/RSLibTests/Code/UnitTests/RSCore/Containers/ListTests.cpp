#include "ListTests.h"

/*
#include "../../HelperFunctions/TestHelperFunctions.h"

#include "../../../../Libraries/RSCore/DataStructures/List.h"
using namespace RSCore;
*/

bool testList(std::string &reportString)
{
  std::string testName = "rsList";
  bool testResult = true;

  testResult &= testListInsert(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testListInsert(std::string &reportString)
{
  std::string testName = "rsListInsert";
  bool testResult = true;

  rsList<int> testList;
  testList.appendItem(1);
  testList.appendItem(2);
  testList.appendItem(3);
  testList.prependItem(0);
  testList.insertSorted(2, false);
  testList.insertSorted(2, true);

  testList.insertSorted(-1, true);
  testList.insertSorted(-2, false);


  rsArray<int> listAsArray = testList.getAsArray();

  int cArray[8] = {-2,-1,0,1,2,2,2,3};
  rsArray<int> targetArray(cArray, 8);

  testResult &= ( listAsArray == targetArray );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
