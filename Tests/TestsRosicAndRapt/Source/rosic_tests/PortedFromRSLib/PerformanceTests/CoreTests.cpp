#include "CoreTests.h"

void testFlagArray(std::string &reportString)
{
  static const rsUint32 numFlags = 1000000;

  ::ProcessorCycleCounter counter;
  double cyclesPerFlag;
  int i;

  counter.init();
  rsFlagArray a(numFlags);
  cyclesPerFlag = (double) counter.getNumCyclesSinceInit() / numFlags;
  printPerformanceTestResult("rsFlagArray, Creation", cyclesPerFlag);

  counter.init();
  rsUint64 dummy = a.getNumTrueFlags();
  cyclesPerFlag = (double) counter.getNumCyclesSinceInit() / numFlags;
  printPerformanceTestResult("rsFlagArray, Count true flags", cyclesPerFlag);

  a.setAllFalse();
  counter.init();
  for(i = 0; i < numFlags; i++)
    a.setFlagTrue(i);
  cyclesPerFlag = (double) counter.getNumCyclesSinceInit() / numFlags;
  printPerformanceTestResult("rsFlagArray, Write", cyclesPerFlag);

  counter.init();
  bool b;
  for(i = 0; i < numFlags; i++)
    b = a.isFlagTrue(i);
  cyclesPerFlag = (double) counter.getNumCyclesSinceInit() / numFlags;
  printPerformanceTestResult("rsFlagArray, Read", cyclesPerFlag);
}

