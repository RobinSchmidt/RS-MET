#include "romos_PerformanceTest.h"
using namespace rsTestRomos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

PerformanceTest::PerformanceTest(const char *testName) 
: UnitTest(testName)
{
  numFramesPerRun = 2000; 
  numRuns         = 20;
}

PerformanceTest::~PerformanceTest()
{
  moduleFactory.deleteModule(moduleToTest);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// running the tests:

bool PerformanceTest::runTest()
{
  return false;
}

rosic::rsString PerformanceTest::runTestsAndGetReport()
{
  rosic::rsString report = name;
  report.padToLength(nameLength, padCharacter);

  // \todo set inputs to random values....

  voiceAllocator.reset();
  voiceAllocator.noteOn(81,  64);
  voiceAllocator.noteOn(93,  64);
  voiceAllocator.noteOn(105, 64);

  setTestPolyphonic(false);
  report += runFrameWiseTestAndGetReport();
  report += runBlockWiseTestAndGetReport();

  setTestPolyphonic(true);
  report += runFrameWiseTestAndGetReport();
  report += runBlockWiseTestAndGetReport();

  int numPlayingVoices = voiceAllocator.getNumPlayingVoices();

  report.printToStandardOutput();
  report += "\n";
  return report;
}

void PerformanceTest::setTestPolyphonic(bool shouldBePolyphonic)
{
  if( moduleToTest->isContainerModule() )
    ((ContainerModule *) moduleToTest)->setPolyphonicRecursively(shouldBePolyphonic);
  else
    moduleToTest->setPolyphonic(shouldBePolyphonic);
}

rosic::rsString PerformanceTest::runFrameWiseTestAndGetReport()
{
  double minCyclesPerFrame = pow(10.0, 100.0);
  for(int runIndex = 1; runIndex <= numRuns; runIndex++)
  {
    counter.init();
    for(int frameIndex = 0; frameIndex < numFramesPerRun; frameIndex++)
      moduleToTest->processSampleFrame();
    double cyclesPerFrame = (double) counter.getNumCyclesSinceInit() / (double) numFramesPerRun;
    if( cyclesPerFrame < minCyclesPerFrame && cyclesPerFrame > 0.0 )
      minCyclesPerFrame = cyclesPerFrame;
  }

  if( moduleToTest->isPolyphonic() )
    minCyclesPerFrame /= voiceAllocator.getNumPlayingVoices();
  char charBuffer[32];
  sprintf(charBuffer, "%.2f", minCyclesPerFrame);
  rosic::rsString cyclesString = rosic::rsString(charBuffer);
  cyclesString.prePadToLength(numberLength, padCharacter);
  return cyclesString;
}

rosic::rsString PerformanceTest::runBlockWiseTestAndGetReport()
{
  int blockSize       = processingStatus.getBufferSize();
  int numBlocksPerRun = numFramesPerRun / blockSize;
  int lastBlockSize   = numFramesPerRun - numBlocksPerRun*blockSize;
  double minCyclesPerFrame = pow(10.0, 100.0);
  for(int runIndex = 1; runIndex <= numRuns; runIndex++)
  {
    counter.init();
    for(int blockIndex = 0; blockIndex < numBlocksPerRun; blockIndex++)
      moduleToTest->processBlockOfSamples(blockSize);
    if( lastBlockSize > 0 )
      moduleToTest->processBlockOfSamples(lastBlockSize);
    double cyclesPerFrame = (double) counter.getNumCyclesSinceInit() / (double) numFramesPerRun;
    if( cyclesPerFrame < minCyclesPerFrame && cyclesPerFrame > 0.0 )
      minCyclesPerFrame = cyclesPerFrame;
  }

  if( moduleToTest->isPolyphonic() )
    minCyclesPerFrame /= voiceAllocator.getNumPlayingVoices();
  char charBuffer[32];
  sprintf(charBuffer, "%.2f", minCyclesPerFrame);
  rosic::rsString cyclesString = rosic::rsString(charBuffer);
  cyclesString.prePadToLength(numberLength, padCharacter);  
  return cyclesString;
}
