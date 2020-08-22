#include "romos_ProcessingTest.h"
using namespace rsTestRomos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

ProcessingTest::ProcessingTest(const char *testName)
: UnitTest(testName)
{
  //tolerance          = 0.0;       // works for msc but not gcc
  tolerance          = 1.e-13;

  moduleToTest       = NULL;      // subclasses are responsible for assigning this pointer
  numFramesToProcess = 647;
  numVoicesToUse     = 3;         // should be <= number of available voices, otherwise the test fails

  fillInputSignalArraysRandomly(1);
  clearOutputSignalArrays();
  clearDesiredOutputSignalArrays();
  fillTimeAxisWithSampleIndices();
  voiceAllocator.reset();  // maybe use a function initialize which is even stronger (resets even more state variables)

  events = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToUse, 12);


  // create modules that we use to feed the test module with input signals and to retrieve the produced output:
  int i;
  for(i = 0; i < maxNumInputs; i++)
  {
    //inputFeederModules[i] = ModuleFactory::createModule(
    //  ModuleTypeRegistry::IDENTITY, rosic::rsString("In") + rosic::rsString(i+1), 0, i, true);

    inputFeederModules[i] = moduleFactory.createModule(
      "Identity", std::string("In") + std::to_string(i+1), 0, i, true);
    // shouldn't we use a proper AudioInput module isntead?
  }
  for(i = 0; i < maxNumOutputs; i++)
  {
    //outputRetrieverModules[i] = ModuleFactory::createModule(ModuleTypeRegistry::IDENTITY,
    //  rosic::rsString("Out") + rosic::rsString(i+1), 0, i, true);

    outputRetrieverModules[i] = moduleFactory.createModule(
      "Identity", std::string("Out") + std::to_string(i+1), 0, i, true);
  }
}

ProcessingTest::~ProcessingTest()
{
  moduleFactory.deleteModule(moduleToTest);
  for(int i = 0; i < maxNumInputs; i++)
    moduleFactory.deleteModule(inputFeederModules[i]);
  for(int i = 0; i < maxNumOutputs; i++)
    moduleFactory.deleteModule(outputRetrieverModules[i]);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// running the tests:

bool ProcessingTest::runTest()
{
  initTest();

  bool monoInFramesPassed = false;
  bool monoInBlocksPassed = false;
  bool polyInFramesPassed = false;
  bool polyInBlocksPassed = false;


  // todo - include a self-test to make sure that the test can actually fail by filling the desiredOutputs only after

  setTestPolyphonic(false);

  clearOutputSignalArrays();
  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();
  monoInFramesPassed = doOutputsMatchDesiredOutputs(false);
  handleTestResult(monoInFramesPassed);

  clearOutputSignalArrays();
  moduleToTest->resetStateForAllVoices();
  processModuleInBlocks();
  monoInBlocksPassed = doOutputsMatchDesiredOutputs(false);
  handleTestResult(monoInBlocksPassed);

  setTestPolyphonic(true);

  clearOutputSignalArrays();
  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();
  polyInFramesPassed = doOutputsMatchDesiredOutputs(true);
  handleTestResult(polyInFramesPassed);

  clearOutputSignalArrays();
  moduleToTest->resetStateForAllVoices();
  processModuleInBlocks();
  polyInBlocksPassed = doOutputsMatchDesiredOutputs(true);
  handleTestResult(polyInBlocksPassed);

  return monoInFramesPassed && monoInBlocksPassed && polyInFramesPassed && polyInBlocksPassed;
}

void ProcessingTest::handleTestResult(bool didTestPass)
{
  if(!didTestPass) {
    //RS_DEBUG_BREAK;
    //plotDesiredAndActualOutput(0, 0, numFramesToProcess, 0);
  }
}

void ProcessingTest::initTest()
{
  connectTestModuleToInputFeederModules();
  connectTestModuleToOutputRetrieverModules();
  fillInputSignalArraysWithTestSignal();
  fillDesiredOutputSignalArrays(false);
}

void ProcessingTest::connectTestModuleToInputFeederModules()
{
  rassert( moduleToTest != NULL ); // this function should be called after the moduleToTest has been allocated
  for(unsigned int i = 0; i < moduleToTest->getNumInputPins(); i++)
  {
    moduleToTest->connectInputPinTo(i, inputFeederModules[i], 0);
    int dummy = 0;
  }
}

void ProcessingTest::connectTestModuleToOutputRetrieverModules()
{
  for(unsigned int i = 0; i < moduleToTest->getNumOutputPins(); i++)
  {
    outputRetrieverModules[i]->connectInputPinTo(0, moduleToTest, i);
    int dummy = 0;
  }
}


bool ProcessingTest::doOutputsMatchDesiredOutputs(bool polyphonic)
{
  int highestVoiceIndexToCheck;
  if( polyphonic )
    highestVoiceIndexToCheck = maxNumVoices - 1;
  else
    highestVoiceIndexToCheck = 0;

  bool outputsMatch = true;
  for(int voiceIndex = 0; voiceIndex <= highestVoiceIndexToCheck; voiceIndex++)
  {
    for(unsigned int pinIndex = 0; pinIndex < moduleToTest->getNumOutputPins(); pinIndex++)
      outputsMatch &= RAPT::rsArrayTools::almostEqual(outputs[voiceIndex][pinIndex], desiredOutputs[voiceIndex][pinIndex],
      numFramesToProcess, tolerance);
  }

  return outputsMatch;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// filling of the internal arrays:

void ProcessingTest::fillInputSignalArraysRandomly(int seed)
{
  RAPT::rsRandomUniform(-1.0, 1.0, seed);
  for(int voiceIndex = 0; voiceIndex < maxNumVoices; voiceIndex++)
  {
    for(int pinIndex = 0; pinIndex < maxNumInputs; pinIndex++)
    {
      for(int frameIndex = 0; frameIndex < maxNumFrames; frameIndex++)
      {
        inputs[voiceIndex][pinIndex][frameIndex] = RAPT::rsRandomUniform(-1.0, 1.0);
      }
    }
  }
}

void ProcessingTest::clearInputSignalArrays()
{
  memset(inputs, 0, maxNumVoices * maxNumInputs * maxNumFrames * sizeof(double));
}

void ProcessingTest::clearOutputSignalArrays()
{
  memset(outputs, 0, maxNumVoices * maxNumOutputs * maxNumFrames * sizeof(double));
}

void ProcessingTest::clearDesiredOutputSignalArrays()
{
  memset(desiredOutputs, 0, maxNumVoices * maxNumOutputs * maxNumFrames * sizeof(double));
}

void ProcessingTest::fillTimeAxisWithSampleIndices()
{
  for(int frameIndex = 0; frameIndex < maxNumFrames; frameIndex++)
    timeAxis[frameIndex] = (double) frameIndex;
}

void ProcessingTest::fillInputSignalArraysWithTestSignal()
{
  fillInputSignalArraysRandomly(1);
}

void ProcessingTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  clearDesiredOutputSignalArrays();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:




//-----------------------------------------------------------------------------------------------------------------------------------------
// processing:

// \todo maybe factor out the common structure of processModuleInFrames and processModuleInBlocks

void ProcessingTest::processModuleInFrames()
{
  voiceAllocator.reset();
  if( events.empty() )
    processModuleInFramesNoEvents(numFramesToProcess, 0);
  else
  {
    int frameIndex       = 0;
    int numEventsHandled = 0;
    while( frameIndex < numFramesToProcess )
    {
      romos::NoteEvent e = events.at(numEventsHandled);
      int numFramesUntilNextEvent = e.getDeltaFrames() - frameIndex;
      processModuleInFramesNoEvents(numFramesUntilNextEvent, frameIndex);        // process chunk until the next event
        // there, process 1st frame, reset the trigger flags and then process other frames


      voiceAllocator.noteOn(e.getKey(), e.getVelocity());
      numEventsHandled++;
      frameIndex += numFramesUntilNextEvent;
      if( numEventsHandled == events.size() )                                    // process tail after all events have been handled
      {
        processModuleInFramesNoEvents(numFramesToProcess - frameIndex, frameIndex);
        frameIndex = numFramesToProcess;
      }

      //voiceAllocator.resetTriggerFlags();  // perhaps we need this also for the block processing function - but then we must process
      //                                     // the 1st sample of the block separately

    }
  }
}

void ProcessingTest::processModuleInBlocks()
{
  voiceAllocator.reset();
  if( events.empty() )
    processModuleInBlocksNoEvents(numFramesToProcess, 0);
  else
  {
    int maxBlockSize     = romos::processingStatus.getBufferSize();
    int blockStart       = 0;
    int numEventsHandled = 0;
    while( blockStart < numFramesToProcess )
    {
      romos::NoteEvent e = events.at(numEventsHandled);
      int numFramesUntilNextEvent = e.getDeltaFrames() - blockStart;
      processModuleInBlocksNoEvents(numFramesUntilNextEvent, blockStart);
      voiceAllocator.noteOn(e.getKey(), e.getVelocity());
      numEventsHandled++;
      blockStart += numFramesUntilNextEvent;
      if( numEventsHandled == events.size() )
      {
        processModuleInBlocksNoEvents(numFramesToProcess - blockStart, blockStart);
        blockStart = numFramesToProcess;
      }
    }
  }
}


void ProcessingTest::processModuleInFramesNoEvents(int numFrames, int startIndex)
{
  int endIndex = RAPT::rsMin(startIndex + numFrames - 1, numFramesToProcess - 1);
  for(int frameIndex = startIndex; frameIndex <= endIndex; frameIndex++)
  {
    processFrame(frameIndex);
    voiceAllocator.resetTriggerFlags();
  }

  //for(int frameIndex = startIndex; frameIndex < startIndex + numFrames; frameIndex++)
  //  processFrame(frameIndex);
}

void ProcessingTest::processModuleInBlocksNoEvents(int numFrames, int startIndex)
{
  int maxBlockSize = romos::processingStatus.getBufferSize();
  int blockStart   = startIndex;
  int endIndex     = RAPT::rsMin(startIndex + numFrames - 1, numFramesToProcess - 1);
  while( blockStart <= endIndex )
  {
    int blockSize = (int) ::round(RAPT::rsRandomUniform(1.0, maxBlockSize));
    if( blockStart + blockSize > endIndex + 1 )
      blockSize = endIndex - blockStart + 1;
    processBlock(blockStart, blockSize);
    blockStart += blockSize;
  }
}

/*
// old:
void ProcessingTest::processModuleInBlocksNoEvents(int numFrames, int startIndex)
{
  int maxBlockSize = romos::processingStatus.getAllocatedBlockSize();
  int blockStart   = startIndex;
  while( blockStart < startIndex + numFrames )
  {
    int blockSize = (int) round(randomUniform(1.0, maxBlockSize));
    if( blockStart + blockSize > startIndex + numFrames )
      blockSize = startIndex + numFrames - blockStart;
    processBlock(blockStart, blockSize);
    blockStart += blockSize;
  }
}
*/

void ProcessingTest::processFrame(int frameIndex)
{
  establishInputBlock(frameIndex, 1);
  moduleToTest->processSampleFrame();
  retrieveOutputBlock(frameIndex, 1);
}

void ProcessingTest::processBlock(int blockStart, int blockSize)
{
  establishInputBlock(blockStart, blockSize);
  moduleToTest->processBlockOfSamples(blockSize);
  retrieveOutputBlock(blockStart, blockSize);
  //Plotter::plotData(numFramesToProcess, timeAxis, desiredOutputs[0][0], outputs[0][0]);
}




// \todo maybe factor out the common structure of establishInputBlock/retrieveOutputBlock - the differences are the direction of copying
// and the pointer and array that is used - maybe this difference can be wrapped into a function
void ProcessingTest::establishInputBlock(int blockStart, int blockSize)
{
  if( !moduleToTest->isPolyphonic() )
  {
    for(unsigned int pinIndex = 0; pinIndex < maxNumInputs; pinIndex++)
    {
      double *outputPointer = inputFeederModules[pinIndex]->getOutputPointer(0);
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
      {
        int offset = inputFeederModules[pinIndex]->getOutputPinMemoryOffset(frameIndex, 0, 0);
        outputPointer[offset] = inputs[0][pinIndex][blockStart+frameIndex];
        double dbg = outputPointer[offset];
        int dummy = 0;
      }
    }
  }
  else
  {
    int numPlayingVoices           = romos::processingStatus.getNumPlayingVoices();
    const int *playingVoiceIndices = processingStatus.getPlayingVoiceIndices();
    for(int currentVoice = 0; currentVoice < numPlayingVoices; currentVoice++)
    {
      int voiceIndex = playingVoiceIndices[currentVoice];
      for(unsigned int pinIndex = 0; pinIndex < maxNumOutputs; pinIndex++)
      {
        double *outputPointer = inputFeederModules[pinIndex]->getOutputPointer(0);
        for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
        {
          int offset = inputFeederModules[pinIndex]->getOutputPinMemoryOffset(frameIndex, voiceIndex, 0);
          outputPointer[offset] = inputs[voiceIndex][pinIndex][blockStart+frameIndex];
          double dbg = outputPointer[offset];
          int dummy = 0;
        }
      }
    }
  }
}

void ProcessingTest::retrieveOutputBlock(int blockStart, int blockSize)
{
  for(int i = 0; i < maxNumOutputs; i++)
    outputRetrieverModules[i]->processBlockOfSamples(blockSize); // copies signals into the retriver-modules' output buffers

  //double *outputPointer = moduleToTest->getOutputPointer(0);

  //double dbg1[12];
  //memcpy(dbg1, outputPointer, 12*sizeof(double));


  if( !moduleToTest->isPolyphonic() )
  {
    for(unsigned int pinIndex = 0; pinIndex < moduleToTest->getNumOutputPins(); pinIndex++)
    {
      double *outputPointer = outputRetrieverModules[pinIndex]->getOutputPointer(0);
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
      {
        int offset = outputRetrieverModules[pinIndex]->getOutputPinMemoryOffset(frameIndex, 0, 0);
        outputs[0][pinIndex][blockStart+frameIndex] = outputPointer[offset];
        double dbg = outputs[0][pinIndex][blockStart+frameIndex];
        int dummy = 0;
      }
    }
  }
  else
  {
    int numPlayingVoices           = romos::processingStatus.getNumPlayingVoices();
    const int *playingVoiceIndices = processingStatus.getPlayingVoiceIndices();
    for(int currentVoice = 0; currentVoice < numPlayingVoices; currentVoice++)
    {
      int voiceIndex = playingVoiceIndices[currentVoice];
      for(unsigned int pinIndex = 0; pinIndex < moduleToTest->getNumOutputPins(); pinIndex++)
      {
        double *outputPointer = outputRetrieverModules[pinIndex]->getOutputPointer(0);
        for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)
        {
          int offset = outputRetrieverModules[pinIndex]->getOutputPinMemoryOffset(frameIndex, voiceIndex, 0);
          outputs[voiceIndex][pinIndex][blockStart+frameIndex] = outputPointer[offset];
          double dbg = outputs[voiceIndex][pinIndex][blockStart+frameIndex];
          int dummy = 0;
        }
      }
    }
  }
  //DEBUG_BREAK;  // check updated code above

  //plotDesiredAndActualOutput(0, 0, numFramesToProcess, 0);
  //int dummy = 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void ProcessingTest::setTestModulePolyphonyRecursively(bool shouldBePolyphonic)
{
  if( moduleToTest->isContainerModule() )
    ((ContainerModule*) moduleToTest)->setPolyphonicRecursively(shouldBePolyphonic);
  else
    moduleToTest->setPolyphonic(shouldBePolyphonic);
}

void ProcessingTest::setTestPolyphonic(bool shouldBePolyphonic)
{
  setTestModulePolyphonyRecursively(shouldBePolyphonic);
  clearDesiredOutputSignalArrays();
  fillDesiredOutputSignalArrays(shouldBePolyphonic);
  for(int i = 0; i < maxNumOutputs; i++)
    outputRetrieverModules[i]->updateInputPointersAndInFrameStrides();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// information output functions:


void ProcessingTest::plotDesiredAndActualOutput(int voiceIndex, int pinIndex, int numFramesToPlot,
  int startFrame)
{
#ifdef RS_DEBUG_PLOTTING
  GNUPlotter plt;
  plt.plot(numFramesToPlot,
    &timeAxis[startFrame],
    desiredOutputs[voiceIndex][pinIndex],
    outputs[voiceIndex][pinIndex]);
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// internal functions:

