#include "AutomaticTests.h"
//using namespace rsTestRomos;

namespace rsTestRomos
{

bool runUnitTests()
{
  UnitTestRunner testRunner;
  bool testsPassed = testRunner.runAllTestsAndPrintResultsToConsole();
  return testsPassed;
}

bool testSorting(bool verboseOutput)
{
  int xMin = 0;
  int xMax = 50;
  int yMin = 0;
  int yMax = 50;
  unsigned int numChildrenToCreate = 50;

  //ContainerModule *testModule = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER);
  ContainerModule* testModule = (ContainerModule*)moduleFactory.createModule("Container");

  unsigned int i;
  int x, y;
  RAPT::rsRandomUniform((double)xMin, (double)xMax, 1);
  for(i = 0; i < numChildrenToCreate; i++)
  {
    //romos::Module *childModule = ModuleFactory::createModule(ModuleTypeRegistry::UNIT_DELAY);
    romos::Module* childModule = moduleFactory.createModule("UnitDelay");
    x = (int)RAPT::rsRandomUniform((double)xMin, (double)xMax);
    y = (int)RAPT::rsRandomUniform((double)yMin, (double)yMax);
    childModule->setModuleName(std::string("Delay ") + std::to_string((int)i)); // to keep track of insertion order
    childModule->setPositionXY(x, y);
    testModule->addChildModule(childModule, false);
  }

  std::vector<romos::Module*> childModulesVectorUnsorted = testModule->getChildModules();
  testModule->sortChildModuleArray();
  std::vector<romos::Module*> childModulesVectorSorted = testModule->getChildModules();

  //rosic::Array<romos::Module*> childModulesArray  = testModule->getChildModules();
  //std::vector<romos::Module*>  childModulesVectorUnsorted; 
  //childModulesArray.toVectorSTL(childModulesVectorUnsorted);

  //testModule->sortChildModuleArray();
  //childModulesArray  = testModule->getChildModules();

  //std::vector<romos::Module*>  childModulesVectorSorted; 
  //childModulesArray.toVectorSTL(childModulesVectorSorted);

  // output to the console (comment out, if not desired):
  if(verboseOutput == true)
  {
    printf("%s", "Unsorted:\n");
    for(i=0; i<childModulesVectorUnsorted.size(); i++)
    {
      x = childModulesVectorUnsorted[i]->getPositionX();
      y = childModulesVectorUnsorted[i]->getPositionY();
      printf("%s %d %s %d %s", "x = ", x, ", y = ", y, "\n");
    }
    printf("%s", "\n");
    printf("%s", "Sorted:\n");
    for(i=0; i<childModulesVectorSorted.size(); i++)
    {
      x = childModulesVectorSorted[i]->getPositionX();
      y = childModulesVectorSorted[i]->getPositionY();
      printf("%s %d %s %d %s", "x = ", x, ", y = ", y, "\n");
    }
  }

  // automatic correctness check:
  bool result = true;
  int x1, y1;
  for(i=0; i<childModulesVectorSorted.size()-1; i++)
  {
    x  = childModulesVectorSorted[i]->getPositionX();
    y  = childModulesVectorSorted[i]->getPositionY();
    x1 = childModulesVectorSorted[i+1]->getPositionX();
    y1 = childModulesVectorSorted[i+1]->getPositionY();
    if(x > x1)
      result = false;
    if(x == x1 && y > y1)
      result = false;
  }

  // cleanup and return:
  moduleFactory.deleteModule(testModule);
  if(result == false)
    printf("%s", "!!! Sorting failed !!!\n");
  else
    printf("%s", "Sorting passed\n");
  return result;
}

bool testModuleTypeRegistry()
{
  //ModuleTypeRegistry* typeRegistry = romos::ModuleTypeRegistry::getSoleInstance();

  //double blah = x[0][0];

  return false;
}

bool testGain(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createGain("Gain", 0, 0, false);
  RAPT::rsArrayTools::multiply(x[0][0], x[0][1], d[0][0], N);  // create desired output
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);

  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "Gain", 0.0);
}

bool testSumDiff(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createSumDiff("SumDiff", 0, 0, false);
  RAPT::rsArrayTools::add(x[0][0], x[0][1], d[0][0], N);
  RAPT::rsArrayTools::subtract(x[0][0], x[0][1], d[0][1], N);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 2, N, "SumDiff", 0.0);
}

bool testWrappedSumDiff(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createSumDiff("WrappedSumDiff", 0, 0, false);
  RAPT::rsArrayTools::add(x[0][0], x[0][1], d[0][0], N);
  RAPT::rsArrayTools::subtract(x[0][0], x[0][1], d[0][1], N);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 2, N, "WrappedSumDiff", 0.0);
}

bool testSummedDiffs(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createSummedDiffs("SummedDiffs", 0, 0, false);
  getDesiredOutputForSummedDiffs(N, px0, pd0);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 4, N, "SummedDiffs", 0.0);
}

bool testMovingAverage(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createMovingAverage("MovingAverage", 0, 0, false);
  getDesiredOutputForMovingAverage(N, x[0][0], x[0][1], x[0][2], d[0][0]);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "MovingAverage", 0.0);
}

bool testLeakyIntegrator(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createLeakyIntegrator("LeakyIntegrator", 0, 0, false);
  getDesiredOutputForLeakyIntegrator(N, x[0][0], x[0][1], d[0][0]);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "LeakyIntegrator", 0.0);
}

bool testLeakyIntegratorDoubleDelay(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createLeakyIntegrator("LeakyIntegratorDoubleDelay", 0, 0, false);
  getDesiredOutputForLeakyIntegratorDoubleDelay(N, x[0][0], x[0][1], d[0][0]);
  //romos::Module *identity = ((ContainerModule*) testModule)->getChildModulesWithTypeOld(ModuleTypeRegistry::IDENTITY).at(0);
  romos::Module* identity = ((ContainerModule*)testModule)->getChildModulesWithType("Identity").at(0);
  identity->setPositionXY(17, 2);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "LeakyIntegratorDoubleDelay", 0.0);
}

bool testTestFilter1(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createTestFilter1("TestFilter1", 0, 0, false);
  getDesiredOutputForTestFilter1(N, x[0][0], x[0][1], x[0][2], x[0][3], d[0][0], d[0][1], d[0][2]);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[2][0], y[2][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 3, N, "TestFilter1", 0.0);
}

bool testBiquadMacro(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createBiquadMacro("BiquadMacro", 0, 0, false);
  getDesiredOutputForBiquad(N, x[0][0], x[0][1], x[0][2], x[0][3], x[0][4], x[0][5], d[0][0]);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  //Plotter::plotData(N, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "BiquadMacro", 0.0);
}

bool testBiquadAtomic(bool verboseOutput)
{
  //romos::Module *testModule = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD);
  romos::Module* testModule = moduleFactory.createModule("Biquad");
  getDesiredOutputForBiquad(N, x[0][0], x[0][1], x[0][2], x[0][3], x[0][4], x[0][5], d[0][0]);
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, N, "BiquadAtomic", 1.e-13); // use test with tolerance
}

bool testContainerizationAddedConstants(bool verboseOutput)
{
  static const int numIterations = 50;

  romos::Module* testModule = TestModuleBuilder::createAddedConstants("AddedConstants", 0, 0, false);

  //testModule->resetState();  
  testModule->resetStateForAllVoices();
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  RAPT::rsArrayTools::copy(y[0][0], d[0][0], N);

  for(int i=0; i<numIterations; i++)
  {
    randomizeContainment(testModule);
    //testModule->resetState();  // resetting the state lets the test fail
    processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
    if(verboseOutput == true)
      printModuleStructure(testModule, 0);
    if(!RAPT::rsArrayTools::equal(y[0][0], d[0][0], N))
    {
      printf("%s", "!!! ContainerizationAddedConstants failed !!!\n");
      //printModuleStructure(testModule, 0);
      //Plotter::plotData(N, t, d[0][0], y[0][0]);
      moduleFactory.deleteModule(testModule);
      return false;
    }
  }

  moduleFactory.deleteModule(testModule);
  printf("%s", "ContainerizationAddedConstants passed\n");
  return true;
}

bool testPinSorting(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createPinSortTest("PinSorting", 0, 0, false);

  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  processModuleInFrames(testModule, N, ppx, ppd, NULL, false);

  // retrieve pointers to some embedded modules (assumes certain order):
  romos::ContainerModule* inner = dynamic_cast<romos::ContainerModule*> (((ContainerModule*)testModule)->getChildModule(3));
  rassert(inner != NULL);
  romos::Module* in1   = inner->getAudioInputModule(0);
  romos::Module* in2   = inner->getAudioInputModule(1);
  romos::Module* in3   = inner->getAudioInputModule(2);
  romos::Module* out1  = inner->getAudioOutputModule(0);
  romos::Module* out2  = inner->getAudioOutputModule(1);
  romos::Module* out3  = inner->getAudioOutputModule(2);

  bool result = true;

  //Plotter::plotData(N, t, d[0][0], y[0][0]);

  // test all possible permutations of the input order:
  exchangeModulePositions(in1, in2);      // input order: 1,2,3 -> 2,1,3
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(in1, in3);      // input order: 2,1,3 -> 2,3,1
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(in1, in2);      // input order: 2,3,1 -> 1,3,2
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(in1, in3);      // input order: 1,3,2 -> 3,1,2
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(in1, in2);      // input order: 3,1,2 -> 3,2,1
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(in1, in3);      // input order: 3,2,1 -> 1,2,3 (again)

  // test all possible permutations of the output order:
  exchangeModulePositions(out1, out2);      // output order: 1,2,3 -> 2,1,3
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(out1, out3);      // output order: 2,1,3 -> 2,3,1
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(out1, out2);      // output order: 2,3,1 -> 1,3,2
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(out1, out3);      // output order: 1,3,2 -> 3,1,2
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(out1, out2);      // output order: 3,1,2 -> 3,2,1
  processModuleInFrames(testModule, N, ppx, ppy, NULL, false);
  result &= checkResult(py0, pd0, 3, N, 0.0);
  exchangeModulePositions(out1, out3);      // output order: 3,2,1 -> 1,2,3 (again)

  //Plotter::plotData(N, t, d[0][0], y[0][0]);

  moduleFactory.deleteModule(testModule);
  if(result == false)
    printf("%s", "!!! PinSorting failed !!!\n");
  else
    printf("%s", "PinSorting passed\n");
  return result;
}

bool testAdderBlock(bool verboseOutput)
{
  //romos::Module *testModule = ModuleFactory::createModule(ModuleTypeRegistry::ADDER);
  romos::Module* testModule = moduleFactory.createModule("Adder");
  RAPT::rsArrayTools::add(x[0][0], x[0][1], d[0][0], maxNumFrames);  // establish desired result
  bool result = checkBlockProcessingAndPrintResult(testModule, ppx, ppy, ppd, maxNumFrames, 50, "AdderBlock", 0.0);
  moduleFactory.deleteModule(testModule);
  return result;
}

bool testBiquadAtomicBlock(bool verboseOutput)
{
  //romos::Module *testModule = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD);
  romos::Module* testModule = moduleFactory.createModule("Biquad");
  getDesiredOutputForBiquad(maxNumFrames, x[0][0], x[0][1], x[0][2], x[0][3], x[0][4], x[0][5], d[0][0]);
  bool result = checkBlockProcessingAndPrintResult(testModule, ppx, ppy, ppd, maxNumFrames, 50, "BiquadAtomicBlock", 1.e-13);
  //Plotter::plotData(200, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return result;
}

bool testBlip(bool verboseOutput)
{
  romos::Module* testModule = TestModuleBuilder::createBlip("Blip", 0, 0, false);

  std::vector<NoteEvent> events = generateNoteOnOffPair(81, 64, 0, maxNumFrames-1);
  processModuleInFrames(testModule, maxNumFrames, ppx, ppy, &events, false);

  getDesiredOutputForFilterBlip(maxNumFrames, 880.0, 20.0, d[0][0]);
  //Plotter::plotData(maxNumFrames, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return checkAndPrintResult(py0, pd0, 1, maxNumFrames, "Blip", 0.0);
}


/*
bool testBlipOneNote(bool verboseOutput)
{
  getDesiredOutputForFilterBlip(maxNumFrames, 880.0, 20.0, d[0][0]);
  fillWithZeros(y[0][0], maxNumFrames);  // to make sure that it doesnt' contain the desired result from a previous test

  romos::ModularSynth *modularSynth = new romos::ModularSynth();
  romos::Module *testModule     = TestModuleBuilder::createBlip("Blip", 20, 4, false);
  romos::Module *outputModuleL  = modularSynth->getTopLevelModule()->getAudioOutputModule(0);

  modularSynth->getTopLevelModule()->addChildModule(testModule, true);
  modularSynth->getTopLevelModule()->addAudioConnection(testModule, 0, outputModuleL, 0);

  double inOutL = 0.0;
  double inOutR = 0.0;
  modularSynth->resetAllVoices();
  modularSynth->noteOn(81, 64);
  for(int n=0; n<maxNumFrames; n++)
  {
    modularSynth->getSampleFrameStereo(&inOutL, &inOutR);
    y[0][0][n] = inOutL;
    y[0][1][n] = inOutR;
  }
  modularSynth->resetAllVoices(); // to clean up the states in the global ProcessingState for subsequent tests

  //Plotter::plotData(maxNumFrames, t, d[0][0], y[0][0]);

  delete modularSynth;
  return checkAndPrintResult(py0, pd0, 1, maxNumFrames, "BlipOneNote", 0.0);
}

bool testBlipTwoNotes(bool verboseOutput)
{
  fillWithZeros(y[0][0], maxNumFrames);  // to make sure that it doesnt' contain the desired result from a previous test

  getDesiredOutputForFilterBlip(maxNumFrames,  880.0, 20.0, d[0][0]);
  getDesiredOutputForFilterBlip(maxNumFrames, 1760.0, 20.0, d[0][1]);
  add(d[0][0], d[0][1], d[0][0], maxNumFrames);

  romos::ModularSynth *modularSynth = new romos::ModularSynth();
  romos::Module *testModule     = TestModuleBuilder::createBlip("Blip", 20, 4, false);

  romos::Module *outputModuleL  = modularSynth->getTopLevelModule()->getAudioOutputModule(0);

  modularSynth->getTopLevelModule()->addChildModule(testModule, true);
  modularSynth->getTopLevelModule()->addAudioConnection(testModule, 0, outputModuleL, 0);


  ((ContainerModule *) testModule)->setPolyphonicRecursively(true);

  std::vector<romos::Module*> childModuleVector = modularSynth->getTopLevelModule()->getChildModules();

  //std::vector<romos::Module*> childModuleVector;
  //modularSynth->getTopLevelModule()->getChildModules().toVectorSTL(childModuleVector);

  double inOutL = 0.0;
  double inOutR = 0.0;
  modularSynth->resetAllVoices();
  modularSynth->noteOn(81, 64);
  modularSynth->noteOn(93, 64);
  for(int n=0; n<maxNumFrames; n++)
  {
    modularSynth->getSampleFrameStereo(&inOutL, &inOutR);
    y[0][0][n] = inOutL;
    y[0][1][n] = inOutR;
  }
  modularSynth->noteOff(81);
  modularSynth->noteOff(93);
  modularSynth->resetAllVoices(); // to clean up the states in the global ProcessingState for subsequent tests

  Plotter::plotData(maxNumFrames, t, d[0][0], y[0][0]);

  delete modularSynth;
  return checkAndPrintResult(py0, pd0, 1, maxNumFrames, "BlipTwoNotes", 0.0);
}
*/

bool testAdderProcessingFunctions(int numVoicesToCheck)
{
  //romos::Module *testModule = ModuleFactory::createModule(ModuleTypeRegistry::ADDER);
  romos::Module* testModule = moduleFactory.createModule("Adder");

  for(int v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::add(x[v][0], x[v][1], d[v][0], maxNumFrames);

  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToCheck, 12);
  bool result = checkProcessingFunctionsAndPrintResults(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0, "Adder", &events);

  moduleFactory.deleteModule(testModule);
  return result;
}

bool testUnitDelayProcessingFunctions(int numVoicesToCheck)
{
  //romos::Module *testModule = ModuleFactory::createModule(ModuleTypeRegistry::UNIT_DELAY);
  romos::Module* testModule = moduleFactory.createModule("UnitDelay");
  for(int v = 0; v < numVoicesToCheck; v++)
    getDesiredOutputForUnitDelay(maxNumFrames, x[v][0], d[v][0]);
  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToCheck, 12);
  bool result = checkProcessingFunctionsAndPrintResults(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "UnitDelay", &events);
  moduleFactory.deleteModule(testModule);
  return result;
}

bool testWrappedAdderProcessingFunctions(int numVoicesToCheck)
{
  romos::Module* testModule =  TestModuleBuilder::createWrappedAdder("WrappedAdder", 0, 0, false);
  for(int v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::add(x[v][0], x[v][1], d[v][0], maxNumFrames);
  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToCheck, 12);
  bool result = checkProcessingFunctionsAndPrintResults(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "WrappedAdder", &events);
  //Plotter::plotData(maxNumFrames, t, d[0][0], y[0][0]);
  moduleFactory.deleteModule(testModule);
  return result;
}

bool testMonoToPoly(int numVoicesToCheck)
{
  romos::Module* testModule =  TestModuleBuilder::createMonoToPoly("MonoToPoly", 0, 0, true);

  for(int v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::fillWithValue(d[v][0], maxNumFrames, -1.0);

  bool result = true;


  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToCheck, 12);

  // the container itself is polyphonic:
  result &= checkProcessingInFramesPolyAndPrintResult(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "MonoToPoly, Container Poly", &events);
  result &= checkProcessingInBlocksPolyAndPrintResult(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "MonoToPoly, Container Poly", &events);

  // now make the container monophonic (and adapt the 0th desired output accordingly):
  RAPT::rsArrayTools::fillWithValue(d[0][0], maxNumFrames, (double)-numVoicesToCheck);
  testModule->setPolyphonic(false);
  result &= checkProcessingInFramesMonoAndPrintResult(testModule, maxNumFrames, ppx, ppy, ppd, 0.0,
    "MonoToPoly, Container Mono", &events);
  result &= checkProcessingInBlocksMonoAndPrintResult(testModule, maxNumFrames, ppx, ppy, ppd, 0.0,
    "MonoToPoly, Container Mono", &events);
    // seems that we need to trigger notes - otherwise the polyphonic modules will output silence


  // introduce a delayed connection and test again (uses frame-wise processing)

  moduleFactory.deleteModule(testModule);
  return result;
}

/*
// obsolete because that conversion now happens automatically?
bool testPolyToMono(int numVoicesToCheck)
{
  romos::Module *testModule =  TestModuleBuilder::createPolyToMono("PolyToMono", 0, 0, true);

  for(int v = 0; v < numVoicesToCheck; v++)
    fillWithValue(d[v][0], maxNumFrames, (double) -numVoicesToCheck);

  bool result = true;

  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToCheck, 12);

  // the container itself is polyphonic:
  result &= checkProcessingInFramesPolyAndPrintResult(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "PolyToMono, Container Poly", &events);
  result &= checkProcessingInBlocksPolyAndPrintResult(testModule, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0,
    "PolyToMono, Container Poly", &events);

  // now make the container monophonic (0th desired output stays the same here):
  testModule->setPolyphonic(false);
  result &= checkProcessingInFramesMonoAndPrintResult(testModule, maxNumFrames, ppx, ppy, ppd, 0.0,
    "PolyToMono, Container Mono", &events);
  result &= checkProcessingInBlocksMonoAndPrintResult(testModule, maxNumFrames, ppx, ppy, ppd, 0.0,
    "PolyToMono, Container Mono", &events);

  moduleFactory.deleteModule(testModule);
  return result;
}
*/

bool testGatedNoteFrequency(int numVoicesToCheck)
{
  romos::ContainerModule* testModule =
    (romos::ContainerModule*) TestModuleBuilder::createGatedNoteFrequency("GatedNoteFrequency", 0, 0, true);

  int noteLengthInFrames = 20;
  std::vector<NoteEvent> events = generateSimultaneousNotes(81, 64, 0, 20, numVoicesToCheck, 12);

  bool result = false;

  /*
  getDesiredOutputForGatedNoteFrequencies(maxNumFrames, &events, ppd, false, false);
  result &= checkProcessingInFramesMonoAndPrintResult(testModule, maxNumFrames, ppx, ppy, ppd, 0.0,
    "GatedNoteFrequency, Mono", &events);

  //... more to come..
  */

  RAPT::rsAssert(false, "plotting code needs update");
  //Plotter::plotData(50, t, d[0][0], y[0][0]);

  moduleFactory.deleteModule(testModule);
  return result;
}

bool testTriggerAndKill(int numVoicesToCheck)
{
  romos::ContainerModule* testModule =
    (ContainerModule*)TestModuleBuilder::createTriggerAndKill("TriggerAndKill", 0, 0, true);

  //int frameIndex, voiceIndex, pinIndex; // frameIndex, voiceIndex, pinIndex

  for(int v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::fillWithValue(d[v][0], maxNumFrames, 0.0);

  /*
  testModule->setPolyphonicRecursively(false);
  processingStatus.resetAllVoices();

  std::vector<NoteEvent> events = generateNoteOnOffPair(81, 64, 10, 20);


  checkProcessingInFramesMonoAndPrintResult(module, numVoicesToCheck, maxNumFrames, ppx, ppy, ppd, 0.0, "TriggerAndKill", events);


  processModuleInFramesMono(testModule, maxNumFrames,

  int nextEventFrame = events[0].deltaFrames;

  for(frameIndex = 0; frameIndex < maxNumFrames; frameIndex++)
    processFrameMono(testModule, px0, py0, frameIndex);
    */


  /*
  //int noteOnIndex1 = 10;

  int noteOnIndices[3] = { 10, 20, 30 };



  for(frameIndex = 0; frameIndex < noteOnIndices[0]; frameIndex++)
    processFrameMono(testModule, px0, py0, frameIndex);

  processingStatus.noteOn(81, 64);

  for(frameIndex = noteOnIndices[0]; frameIndex < maxNumFrames; frameIndex++)
    processFrameMono(testModule, px0, py0, frameIndex);
    */


  //Plotter::plotData(maxNumFrames, t, y[0][0]);
  //Plotter::plotData(50, t, y[0][0]);







  bool result = false;  // preliminary


  moduleFactory.deleteModule(testModule);
  return result;
}


/*
bool testBlipBlock(bool verboseOutput)
{
  romos::Module *testModule = romos::createBlip(0, 0, "Blip");
  romos::Environment::getSoleInstance()->noteTriggered(0, 81, 64);
  processModule(testModule, maxNumFrames, px, py);
  getDesiredOutputForFilterBlip(maxNumFrames, 880.0, 20.0, d[0]);
  //Plotter::plotData(maxNumFrames, t, d[0], y[0]);
  bool result = checkBlockProcessingAndPrintResult(testModule, px, py, pd, maxNumFrames, 50, "BlipBlock", 0.0);
  delete testModule;
  return result;
}
*/

/*
bool testContainerizationWithConstant()
{
  ContainerModule *testModule = new ContainerModule(NULL);
  testModule->setModuleName("Test Module");

  ConstantModule *unity = new ConstantModule(NULL);
  unity->setModuleName(2.5);
  unity->setPositionXY(2, 5);
  testModule->addChildModule(unity, false);

  rosic::Array<romos::Module*> selectedModules;
  selectedModules.appendElement(unity);
  testModule->containerizeModules(selectedModules);

  delete testModule;
  return true;
}
*/

} // end of namespace