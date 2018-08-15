#include "romos_ConcreteProcessingTests.h"
using namespace romos;

/*
PolyphonicProcessingTest::PolyphonicProcessingTest(const char *testName)
: ProcessingTest(testName)
{
  numVoicesToUse = 3; // should be <= number of available voices, otherwise the test fails
  events = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToUse, 12);
}
*/



IdentityTest::IdentityTest()
: ProcessingTest("Identity")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::IDENTITY);
}
void IdentityTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    copyBuffer(inputs[v][0], desiredOutputs[v][0], numFramesToProcess);
}


AdderTest::AdderTest()
: ProcessingTest("Adder")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::ADDER);
}
void AdderTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    add(inputs[v][0], inputs[v][1], desiredOutputs[v][0], numFramesToProcess);
}


Adder3Test::Adder3Test()
: ProcessingTest("Adder3")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::ADDER_3);
}
void Adder3Test::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    for(int n = 0; n < numFramesToProcess; n++)
      desiredOutputs[v][0][n] = inputs[v][0][n] + inputs[v][1][n] + inputs[v][2][n];
  }
}


Adder4Test::Adder4Test()
: ProcessingTest("Adder4")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::ADDER_4);
}
void Adder4Test::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    for(int n = 0; n < numFramesToProcess; n++)
      desiredOutputs[v][0][n] = inputs[v][0][n] + inputs[v][1][n] + inputs[v][2][n] + inputs[v][3][n]; 
  }
}


Adder5Test::Adder5Test()
: ProcessingTest("Adder5")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::ADDER_5);
}
void Adder5Test::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    for(int n = 0; n < numFramesToProcess; n++)
      desiredOutputs[v][0][n] = inputs[v][0][n] + inputs[v][1][n] + inputs[v][2][n] + inputs[v][3][n] + inputs[v][4][n]; 
  }
}



WrappedAdderTest::WrappedAdderTest()
: ProcessingTest("WrappedAdder")
{
  moduleToTest = TestModuleBuilder::createWrappedAdder("WrappedAdder", 0, 0, true);
}
void WrappedAdderTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    add(inputs[v][0], inputs[v][1], desiredOutputs[v][0], numFramesToProcess);
}


SubtractorTest::SubtractorTest()
: ProcessingTest("Subtractor")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::SUBTRACTOR);
}
void SubtractorTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    subtract(inputs[v][0], inputs[v][1], desiredOutputs[v][0], numFramesToProcess);
}


UnitDelayTest::UnitDelayTest()
: ProcessingTest("UnitDelay")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::UNIT_DELAY);
}
void UnitDelayTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forUnitDelay(numFramesToProcess, inputs[v][0], desiredOutputs[v][0]);
}



NoiseGeneratorTest::NoiseGeneratorTest()
: ProcessingTest("NoiseGeneratorTest")
{
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::WHITE_NOISE);
}
void NoiseGeneratorTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forWhiteNoiseUniform(numFramesToProcess, desiredOutputs[v][0], 0);
}



SumDiffProdTest::SumDiffProdTest()
: ProcessingTest("SumDiffProd")
{
  moduleToTest = TestModuleBuilder::createSumDiffProd("SumDiffProd", 0, 0, true);
}
void SumDiffProdTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    add(     inputs[v][0], inputs[v][1], desiredOutputs[v][0], numFramesToProcess);
    subtract(inputs[v][0], inputs[v][1], desiredOutputs[v][1], numFramesToProcess);
    multiply(inputs[v][0], inputs[v][1], desiredOutputs[v][2], numFramesToProcess);
  }
}


WrappedSumDiffProdTest::WrappedSumDiffProdTest()
: ProcessingTest("WrappedSumDiffProd")
{
  moduleToTest = TestModuleBuilder::createWrappedSumDiffProd("WrappedSumDiffProd", 0, 0, true);
}
void WrappedSumDiffProdTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    add(     inputs[v][0], inputs[v][1], desiredOutputs[v][0], numFramesToProcess);
    subtract(inputs[v][0], inputs[v][1], desiredOutputs[v][1], numFramesToProcess);
    multiply(inputs[v][0], inputs[v][1], desiredOutputs[v][2], numFramesToProcess);
  }
}


WrappedAdderNTest::WrappedAdderNTest()
: ProcessingTest("WrappedAdderN")
{
  moduleToTest = TestModuleBuilder::createWrappedAdderN("WrappedAdderN", 0, 0, true);
}
void WrappedAdderNTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    copyBuffer(inputs[v][0],  desiredOutputs[v][0], numFramesToProcess);
    scale(desiredOutputs[v][0], desiredOutputs[v][0], numFramesToProcess, (double) getAdderNumConnectedInputPins());
  }
}  
bool WrappedAdderNTest::runTest()
{
  bool result = true;

  result &= getAdderNumInputPins() == 11;
  result &= ProcessingTest::runTest();

  // remove a connection from the middle - should not remove the input pin:
  removeConnection(5);  
  result &= getAdderNumInputPins() == 11;
  result &= ProcessingTest::runTest();

  // remove last connection (index 9) - should remove the last input pin (which is currently at index 9+1 = 10) and reduce the number of 
  // pins to input pins to 10:
  removeConnection(9);  
  result &= getAdderNumInputPins() == 10;
  result &= ProcessingTest::runTest();

  int dummy = 0;
  return result;
}   
void WrappedAdderNTest::removeConnection(int index)
{
  romos::Module* audioInput = ((ModuleContainer*)moduleToTest)->getChildModule(0);
  romos::Module* adderN     = ((ModuleContainer*)moduleToTest)->getChildModule(1);
  ((ModuleContainer*)moduleToTest)->deleteAudioConnection(audioInput, 0, adderN, index);
}
int WrappedAdderNTest::getAdderNumInputPins()
{
  romos::Module* adderN = ((ModuleContainer*)moduleToTest)->getChildModule(1);
  return adderN->getNumInputPins();
}
int WrappedAdderNTest::getAdderNumConnectedInputPins()
{
  romos::Module* adderN = ((ModuleContainer*)moduleToTest)->getChildModule(1);
  int result = 0;
  for(unsigned int i = 0; i < adderN->getNumInputPins(); i++)
  {
    if( adderN->isInputPinConnected(i) )
      result++;
  }
  return result;
}

SummedDiffsTest::SummedDiffsTest()
: ProcessingTest("SummedDiffs")
{
  moduleToTest = TestModuleBuilder::createSummedDiffs(name, 0, 0, false);
}
void SummedDiffsTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  double d1, d2, d3, d4;
  for(int v = 0; v < numVoicesToUse; v++) 
  {
    for(int n = 0; n < numFramesToProcess; n++)
    {
      d1 = inputs[v][0][n] - inputs[v][1][n];
      d2 = inputs[v][1][n] - inputs[v][0][n];
      d3 = inputs[v][2][n] - inputs[v][1][n];
      d4 = inputs[v][1][n] - inputs[v][2][n];
      desiredOutputs[v][0][n] = d1 + d2 + d4;
      desiredOutputs[v][1][n] = d2 + d4;
      desiredOutputs[v][2][n] = d1 + d3;
      desiredOutputs[v][3][n] = d1 + d3 + d4;
    }
  }
}


MovingAverageTest::MovingAverageTest()
: ProcessingTest("MovingAverage")
{
  moduleToTest = TestModuleBuilder::createMovingAverage(name, 0, 0, false);
}
void MovingAverageTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    GenerateDesiredOutput::forMovingAverage(numFramesToProcess, inputs[v][0], inputs[v][1], inputs[v][2], desiredOutputs[v][0]);
}



DelayedConnectionTest::DelayedConnectionTest()
: ProcessingTest("DelayedConnection")
{
  moduleToTest = TestModuleBuilder::createDelayedConnection("DelayedConnection", 0, 0, false);
}
void DelayedConnectionTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    GenerateDesiredOutput::forUnitDelay(numFramesToProcess, inputs[v][0], desiredOutputs[v][0]);
}


LeakyIntegratorTest::LeakyIntegratorTest()
: ProcessingTest("LeakyIntegrator")
{
  moduleToTest = TestModuleBuilder::createLeakyIntegrator("LeakyIntegrator", 0, 0, false);
}
void LeakyIntegratorTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    GenerateDesiredOutput::forLeakyIntegrator(numFramesToProcess, inputs[v][0], inputs[v][1], desiredOutputs[v][0]);
}


LeakyIntegratorDoubleDelayTest::LeakyIntegratorDoubleDelayTest()
: ProcessingTest("LeakyIntegratorDoubleDelay")
{
  moduleToTest = TestModuleBuilder::createLeakyIntegrator(name, 0, 0, false);
  romos::Module *identity = ((ModuleContainer*) moduleToTest)->getChildModulesWithTypeOld(ModuleTypeRegistry::IDENTITY).at(0);
  identity->setPositionXY(17, 2);
}
void LeakyIntegratorDoubleDelayTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forLeakyIntegratorDoubleDelay(numFramesToProcess, inputs[v][0], inputs[v][1], desiredOutputs[v][0]);
}


TestFilter1Test::TestFilter1Test()
: ProcessingTest("TestFilter1")
{
  moduleToTest = TestModuleBuilder::createTestFilter1(name, 0, 0, false);
}
void TestFilter1Test::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forTestFilter1(numFramesToProcess, inputs[v][0], inputs[v][1], inputs[v][2], inputs[v][3], 
                                          desiredOutputs[v][0], desiredOutputs[v][1], desiredOutputs[v][2]);
}




BiquadMacroTest::BiquadMacroTest()
: ProcessingTest("BiquadMacro")
{
  tolerance    = 1.e-14;
  moduleToTest = TestModuleBuilder::createBiquadMacro(name, 0, 0, false);
}
void BiquadMacroTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forBiquad(numFramesToProcess, inputs[v][0], inputs[v][1], inputs[v][2], inputs[v][3], inputs[v][4], inputs[v][5], 
                                     desiredOutputs[v][0]);
}

BiquadAtomicTest::BiquadAtomicTest()
: ProcessingTest("BiquadAtomic")
{
  tolerance    = 1.e-14;
  moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD);
}
void BiquadAtomicTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++)
    GenerateDesiredOutput::forBiquad(numFramesToProcess, inputs[v][0], inputs[v][1], inputs[v][2], inputs[v][3], inputs[v][4], inputs[v][5], 
                                     desiredOutputs[v][0]);
}


BlipTest::BlipTest()
: ProcessingTest("Blip")
{
  tolerance          = 1.e-14;
  moduleToTest       = TestModuleBuilder::createBlip("Blip", 0, 0, false);
  numFramesToProcess = maxNumFrames;

  //events             = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, maxNumFrames-1);

  events             = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, 3, 12);


}
void BlipTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  int octaveMultiplier = 1;
  for(int v = 0; v < numVoicesToUse; v++)
  {
    GenerateDesiredOutput::forFilterBlip(numFramesToProcess, octaveMultiplier*880.0, 20.0, desiredOutputs[v][0]);
    octaveMultiplier *= 2;
  }
}


MonoToPolyTest::MonoToPolyTest()
: ProcessingTest("MonoToPoly")
{
  moduleToTest = TestModuleBuilder::createMonoToPoly("MonoToPoly", 0, 0, false);
}
void MonoToPolyTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  for(int v = 0; v < numVoicesToUse; v++) 
    fillWithValue(desiredOutputs[v][0], numFramesToProcess, -1.0);
}



VoiceCombinerTest::VoiceCombinerTest()
: ProcessingTest("VoiceCombiner")
{
  moduleToTest = TestModuleBuilder::createVoiceCombiner("VoiceCombiner", 0, 0, false);
}
void VoiceCombinerTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  if( testModuleIsPolyphonic )
  {
    // 0-th voice output contains the sum ( == numvoicesToUse), other voices are referred to the same value (because our output
    // is monophonic) - so all voices contain numVoicesToUse in this case:
    for(int v = 0; v < numVoicesToUse; v++) 
      fillWithValue(desiredOutputs[v][0], numFramesToProcess, (double) numVoicesToUse);
  }
  else
  {
    // 0-th voice output contains the passed through 0-th voice (unity), other voices are referred to the same value (because our output
    // is monophonic) - so all voices contaim unity in this case:
    for(int v = 0; v < numVoicesToUse; v++) 
      fillWithValue(desiredOutputs[v][0], numFramesToProcess, 1.0);
  }
}

GateAndKillTest::GateAndKillTest()
: ProcessingTest("GateAndKill")
{
  moduleToTest       = TestModuleBuilder::createGateAndKill("GateAndKill", 0, 0, true);
  numFramesToProcess = 400;
  
  events = TestEventGenerator::generateNoteOnOffPair(1, 64, 10, 100);
  //events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(2, 64,  25, 100));
  //events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(3, 64,  50, 100));
  //events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(4, 64, 120, 100));
}
void GateAndKillTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{

}
/*
void GateAndKillTest::runTest()
{

}
*/
