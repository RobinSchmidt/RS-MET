#include "romos_ConcreteModularSystemTests.h"
using namespace romos;




BypassTest::BypassTest(const char *testName)
: ModularSystemTest(testName)
{    
  tolerance = 1.e-7;  // weird - this high tolerance is needed in release builds - why so huge? it's just a bypass....
}
void BypassTest::createAndConnectTestChildModules()
{
  topLevelModule->addAudioConnection(inModuleL, 0, outModuleL, 0);
  topLevelModule->addAudioConnection(inModuleR, 0, outModuleR, 0);

}
void BypassTest::fillDesiredOutputSignalArrays()
{
  RAPT::rsArrayTools::copy(inputs[0], desiredOutputs[0], signalLength);
  RAPT::rsArrayTools::copy(inputs[1], desiredOutputs[1], signalLength);
}



BypassWithChildTest::BypassWithChildTest()
: BypassTest("BypassWithChildTest")
{  

}
void BypassWithChildTest::createAndConnectTestChildModules()
{
  BypassTest::createAndConnectTestChildModules();
  //topLevelModule->addChildModule(getTypeId("UnitDelay"), "D",  10,  4, false, true);
  topLevelModule->addChildModule("UnitDelay", "D",  10,  4, false, true);
}





PolyBlipStereoTest::PolyBlipStereoTest()
: ModularSystemTest("PolyBlipStereoTest")
{  
  //events = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToUse, 12);
  events    = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, signalLength-1, 1, 12);
  tolerance = 1.e-8;
}
void PolyBlipStereoTest::createAndConnectTestChildModules()
{
  polyBlipStereo = (ContainerModule*) TestModuleBuilder::createPolyBlipStereo("PolyBlipStereo", 10, 10, false);

  topLevelModule->addChildModule(polyBlipStereo, true);
  topLevelModule->addAudioConnection(polyBlipStereo, 0, outModuleL, 0);
  topLevelModule->addAudioConnection(polyBlipStereo, 1, outModuleR, 0);

  // this action is here to expose a bug that occurred in ModuleConatiner updateHasDelayedConnectionsFlag() - the flag was set to true even 
  // when weren't any delayed connections:
  ContainerModule *in1Out2 = (ContainerModule*) polyBlipStereo->getChildModule(1);
  in1Out2->sortChildModuleArray();
}
void PolyBlipStereoTest::fillDesiredOutputSignalArrays()
{
  GenerateDesiredOutput::forFilterBlip(signalLength, 880.0, 20.0, desiredOutputs[0]);
  GenerateDesiredOutput::forFilterBlip(signalLength, 880.0, 20.0, desiredOutputs[1]);
}


NoiseFluteTest::NoiseFluteTest()
: ModularSystemTest("NoiseFluteTest")
{  
  //events = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, maxNumFrames-1, numVoicesToUse, 12);
  events    = TestEventGenerator::generateSimultaneousNotes(81, 64, 0, signalLength, 1, 12);  
  tolerance = 2.e-8;
}
void NoiseFluteTest::createAndConnectTestChildModules()
{
  noiseFlute = (ContainerModule*) TestModuleBuilder::createNoiseFlute("NoiseFlute", 10, 10, false);

  topLevelModule->addChildModule(noiseFlute, true);
  topLevelModule->addAudioConnection(noiseFlute, 0, outModuleL, 0);
  topLevelModule->addAudioConnection(noiseFlute, 1, outModuleR, 0);

}
void NoiseFluteTest::fillDesiredOutputSignalArrays()
{
  GenerateDesiredOutput::forWhiteNoiseUniform(signalLength, desiredOutputs[0], 0);
  GenerateDesiredOutput::forWhiteNoiseUniform(signalLength, desiredOutputs[1], 1);
  double coeffs[5];  
  romos::biquadBandpassConstSkirtCoeffs(coeffs, 880.0, 50.0);
  GenerateDesiredOutput::forBiquadWithFixedCoeffs(signalLength, desiredOutputs[0], coeffs[0],  coeffs[1],  coeffs[2],  coeffs[3],  
                                                  coeffs[4], desiredOutputs[0]);
  GenerateDesiredOutput::forBiquadWithFixedCoeffs(signalLength, desiredOutputs[1], coeffs[0],  coeffs[1],  coeffs[2],  coeffs[3],  
                                                  coeffs[4], desiredOutputs[1]);
  RAPT::rsArrayTools::scale(desiredOutputs[0], desiredOutputs[0], signalLength, 0.125);
  RAPT::rsArrayTools::scale(desiredOutputs[1], desiredOutputs[1], signalLength, 0.125);
}

