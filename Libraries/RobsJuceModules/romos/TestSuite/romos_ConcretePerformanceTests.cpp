#include "romos_ConcretePerformanceTests.h"
using namespace rsTestRomos;  // try to get rid



BiquadMacroPerformanceTest::BiquadMacroPerformanceTest() 
: PerformanceTest("BiquadMacro")
{
  moduleToTest = TestModuleBuilder::createBiquadMacro(name, 0, 0, false);
}


IdentityChainPerformanceTest::IdentityChainPerformanceTest() 
: PerformanceTest("IdentityChain")
{
  moduleToTest = TestModuleBuilder::createIdentityChain(name, 0, 0, false);
}


IdentityChainWithFeedbackPerformanceTest::IdentityChainWithFeedbackPerformanceTest() 
: PerformanceTest("IdentityChainWithFeedback")
{
  moduleToTest = TestModuleBuilder::createIdentityChainWithFeedback(name, 0, 0, false);
}


AdderChainPerformanceTest::AdderChainPerformanceTest() 
: PerformanceTest("AdderChain")
{
  moduleToTest = TestModuleBuilder::createAdderChain(name, 0, 0, false);
}


AdderChainWithFeedbackPerformanceTest::AdderChainWithFeedbackPerformanceTest() 
: PerformanceTest("AdderChainWithFeedback")
{
  moduleToTest = TestModuleBuilder::createAdderChainWithFeedback(name, 0, 0, false);
}




FirstOrderFilterPerformanceTest::FirstOrderFilterPerformanceTest() 
: PerformanceTest("FirstOrderFilter")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::FIRST_ORDER_FILTER);
  moduleToTest = moduleFactory.createModule("FirstOrderFilter");
}

BiquadPerformanceTest::BiquadPerformanceTest() 
: PerformanceTest("Biquad")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD);
  moduleToTest = moduleFactory.createModule("Biquad");
}

BiquadDesignerPerformanceTest::BiquadDesignerPerformanceTest() 
: PerformanceTest("BiquadDesigner")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD_DESIGNER);
  moduleToTest = moduleFactory.createModule("BiquadDesigner");
}

BandlimitedImpulseTrainPerformanceTest::BandlimitedImpulseTrainPerformanceTest() 
: PerformanceTest("BandlimitedImpulseTrain")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::BANDLIMITED_IMPULSE_TRAIN);
  moduleToTest = moduleFactory.createModule("BandlimitedImpulseTrain");
}

SawOscillatorPerformanceTest::SawOscillatorPerformanceTest() 
: PerformanceTest("SawOscillator")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::BLIT_SAW_OSCILLATOR);
  moduleToTest = moduleFactory.createModule("BlitOscillator");
}
