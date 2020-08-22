#include "romos_ContainerManipulationTests.h"
//using namespace romos;
//using namespace rsTestRomos;

namespace rsTestRomos
{

Containerize01::Containerize01() : ContainerManipulationTest("ContainerizationTest01")
{
  moduleToTest = TestModuleBuilder::createContainerize01("ContainerizationTest01", 0, 0, false);
}
bool Containerize01::runTest()
{
  bool result = true;

  // retrieve pointers to the child modules:
  romos::ContainerModule* containerLevel0  = (romos::ContainerModule*) moduleToTest;
  romos::Module* inputLevel0               = containerLevel0->getChildModule(0);
  romos::Module* minus                     = containerLevel0->getChildModule(1);
  romos::Module* output1Level0             = containerLevel0->getChildModule(2);
  romos::Module* output2Level0             = containerLevel0->getChildModule(3);

  // check if child modules have the right types:
  result &= inputLevel0->getTypeName() == "AudioInput";
  result &= minus->getTypeName() == "UnaryMinus";
  result &= output1Level0->getTypeName() == "AudioOutput";
  result &= output2Level0->getTypeName() == "AudioOutput";

  // we try to containerize them all -> should only containerize the UnaryMinus (I/O modules are excluded from containerizing):
  std::vector<romos::Module*> modulesToContainerize = containerLevel0->getChildModules();
  containerLevel0->containerizeModules(modulesToContainerize);

  romos::ContainerModule* containerLevel1 = (romos::ContainerModule*) containerLevel0->getChildModule(1);
  //result &= containerLevel1->getTypeIdentifierOld()  == ModuleTypeRegistry::CONTAINER;
  result &= containerLevel1->getTypeName() == "Container";

  // check number of children for inner and outer container
  result &= containerLevel0->getNumChildModules() == 4;
  result &= containerLevel1->getNumChildModules() == 3;

  // check if parent-relation is correct:
  result &= containerLevel1->getParentModule()    == containerLevel0;
  result &= inputLevel0->getParentModule()    == containerLevel0;
  result &= minus->getParentModule()    == containerLevel1;
  result &= output1Level0->getParentModule()    == containerLevel0;
  result &= output2Level0->getParentModule()    == containerLevel0;

  // check the created container's external connections:
  result &= containerLevel1->inputPins[0].sourceModule == inputLevel0;
  result &= output1Level0->inputPins[0].sourceModule == containerLevel1;
  result &= output1Level0->inputPins[0].outputIndex  == 0;
  result &= output2Level0->inputPins[0].sourceModule == containerLevel1;
  result &= output2Level0->inputPins[0].outputIndex  == 0;

  // check the created container's internal connections:
  romos::Module* inputLevel1  = containerLevel1->getChildModule(0);
  romos::Module* outputLevel1 = containerLevel1->getChildModule(2);
  result &= minus->inputPins[0].sourceModule == inputLevel1;
  result &= outputLevel1->inputPins[0].sourceModule == minus;



  //printModuleStructure(moduleToTest);

  return result;
}


Containerize02::Containerize02() : ContainerManipulationTest("ContainerizationTest02")
{
  moduleToTest = TestModuleBuilder::createContainerize02("ContainerizationTest02", 0, 0, false);
}
bool Containerize02::runTest()
{
  bool result = true;

  // retrieve child pointers:
  romos::ContainerModule* containerLevel0  = (romos::ContainerModule*) moduleToTest;
  romos::Module* constant1                 = containerLevel0->getChildModule(0);
  romos::Module* constant2                 = containerLevel0->getChildModule(1);
  romos::Module* constant3                 = containerLevel0->getChildModule(2);
  romos::Module* biquadDesigner            = containerLevel0->getChildModule(3);
  romos::Module* adder                     = containerLevel0->getChildModule(4);
  romos::Module* multiplier                = containerLevel0->getChildModule(5);
  romos::Module* output1Level0             = containerLevel0->getChildModule(6);
  romos::Module* output2Level0             = containerLevel0->getChildModule(7);
  romos::Module* output3Level0             = containerLevel0->getChildModule(8);

  // containerize subset of the child modules:
  std::vector<romos::Module*> modulesToContainerize;
  modulesToContainerize.push_back(constant3);
  modulesToContainerize.push_back(biquadDesigner);
  modulesToContainerize.push_back(adder);
  containerLevel0->containerizeModules(modulesToContainerize);

  // retrieve the created container:
  romos::ContainerModule* containerLevel1 = (romos::ContainerModule*) containerLevel0->getChildModule(2);
  //result &= containerLevel1->getTypeIdentifierOld()  == ModuleTypeRegistry::CONTAINER;
  result &= containerLevel1->getTypeName()  == "Container";

  // check number of children for inner and outer container
  result &= containerLevel0->getNumChildModules() == 7;
  result &= containerLevel1->getNumChildModules() == 8;

  // retrieve I/O modules of the created container:
  romos::Module* input1Level1           = containerLevel1->getChildModule(0);
  romos::Module* input2Level1           = containerLevel1->getChildModule(1);
  romos::Module* output1Level1          = containerLevel1->getChildModule(5);
  romos::Module* output2Level1          = containerLevel1->getChildModule(6);
  romos::Module* output3Level1          = containerLevel1->getChildModule(7);

  // check, if outer container is parent of the modules when it should be:
  result &= constant1->getParentModule() == containerLevel0;
  result &= constant2->getParentModule() == containerLevel0;
  result &= containerLevel1->getParentModule() == containerLevel0;
  result &= multiplier->getParentModule() == containerLevel0;
  result &= output1Level0->getParentModule() == containerLevel0;
  result &= output2Level0->getParentModule() == containerLevel0;
  result &= output3Level0->getParentModule() == containerLevel0;

  // check, if inner container is parent of the modules when it should be:
  result &= constant3->getParentModule() == containerLevel1;
  result &= biquadDesigner->getParentModule() == containerLevel1;
  result &= adder->getParentModule() == containerLevel1;

  // check the created container's external connections:
  result &= containerLevel1->inputPins[0].sourceModule == constant1;
  result &= containerLevel1->inputPins[1].sourceModule == constant2;
  result &= multiplier->inputPins[0].sourceModule == containerLevel1;
  result &= multiplier->inputPins[0].outputIndex  == 2;
  result &= output1Level0->inputPins[0].sourceModule == containerLevel1;
  result &= output1Level0->inputPins[0].outputIndex  == 0;
  result &= output2Level0->inputPins[0].sourceModule == containerLevel1;
  result &= output2Level0->inputPins[0].outputIndex  == 1;
  result &= output3Level0->inputPins[0].sourceModule == constant2;
  result &= output3Level0->inputPins[0].outputIndex  == 0;

  // check the created container's internal connections:
  result &= biquadDesigner->inputPins[0].sourceModule  == input1Level1;
  result &= biquadDesigner->inputPins[1].sourceModule  == constant3;
  result &= biquadDesigner->inputPins[2].sourceModule  == input2Level1;
  result &= adder->inputPins[0].sourceModule  == biquadDesigner;
  result &= adder->inputPins[0].outputIndex   == 3;
  result &= adder->inputPins[0].outputPointer == biquadDesigner->audioOutputs + 3;
  result &= output1Level1->inputPins[0].sourceModule  == biquadDesigner;
  result &= output1Level1->inputPins[0].outputIndex   == 0;
  result &= output1Level1->inputPins[0].outputPointer == biquadDesigner->audioOutputs;
  result &= output2Level1->inputPins[0].sourceModule  == biquadDesigner;
  result &= output2Level1->inputPins[0].outputIndex   == 1;
  result &= output2Level1->inputPins[0].outputPointer == biquadDesigner->audioOutputs + 1;
  result &= output3Level1->inputPins[0].sourceModule  == biquadDesigner;
  result &= output3Level1->inputPins[0].outputIndex   == 2;
  result &= output3Level1->inputPins[0].outputPointer == biquadDesigner->audioOutputs + 2;

  //printModuleStructure(moduleToTest);

  return result;
}




OutputModuleDeletion::OutputModuleDeletion() : ContainerManipulationTest("OutputModuleDeletion")
{
  moduleToTest = TestModuleBuilder::createOutputDeletion("OutputDeletion", 0, 0, false);
}
bool OutputModuleDeletion::runTest()
{
  bool result = true;

  ContainerModule* outerContainer = (ContainerModule*)moduleToTest;
  ContainerModule* innerContainer = (ContainerModule*)outerContainer->getChildModule(0);
  romos::Module* out1           = innerContainer->getChildModule(3);
  romos::Module* out2           = innerContainer->getChildModule(4);
  romos::Module* out3           = innerContainer->getChildModule(5);
  romos::Module* minus1         = outerContainer->getChildModule(1);
  romos::Module* minus2         = outerContainer->getChildModule(2);
  romos::Module* minus3         = outerContainer->getChildModule(3);

  // delete the middle output module of the 3:
  innerContainer->deleteChildModule(out2);

  // check if connection of minus2 was disconnected:
  result &= minus2->inputPins[0].sourceModule      == NULL;
  result &= minus2->inputPins[0].outputIndex       == 0;
  result &= minus2->inputPins[0].outputFrameSize   == 0;
  result &= minus2->inputPins[0].outputVoiceStride == 0;
  result &= minus2->inputPins[0].outputPointer     == &(minus2->inputPins[0].defaultValue);

  // check if connection of minus3 was moved up:
  result &= minus3->inputPins[0].outputIndex       == 1;

  // delete the bottom output module of the remaining two:
  innerContainer->deleteChildModule(out3);

  // check if connection of minus3 was disconnected:
  result &= minus3->inputPins[0].sourceModule      == NULL;
  result &= minus3->inputPins[0].outputIndex       == 0;
  result &= minus3->inputPins[0].outputFrameSize   == 0;
  result &= minus3->inputPins[0].outputVoiceStride == 0;
  result &= minus3->inputPins[0].outputPointer     == &(minus3->inputPins[0].defaultValue);

  // delete the remaining output module:
  innerContainer->deleteChildModule(out1);

  // check if connection of minus1 was disconnected:
  result &= minus1->inputPins[0].sourceModule      == NULL;
  result &= minus1->inputPins[0].outputIndex       == 0;
  result &= minus1->inputPins[0].outputFrameSize   == 0;
  result &= minus1->inputPins[0].outputVoiceStride == 0;
  result &= minus1->inputPins[0].outputPointer     == &(minus1->inputPins[0].defaultValue);

  return result;
}













ContainerizationAddedConstantsTest::ContainerizationAddedConstantsTest()
  : ProcessingTest("ContainerizationAddedConstants")
{
  moduleToTest = TestModuleBuilder::createAddedConstants("AddedConstants", 0, 0, false);
}
void ContainerizationAddedConstantsTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  RAPT::rsArrayTools::fillWithValue(desiredOutputs[0][0], numFramesToProcess, 136.0);
}
bool ContainerizationAddedConstantsTest::runTest()
{
  initTest();
  fillDesiredOutputSignalArrays(false);


  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();
  if(!RAPT::rsArrayTools::equal(outputs[0][0], desiredOutputs[0][0], numFramesToProcess))
    return false;

  containerizeSum();
  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();
  if(!RAPT::rsArrayTools::equal(outputs[0][0], desiredOutputs[0][0], numFramesToProcess))
    return false;

  unContainerizeSum();
  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();
  if(!RAPT::rsArrayTools::equal(outputs[0][0], desiredOutputs[0][0], numFramesToProcess))
    return false;

  static const int numIterations = 50;
  bool verboseOutput = false;
  for(int i = 1; i <= numIterations; i++)
  {
    randomizeContainment();
    moduleToTest->resetStateForAllVoices();
    processModuleInFrames();

    //printModuleStructure(moduleToTest, 0);

    if(!RAPT::rsArrayTools::equal(outputs[0][0], desiredOutputs[0][0], numFramesToProcess))
    {
      printModuleStructure(moduleToTest, 0);
      RAPT::rsAssert(false, "plotting code needs update");
      //Plotter::plotData(numFramesToProcess, timeAxis, desiredOutputs[0][0], outputs[0][0]);
      return false;
    }
  }

  return true;
}
void ContainerizationAddedConstantsTest::containerizeSum()
{
  romos::ContainerModule* container = dynamic_cast<romos::ContainerModule*> (moduleToTest);
  std::vector<romos::Module*> childModules = container->getChildModules();
  std::vector<romos::Module*> sumInVector;
  sumInVector.push_back(childModules[16]);
  container->containerizeModules(sumInVector);
}
void ContainerizationAddedConstantsTest::unContainerizeSum()
{
  romos::ContainerModule* container = dynamic_cast<romos::ContainerModule*> (moduleToTest);
  std::vector<romos::Module*> childModules = container->getChildModules();
  std::vector<romos::Module*> containerizedSumInVector;
  containerizedSumInVector.push_back(childModules[16]);
  container->unContainerizeModules(containerizedSumInVector);
}
void ContainerizationAddedConstantsTest::randomizeContainment()
{
  romos::ContainerModule* container = dynamic_cast<romos::ContainerModule*> (moduleToTest);
  if(container != NULL)
  {
    std::vector<romos::Module*> childModules    = container->getChildModules();
    std::vector<romos::Module*> toBeContainerized;
    std::vector<romos::Module*> toBeUnContainerized;
    for(unsigned int i = 0; i < childModules.size(); i++)
    {
      //if( childModules[i]->getTypeIdentifierOld() == ModuleTypeRegistry::CONSTANT || childModules[i]->isContainerModule() )
      if(childModules[i]->getTypeName() == "Constant" || childModules[i]->isContainerModule())
      {
        double randomNumber = rsTestRomos::random(0.0, 1.0);
        if(randomNumber < 1.0/3.0)
          toBeContainerized.push_back(childModules[i]);
        else if(randomNumber < 2.0/3.0)
          toBeUnContainerized.push_back(childModules[i]);
        else
        {
          // do nothing
        }
      }
    }
    container->containerizeModules(toBeContainerized);
    container->unContainerizeModules(toBeUnContainerized);
  }
}


PinSortingTest::PinSortingTest()
  : ProcessingTest("PinSorting")
{
  moduleToTest = TestModuleBuilder::createPinSortTest("PinSorting", 0, 0, false);
}
void PinSortingTest::fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic)
{
  processModuleInFrames();
  RAPT::rsArrayTools::copy(outputs[0][0], desiredOutputs[0][0], numFramesToProcess);
  RAPT::rsArrayTools::copy(outputs[0][1], desiredOutputs[0][1], numFramesToProcess);
  RAPT::rsArrayTools::copy(outputs[0][2], desiredOutputs[0][2], numFramesToProcess);
}
bool PinSortingTest::runTest()
{
  initTest();
  //fillDesiredOutputSignalArrays();

  // retrieve pointers to some embedded modules (assumes certain order):
  romos::ContainerModule* inner = dynamic_cast<romos::ContainerModule*> (((ContainerModule*)moduleToTest)->getChildModule(3));
  rassert(inner != NULL);
  romos::Module* in1   = inner->getAudioInputModule(0);
  romos::Module* in2   = inner->getAudioInputModule(1);
  romos::Module* in3   = inner->getAudioInputModule(2);
  romos::Module* out1  = inner->getAudioOutputModule(0);
  romos::Module* out2  = inner->getAudioOutputModule(1);
  romos::Module* out3  = inner->getAudioOutputModule(2);

  bool result = true;

  //Plotter::plotData(N, t, d[0][0], y[0][0]);

  // test all possible permutations of the input order (wrap into function):
  exchangeModulePositions(in1, in2);      // input order: 1,2,3 -> 2,1,3
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(in1, in3);      // input order: 2,1,3 -> 2,3,1
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(in1, in2);      // input order: 2,3,1 -> 1,3,2
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(in1, in3);      // input order: 1,3,2 -> 3,1,2
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(in1, in2);      // input order: 3,1,2 -> 3,2,1
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(in1, in3);      // input order: 3,2,1 -> 1,2,3 (again)

  // test all possible permutations of the output order (wrap into function):
  exchangeModulePositions(out1, out2);      // output order: 1,2,3 -> 2,1,3
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(out1, out3);      // output order: 2,1,3 -> 2,3,1
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(out1, out2);      // output order: 2,3,1 -> 1,3,2
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(out1, out3);      // output order: 1,3,2 -> 3,1,2
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(out1, out2);      // output order: 3,1,2 -> 3,2,1
  processModuleInFrames();
  result &= doOutputsMatchDesiredOutputs(false);
  exchangeModulePositions(out1, out3);      // output order: 3,2,1 -> 1,2,3 (again)

  return result;
}
void PinSortingTest::exchangeModulePositions(romos::Module* module1, romos::Module* module2)
{
  int x1 = module1->getPositionX();
  int y1 = module1->getPositionY();
  module1->setPositionXY(module2->getPositionX(), module2->getPositionY());
  module2->setPositionXY(x1, y1);
}

}