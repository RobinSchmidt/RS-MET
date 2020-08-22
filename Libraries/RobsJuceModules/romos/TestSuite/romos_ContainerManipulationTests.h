#ifndef romos_ContainerManipulationTests_h
#define romos_ContainerManipulationTests_h

//#include "romos_ProcessingTest.h"
//#include "romos_TestModuleBuilder.h"
//#include "romos_TestEventGenerator.h"
//#include "romos_GenerateDesiredOutput.h"


namespace rsTestRomos
{

  /**

  This file contains test-classes that involve manipulations of container modules such as moving around child-modules, 
  containerizing/uncontainerizing, etc. 

  */

  class ContainerManipulationTest : public UnitTest
  {
  public:
    ContainerManipulationTest(const char *testName) : UnitTest(testName) 
    { 
      moduleToTest = NULL;
    }
    virtual ~ContainerManipulationTest()
    {
      romos::moduleFactory.deleteModule(moduleToTest);
    }

  protected:
    romos::Module *moduleToTest;

    // maybe have a data structure to store the structure of the moduleToTest
    // -> containerize -> uncontainerize -> compare
  };


  /** Uses a container that has an Input connected to a UnaryMinus connected to two Outputs - does some conatinerization actions
  and checks if the resulting structure is as desired. */
  class Containerize01 : public ContainerManipulationTest
  {
  public:
    Containerize01();
  protected:
    virtual bool runTest();
  };

  /** Does containerization for a slightly more compicated case. */
  class Containerize02 : public ContainerManipulationTest
  {
  public:
    Containerize02();
  protected:
    virtual bool runTest();
  };



  /** Creates a container module with 3 constants that go to 3 output modules and an outlying container that picks up these three output 
  signals (feeding them into UnaryMinus modules). In the test, the outputs of the inner container are deleted and we check if the minusses 
  in the outer container properly disconnect themselves. */
  class OutputModuleDeletion : public ContainerManipulationTest
  {
  public:
    OutputModuleDeletion();
  protected:
    virtual bool runTest();
  };



  /** Uses a module that contains a bunch of constant modules that are summed at the output. Within a loop, a pseudo-random subset of the 
  child-modules is selected and containerized and another subset is uncontainerized. The output sum should always stay the same regardless 
  of the containment structure. */
  class ContainerizationAddedConstantsTest : public ProcessingTest
  {
  public:
    ContainerizationAddedConstantsTest();
  protected:
    virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
    virtual bool runTest();
    virtual void containerizeSum();
    virtual void unContainerizeSum();
    virtual void randomizeContainment();
  };

  /** Uses a module that has an inner and an outer container. The input and output modules of the inner container are moved so as to change 
  the position of the corresponding pin of the inner container. After each such move, the outer container should update its connections 
  that go into/out-of the inner container so as to reflect the new ordering of the pins. The test assures that the output container 
  correctly does this connection update by verifying that the I/O function of the outer container remains the same, regardless of the order 
  of the inner container's I/O modules. */
  class PinSortingTest : public ProcessingTest
  {
  public:
    PinSortingTest();
  protected:
    virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic = false);
    virtual bool runTest();
    //virtual bool testInputPermutation();
    //virtual bool testOutputPermutation();
    virtual void exchangeModulePositions(romos::Module *module1, romos::Module *module2);
  };






} // end namespace romos

#endif 
