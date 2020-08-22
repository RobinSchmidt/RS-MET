//#include "PerformanceTests.h"
#include "TestHelpers.h"
using namespace rosic;
using namespace romos;


void testCodeGenerator()
{
  // here we choose what kind of module to create (and genreate code for):
  //romos::Module *testModule = createSumDiffModule(NULL);
  //romos::Module *testModule = createWrappedSumDiffModule(NULL);
  romos::Module *testModule = rsTestRomos::TestModuleBuilder::createTestFilter1("TestFilter1", 0, 0, false);




  rosic::rsString codeForModule = rsTestRomos::ModuleBuildCodeGenerator::getCodeForModule(testModule);
  codeForModule.printToStandardOutput();
  moduleFactory.deleteModule(testModule);
}
