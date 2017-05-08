#include "romos_UnitTest.h"
using namespace romos;

#include <cstring>
#include "stdio.h"

UnitTest::UnitTest(const char *testName)
{
  int nameLength = strlen(testName);
  name = new char[nameLength+1];
  strcpy(name, testName);
}

UnitTest::~UnitTest()
{
  delete[] name;
}

bool UnitTest::runTestAndPrintResultToConsole()
{
  printf("%s %s %s", "Running ", name, "...\n");



  bool testPassed = runTest();

  if( testPassed )
    printf("%s", " passed.\n");
  else
    printf("%s", " FAILED !!!\n");

  return testPassed;
}

void UnitTest::printModuleStructure(romos::Module *module, int indent)
{
  char *spaces = new char[indent+1];
  for(int i=0; i<indent+1; i++)
    spaces[i] = ' ';
  spaces[indent] = '\0';
  printf("%s", spaces);
  delete[] spaces;

  printf("%s", module->getName().getRawString() );
  /*
  printf("%s", " - #AI: " );
  printf("%d", module->getNumAudioInputs() );
  printf("%s", ", #AO: " );
  printf("%d", module->getNumAudioOutputs() );
  printf("%s", ", #EI: " );
  printf("%d", module->getNumEventInputs() );
  printf("%s", ", #EO: " );
  printf("%d", module->getNumEventOutputs() );
  */
  printf("%s", "\n" );

  romos::ModuleContainer *container = dynamic_cast<romos::ModuleContainer*> (module); 
  if( container != NULL )
  {
    for(unsigned int i=0; i<container->getNumChildModules(); i++)
      printModuleStructure(container->getChildModule(i), indent+1);
  }
}

