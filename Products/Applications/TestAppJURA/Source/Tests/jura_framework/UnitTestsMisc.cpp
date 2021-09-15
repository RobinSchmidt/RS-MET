#include "UnitTestsMisc.h"
using namespace juce;
using namespace jura;

UnitTestMisc::UnitTestMisc() 
  : juce::UnitTest("Misc", "Misc")
{

}

void UnitTestMisc::runTest()
{
  runTestColor();
}

void UnitTestMisc::runTestColor()
{
  beginTest("Color");
}