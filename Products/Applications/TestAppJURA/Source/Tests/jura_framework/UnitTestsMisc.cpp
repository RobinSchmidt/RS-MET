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

  using HSL = jura::ColourAHSL;
  using RGB = juce::Colour;


  auto equal = [](RGB c1, RGB c2) { return c1 == c2; };

  // The target RGB colors in the following tests were obtained by running the code on a windows 
  // machine on which the conversion is known to work well:
  HSL hsl(0.0f, 0.5f, 0.5f);
  RGB rgb = hsl.getAsJuceColour();          // b=63, g=63, r=191
  expect(equal(rgb, RGB(191, 63, 63)));

  //hsl.set(0.1f, 0.5f, 0.5f);  // todo: implement set()

  int dummy = 0;
}