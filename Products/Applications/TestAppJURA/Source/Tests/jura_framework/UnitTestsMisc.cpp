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

  // Tests, if the given hsl values convert correctly to the given target rgb values:
  auto testHSL2RGB = [&](float h, float s, float l, uint8 r, uint8 g, uint8 b)
  {
    HSL hsl(h, s, l);
    RGB rgb = hsl.getAsJuceColour();
    expect(rgb == RGB(r, g, b));
  };

  // Create HSL colors with S=0.5, L=0.5 and H going from 0 to 1 in 0.1 steps and check the result
  // of converting to RGB as juce::Color object. The target RGB colors in the following tests were 
  // obtained by running the code on a windows machine on which the conversion is known to work:
  testHSL2RGB(0.0f, 0.5f, 0.5f, 191,  63,  63);
  testHSL2RGB(0.1f, 0.5f, 0.5f, 191, 140,  63);
  testHSL2RGB(0.2f, 0.5f, 0.5f, 165, 191,  63);
  testHSL2RGB(0.3f, 0.5f, 0.5f,  89, 191,  63);
  testHSL2RGB(0.4f, 0.5f, 0.5f,  63, 191, 114);
  testHSL2RGB(0.5f, 0.5f, 0.5f,  63, 191, 191);
  testHSL2RGB(0.6f, 0.5f, 0.5f,  63, 114, 191);
  testHSL2RGB(0.7f, 0.5f, 0.5f,  89,  63, 191);
  testHSL2RGB(0.8f, 0.5f, 0.5f, 165,  63, 191);
  testHSL2RGB(0.9f, 0.5f, 0.5f, 191,  63, 140);
  testHSL2RGB(1.0f, 0.5f, 0.5f, 191,  63,  63);

  int dummy = 0;
}