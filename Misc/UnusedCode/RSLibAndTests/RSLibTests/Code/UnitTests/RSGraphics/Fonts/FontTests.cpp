#include "FontTests.h"

bool testFonts(std::string &reportString)
{
  std::string testName = "Font";
  bool testResult = true;

  testResult &= testFontCreation(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFontCreation(std::string &reportString)
{
  std::string testName = "FontCreation";
  bool testResult = true;

  rsPixelFont *font;

  font = rsGlobalFontInstances::getPixelFontRoundedA7D0();
  testResult &= font->getFontHeight() == 7;

  font = rsGlobalFontInstances::getPixelFontRoundedBoldA10D0();
  testResult &= font->getFontHeight() == 10;

  font = rsGlobalFontInstances::getPixelFontRoundedBoldA16D0();
  testResult &= font->getFontHeight() == 16;

  rsGlobalFontInstances::cleanUpFonts();

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
