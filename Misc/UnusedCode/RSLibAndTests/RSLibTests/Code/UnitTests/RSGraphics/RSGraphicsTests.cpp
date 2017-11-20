#include "RSGraphicsTests.h"

#include "Fonts/FontTests.cpp"
#include "Rendering/RenderingTests.cpp"

// maybe move to somewhere else
bool testImage(std::string &reportString)
{
  std::string testName = "Image";
  bool testResult = true;


  rsImageRGBA img1(40, 30);
  testResult &= img1.getByteSize() == 40*30*4;

  // test copy constructor (should create a deep copy):
  rsImageRGBA img2(img1);
  testResult &= img2.getWidth()  == 40;
  testResult &= img2.getHeight() == 30;
  testResult &= img2.getPointerToPixel(0, 0) != img1.getPointerToPixel(0, 0);

  // test assignment operator (should also create a deep copy):
  rsImageRGBA img3 = img1;
  testResult &= img3.getWidth()  == 40;
  testResult &= img3.getHeight() == 30;
  testResult &= img3.getPointerToPixel(0, 0) != img1.getPointerToPixel(0, 0);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
