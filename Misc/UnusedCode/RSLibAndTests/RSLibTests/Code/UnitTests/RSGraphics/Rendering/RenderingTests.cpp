#include "RenderingTests.h"

bool testRendering(std::string &reportString)
{
  std::string testName = "Rendering";
  bool testResult = true;

  testResult &= testDrawOutlinedPixelRectangle(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testDrawOutlinedPixelRectangle(std::string &reportString)
{
  std::string testName = "DrawOutlinedPixelRectangle";
  bool testResult = true;

  rsImageRGBA image(40, 30, true, rsColorRGBA(0, 0, 0, 0));
  rsImageRegionRGBA imageRegion(&image, 0, 0, image.getWidth(), image.getHeight());
  rsGraphicsRenderer2DImage r(imageRegion);


  r.drawOutlinedPixelRectangle(10, 10, 100, 30);

  //renderer.drawOutlinedPixelRectangle(0, 0, image.getWidth(),   image.getHeight());
  //renderer.drawOutlinedPixelRectangle(1, 1, image.getWidth()-1, image.getHeight()-1);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
