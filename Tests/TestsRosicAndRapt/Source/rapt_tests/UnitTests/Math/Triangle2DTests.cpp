#include "Triangle2DTests.h"

bool testTriangle2D()
{
  std::string testName = "Triangle2D";
  std::string dummy;
  bool testResult = true;

  testResult &= testPointInsideTriangle2D(dummy);

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testPointInsideTriangle2D(std::string &reportString)
{
  std::string testName = "PointInsideTriangle2D";
  bool testResult = true;

  //typedef rsPoint2D<double>    Point;
  typedef rsVector2D<double>   Point;
  typedef rsTriangle2D<double> Triangle;

  Point a(1, 1);
  Point b(6, 5);
  Point c(5, 0);
  Triangle t(a, b, c);

  // test some points inside the triangle:
  testResult &=  t.containsPoint( Point(1, 1) );
  testResult &=  t.containsPoint( Point(2, 1) );
  testResult &=  t.containsPoint( Point(3, 1) );
  testResult &=  t.containsPoint( Point(4, 1) );
  testResult &=  t.containsPoint( Point(3, 2) );
  testResult &=  t.containsPoint( Point(4, 2) );
  testResult &=  t.containsPoint( Point(5, 2) );
  testResult &=  t.containsPoint( Point(4, 3) );
  testResult &=  t.containsPoint( Point(5, 3) );
  testResult &=  t.containsPoint( Point(5, 4) );

  // test some points outside the triangle:
  testResult &= !t.containsPoint( Point(1, 2) );
  testResult &= !t.containsPoint( Point(2, 2) );
  testResult &= !t.containsPoint( Point(3, 3) );
  testResult &= !t.containsPoint( Point(4, 4) );
  testResult &= !t.containsPoint( Point(5, 5) );
  testResult &= !t.containsPoint( Point(6, 4) );
  testResult &= !t.containsPoint( Point(6, 3) );
  testResult &= !t.containsPoint( Point(6, 2) );
  testResult &= !t.containsPoint( Point(6, 1) );
  testResult &= !t.containsPoint( Point(4, 0) );
  testResult &= !t.containsPoint( Point(3, 0) );
  testResult &= !t.containsPoint( Point(2, 0) );
  testResult &= !t.containsPoint( Point(1, 0) );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
