#include "Polygon2DTests.h"


bool testPolygon2D()
{
  std::string testName = "Polygon2D";
  std::string dummy; 
  bool r = true;


  //r &= testRegularPolygonCreation2D(dummy);
  //r &= testPointInsidePolygon2D(dummy);
  r &= convexPolygonClipping(dummy);

  //appendTestResultToReport(reportString, testName, testResult);
  return r;
}

bool testRegularPolygonCreation2D(std::string &reportString)
{
  std::string testName = "RegularPolygonCreation2D";
  bool testResult = true;

  /*
  typedef rsPoint2D<double>   Point;
  typedef rsPolygon2D<double> Polygon;

  Point c1(0, 0);

  // equilateral triangle centered at origin, pointing to the right:
  Polygon p1(3, 1, c1, 0);
  */

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testPointInsidePolygon2D(std::string &reportString)
{
  std::string testName = "PointInsidePolygon2D";
  bool testResult = true;

  /*
  typedef rsPoint2D<double>   Point;
  typedef rsPolygon2D<double> Polygon;


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
  */

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


bool convexPolygonClipping(std::string &reportString)
{
  bool r = true;

  typedef rsVector2DF Vec2;       // 2D vector (for poylgon vertices
  typedef std::vector<Vec2> Poly; // a polygon is an array of vertices

  // a square:
  Vec2 
    s0(1.f, 1.f), 
    s1(5.f, 1.f),
    s2(5.f, 5.f),
    s3(1.0, 5.f);
  Poly square = { s0, s1, s2, s3 };  // clockwise
  //Poly square = { s3, s2, s1, s0 };    // counter clockwise

  // a triangle:
  Vec2 
    t0(  3.0f,  0.0f), 
    t1(  7.5f,  4.5f), 
    t2(  0.0f,  3.0f);
  Poly triangle = { t0, t1, t2 };  // clockwise
  //Poly triangle = { t2, t1, t0 };  // counter clockwise


  // test edge function:
  float d;
  d = edgeFunction(t0, t2, s0);        r &= d >  0;
  d = edgeFunction(t0, t2, t1);        r &= d <  0;
  d = edgeFunction(t0, t2, Vec2(2,1)); r &= d == 0;
  d = edgeFunction(t0, t2, Vec2(1,2)); r &= d == 0;

  
  // test line intersection:
  Vec2 v = lineIntersection(Vec2(3,0), Vec2(4,1), Vec2(0,3), Vec2(2.5,3.5));
  r &= v == Vec2(7.5, 4.5);

  Poly clipped = clipConvexPolygons(square, triangle); 
   // still wrong - but some vertices have the right coords already

  // should result in:  (2,1), (4,1), (5,2), (5,4), (1,3.4), (1,2)

  // doesn't work - problem with vertex order and/or definition of inside/orientation?
  


  // clipping should be commutative, i.e. the result independent from order of arguments
  // -> test this






  return r;
}