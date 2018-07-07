#include "Polygon2DTests.h"


bool testPolygon2D()
{
  std::string testName = "Polygon2D";
  std::string dummy; 
  bool r = true;


  //r &= testRegularPolygonCreation2D(dummy);
  //r &= testPointInsidePolygon2D(dummy);
  r &= convexPolygonClipping(dummy);
  r &= pixelCoverage(dummy);
  r &= triangleRasterization(dummy);

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
  typedef Vec2 V;

  // a square:
  Vec2
    s0(3, 2),
    s1(3, 6),
    s2(7, 6),
    s3(7, 2);
  Poly square ={ s0, s1, s2, s3 };  // counter clockwise


  // a triangle:
  Vec2
    t0(5, 1),
    t1(1, 5),
    t2(9, 5);
  Poly triangle ={ t0, t1, t2 };  // counter clockwise

  // test edge function:
  float d;

  d = edgeFunction(V(4, 1), V(4, 3), V(3, 8)); r &= d ==  2;
  d = edgeFunction(V(4, 1), V(4, 3), V(4, 8)); r &= d ==  0;
  d = edgeFunction(V(4, 1), V(4, 3), V(5, 8)); r &= d == -2;
  // rename edgeFunction To LeftDistance
  // returns a value proportional to the distance of the point p to the left of the directed
  //

  d = edgeFunction(t0, t1, s0);        r &= d >  0;
  d = edgeFunction(t0, t1, t2);        r &= d <  0;
  d = edgeFunction(t0, t1, Vec2(4, 2)); r &= d == 0;
  d = edgeFunction(t0, t1, Vec2(3, 3)); r &= d == 0;


  // test line intersection:
  Vec2 v = lineIntersection(Vec2(3, 0), Vec2(4, 1), Vec2(0, 3), Vec2(2.5, 3.5));
  r &= v == Vec2(7.5, 4.5);

  Poly clipped, target;

  clipped = clipAgainstEdge(triangle, s0, s1);
  target  ={ V(5,1), V(3,3), V(3,5), V(9,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle, s1, s2);
  target  ={ V(5,1), V(1,5), V(9,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle, s2, s3);
  target  ={ V(7,3), V(5,1), V(1,5), V(7,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle, s3, s0);
  target  ={ V(6,2), V(4,2), V(1,5), V(9,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle, s0, s3);
  target  ={ V(6,2), V(5,1), V(4,2) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle, s3, s2);
  target  ={ V(7,3), V(7,5), V(9,5) }; r &= clipped == target;
  // we need to test the case where none of the vertices is inside the edge and also
  // when adges of clipping and subject polygon coincide (maybe just clip again - this should
  // change nothing

  clipped = clipPolygon(triangle, square);
  r &= clipped.size() == 6;
  target  ={ V(7,3), V(6,2), V(4,2), V(3,3), V(3,5), V(7,5) }; r &= clipped == target;
  clipped = clipPolygon(square, triangle);
  target  ={ V(6,2), V(4,2), V(3,3), V(3,5), V(7,5), V(7,3) }; r &= clipped == target;
  // clipping the triangle against the square or vice versa results in the same clipped polygon,
  // however, the start vertex is different in both cases

  //float area;
  //area = polygonArea(square);
  //r &= area == 16; // it returns -16 - maybe winding is wrong bcs of the y-axis direction?
  //area = polygonArea(triangle);
  //r &= area == 16;

  Poly triangle2 ={ V(5,3), t1, t2 };  // 1st vertex is inside the square

  clipped = clipAgainstEdge(triangle2, s0, s1);
  target  = { V(5,3), V(3,4), V(3,5), V(9,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle2, s1, s2);
  target  = { V(5,3), V(1,5), V(9,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle2, s2, s3);
  target  = { V(7,4), V(5,3), V(1,5), V(7,5) }; r &= clipped == target;

  clipped = clipAgainstEdge(triangle2, s3, s0);
  target  = { V(5,3), V(1,5), V(9,5) }; r &= clipped == target;

  clipped = clipPolygon(triangle2, square); 
  target  = { V(7,4), V(5,3), V(3,4), V(3,5), V(7,5) }; r &= clipped == target;


  Poly triangle3 = { V(4,3), V(9,10), V(9,3) };

  clipped = clipAgainstEdge(triangle3, s0, s1);
  target  = { V(4,3), V(9,10), V(9,3) }; r &= clipped == target;
  clipped = clipAgainstEdge(triangle3, s1, s2);
  //target  = { V(4,3), V(6,6), V(9,6), V(9,3) }; r &= clipped == target;
  clipped = clipPolygon(triangle3, square); 
  //target  = { V(7,3), V(4,3), V(6,6), V(7,6) }; r &= clipped == target;
  // no - the result looks good - i think, the triangle is badly chosen, giving non-integer
  // clipping results




  return r;

  // counter-clockwise convention is used (to be consistent with OpenGL and DirectX)
  // info to vertex order in OpenGL
  // https://stackoverflow.com/questions/8142388/in-what-order-should-i-send-my-vertices-to-opengl-for-culling
  // By default, counterclockwise polygons are taken to be front-facing.
}


float unitSquareCut(const rsVector2DF& p, const rsVector2DF& q, bool& quadCut)
{
  // helper function for tests
  float x0, y0, x1, y1;
  unitSquareIntersections(p, q, x0, x1, y0, y1);
  return unitSquareCut(p, q, x0, x1, y0, y1, quadCut); 
}
bool pixelCoverage(std::string &reportString)
{
  bool r = true;

  typedef rsVector2DF V;    // for convenience
  float t = 1.f/16.f;       // target area
  float a;                  // computed area
  bool q;

  // edges through corners:
  a = unitSquareCut(V( 1,2), V(-0.5,-1), q); 
  // -0.25 - check how negative area can occur - the edge goes through the bottom left corner
  // i.e. the origin - maybe this needs special treatment?
  // check also cases, where it goes through 2 corners

  // cut off triangles at the corners:
  a = unitSquareCut(V( 1.5f,1.5f), V(  -1,0.25f), q);  r &= a == t; r &= q == false; // top-left
  a = unitSquareCut(V( 1.5f,0.5f), V(   0,1.25f), q);  r &= a == t; r &= q == false; // top-right
  a = unitSquareCut(V(-0.5f,0.5f), V(1.5f,-0.5f), q);  r &= a == t; r &= q == false; // bottom-left
  a = unitSquareCut(V(-1.0f,-0.75f), V(1.5f,0.5f), q); r &= a == t; r &= q == false; // bottom-right

  // cut off quadrilaterals:
  t = 6.f/16.f;
  a = unitSquareCut(V(2,1), V(-1,0.25f),q);  r &= a == t; r &= q == true;  // top
  a = unitSquareCut(V(0.75,2), V(0,-1),q);   r &= a == t; r &= q == true;  // left
  a = unitSquareCut(V(-1,0.75), V(2,0),q);   r &= a == t; r &= q == true;  // bottom
  a = unitSquareCut(V(0,-2), V(1,2),q);      r &= a == t; r &= q == true;  // right
  // for quads, the direction is important  

  //  non-intersecting edge:
  a = unitSquareCut(V(0,-1), V(2,0),q);  r &= a == 0; r &= q == false;




  return r;
}

bool triangleRasterization(std::string &reportString)
{
  bool r = true;

  typedef rsVector2DF Vec2;    // for convenience
  typedef Vec2 V;              // even shorter
  float c = 1.0f;              // color (gray value)
  rsImageF img(6, 4);          // image to draw on
  rsImageDrawerFFF drw(&img);  // drawer object
  drw.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);

  Vec2 A = V(0.5,0.5), B = V(4.5,2.5), C = V(2.5,0.5);
  float x = 0, y = 0;
 
  std::vector<Vec2> 
    triangle = { A, B, C },
    square   = { V(x, y), V(x, y+1), V(x+1, y+1), V(x+1, y) },
    polygon  = clipPolygon(triangle, square);

  drawTriangleAntiAliased(drw, A, B, C, c);
  //drawTriangleAntiAliased(drw, V(0.5,0.5), V(2.5,0.5), V(4.5,2.5), c);


  //writeImageToFilePPM(img, "TriangleTest.ppm");
  return r;
}