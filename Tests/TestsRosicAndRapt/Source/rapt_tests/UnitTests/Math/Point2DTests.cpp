
bool testPoint2DOperators(std::string &reportString)
{
  std::string testName = "rsPoint2DOperators";
  bool testResult = true;

  rsPoint2D<double> p0(2.0, 1.0);
  rsPoint2D<double> p1(2.0, 1.0);
  rsPoint2D<double> p2(4.0, 2.0);
  rsPoint2D<double> p3(5.0, 2.0);

  testResult &= (  p0 == p1 );
  testResult &= (  p1 != p3 );
  testResult &= ( -p1 == rsPoint2D<double>(-2.0, -1.0) );

  rsPoint2D<double> p4 = p1;  // (2,1)
  p4 += p2;                 // (2,1)+(4,2) = (6,3)
  testResult &= ( p4     == rsPoint2D<double>(6.0, 3.0) );
  p4 -= p2;                 // (6,3)-(4,2) = (2,1)
  testResult &= ( p4     == rsPoint2D<double>(2.0, 1.0) );
  p4 *= 2.0;                // (2,1)*2 = (4,2)
  testResult &= ( p4     == rsPoint2D<double>(4.0, 2.0) );
  p4 /= 2.0;                // (4,2)/2 = (2,1);
  testResult &= ( p4     == rsPoint2D<double>(2.0, 1.0) );

  testResult &= ( p4*2.0 == rsPoint2D<double>(4.0, 2.0) );  // (2,1)*2 = (4,2)
  testResult &= ( 2.0*p4 == rsPoint2D<double>(4.0, 2.0) );  // 2*(2,1) = (4,2)
  testResult &= ( p4/2.0 == rsPoint2D<double>(1.0, 0.5) );  // (2,1)/2 = (1,0.5)
  testResult &= ( p1+p2  == rsPoint2D<double>(6.0, 3.0) );  // (2,1)+(4,2) = (6,3)
  testResult &= ( p3-p2  == rsPoint2D<double>(1.0, 0.0) );  // (5,2)-(4,2) = (1,0)

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testPoint2DTransformations(std::string &reportString)
{
  std::string testName = "rsPoint2DTransformations";
  bool testResult = true;

  // some aliases to ease the template-instantiation:
  typedef rsPoint2D<double>           Point;
  typedef rsAffineTransform2D<double> Trafo;

  Point p1, p2, p3;              // points to be transformed
  Point q1, q2, q3;              // image points
  Trafo t1, t2, t3, t4, t5, tc;  // transformations, tc: composed/chained transformation

  // test primitive transforms:
  t1 = Trafo::rotation(PI/2);
  q1 = t1.getTransformedPoint(Point(1, 0));
  testResult &= ( q1.isCloseTo(Point(0, 1), 0.000001) );
  t1 = Trafo::rotation(PI/4);
  q1 = t1.getTransformedPoint(Point(1, 0));
  testResult &= ( q1.isCloseTo(Point(SQRT2_INV, SQRT2_INV), 0.000001) );
  t1 = Trafo::rotationByHalfPi();
  q1 = t1.getTransformedPoint(Point(1, 1));
  testResult &= ( q1  == Point(-1, 1) );
  t1 = Trafo::rotation(PI/2, Point(1, 1));
  q1 = t1.getTransformedPoint(Point(2, 1));
  testResult &= ( q1.isCloseTo(Point(1, 2), 0.000001) );
  t1 = Trafo::scaling(2, 3);
  q1 = t1.getTransformedPoint(Point(2, 3));
  testResult &= ( q1.isCloseTo(Point(4, 9), 0.000001) );
  t1 = Trafo::translation(2, 3);
  q1 = t1.getTransformedPoint(Point(2, 3));
  testResult &= ( q1.isCloseTo(Point(4, 6), 0.000001) );

  t1 = Trafo::reflection(-1, +1, +1);  // defines the line -1*x + 1*y + 1 = 0 -> y = x-1
  q1 = t1.getTransformedPoint(Point(2, -3));
  testResult &= ( q1  == Point(-2, 1) );

  t1 = Trafo::reflection(1, -1);       // should result in the same transfrom as above
  q1 = t1.getTransformedPoint(Point(2, -3));
  testResult &= ( q1  == Point(-2, 1) );

  // test chaining:
  t1 = Trafo::translation(1, 1); // translate x and y by 1
  t2 = Trafo::scaling(2, 2);     // scale x and y by 2
  q1 = t1.getTransformedPoint(Point(0, 0));
  testResult &= ( q1  == Point(1, 1) );
  q1 = t2.getTransformedPoint(q1);
  testResult &= ( q1  == Point(2, 2) );
  tc = t1.followedBy(t2);   // translate (0,0) to (1,1) then scale by 2
  q1 = tc.getTransformedPoint(Point(0, 0));
  testResult &= ( q1  == Point(2, 2) );
  tc = t2.followedBy(t1);   // scale (0,0) by 2, then translate to (1,1)
  q1 = tc.getTransformedPoint(Point(0, 0));
  testResult &= ( q1  == Point(1, 1) );
  t1 = Trafo::rotation(PI/2);
  t2 = Trafo::translation(1,2);
  t3 = Trafo::scaling(2,3);
  p1 = Point(2,3);
  q1 = t1.getTransformedPoint(p1);
  q1 = t2.getTransformedPoint(q1);
  q1 = t3.getTransformedPoint(q1);
  tc = t1.followedBy(t2).followedBy(t3);
  q2 = tc.getTransformedPoint(p1);
  testResult &= ( q1.isCloseTo(q2, 0.000001) );  // fail with gcc

  // test inversion:
  t1 = Trafo(1,2, 3,4, 5,6);
  t2 = t1.inverted();
  p1 = Point(2,3);
  q1 = t1.getTransformedPoint(p1);
  q2 = t2.getTransformedPoint(q1);
  testResult &= ( p1 == q2 );
  t3 = t2.followedBy(t1);
  q3 = t3.getTransformedPoint(p1);
  testResult &= ( p1 == q3 );

  // test construction from 3 points and their images - we choose to construct a reflection about the line y=x-1:
  p1 = Point(-2,1);
  q1 = Point(2,-3);
  p2 = Point(0,1);
  q2 = Point(2,-1);
  p3 = Point(2,2);
  q3 = Point(3,1);
  t1 = Trafo::from3PointsAndImages(p1, p2, p3, q1, q2, q3);
  t2 = Trafo::reflection(1, -1);
  testResult &= ( t1 == t2 );


  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testPoint2D()
{
  std::string testName = "rsPoint2D";
  std::string dummy;
  bool testResult = true;

  testResult &= testPoint2DOperators(dummy);
  testResult &= testPoint2DTransformations(dummy);

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
