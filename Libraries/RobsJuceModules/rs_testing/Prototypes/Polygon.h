#pragma once

// Algorithms for polygons (clipping, convex hull, triangulation, etc.)

// maybe get rid:
typedef RAPT::rsVector2D<float> rsVector2DF;  
//typedef RAPT::rsConicSection<float> rsConicSectionF;
//typedef RAPT::rsEllipse<float> rsEllipseF;
typedef RAPT::rsLine2D<float> rsLine2DF;


//-------------------------------------------------------------------------------------------------
// Utilities:

/** Returns a value proportional to how much the given point p is to the left of the directed edge
from a to b. The proportionlatity constant is the (signed) length of b-a. A point p to the left of 
the edge will return positive values, a point p to the right, negative values and a point on the 
edge zero. */
float edgeFunction(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& p);
// rename to leftOf (or rightTo) a->b is a directed edge or leftDistance

bool isInsideEdge(const rsVector2DF& p, const rsVector2DF& e0, const rsVector2DF& e1);

/** Intersection point between lines through p0,p1 and q0,q1 */
rsVector2DF lineIntersection(const rsVector2DF& p0, const rsVector2DF& p1,
  const rsVector2DF& q0, const rsVector2DF& q1);

/** Computes the signed area of the given polygon. The sign is positive for counterclockwise 
winding order of the vertices and negative otherwise (verify). */
float polygonSignedArea(const std::vector<rsVector2DF>& p);

// todo: bool isClockwise, centerOfMass

//-------------------------------------------------------------------------------------------------
// Clipping:

/** Clips a subject polygon against a convex clipping polygon. The subject needs to to be 
convex. */
std::vector<rsVector2DF> clipPolygon(const std::vector<rsVector2DF>& subject, 
  const std::vector<rsVector2DF>& clipper);
// what about winding? does it matter - probably not as long as its consistent for the two 
// polygons


/** Clips the given polygon p against the edge from e0 to e1. */
std::vector<rsVector2DF> clipAgainstEdge(const std::vector<rsVector2DF>& p,
  const rsVector2DF& e0, const rsVector2DF& e1);
// rename to clipPolygonAgainstedge

std::vector<rsVector2DF> clipConvexPolygons2(const std::vector<rsVector2DF>& p, 
  const std::vector<rsVector2DF>& c);
// obsolete



// new experimental pixel coverage computations - they are based on clipping a triangle to a 
// unit square, so they algorithmically belong to polygon clipping - but application-wise, they 
// belong to drawing/rasterization:
void unitSquareIntersections(const rsVector2DF& p, const rsVector2DF& q, 
  float& x0, float& x1, float& y0, float& y1); 
float unitSquareCut(const rsVector2DF& p, const rsVector2DF& q, 
  float& x0, float& x1, float& y0, float& y1, bool& quadCut); 
float pixelCoverage2(float x, float y, rsVector2DF a, rsVector2DF b, rsVector2DF c);

