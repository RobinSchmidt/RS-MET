#pragma once

// Algorithms for polygons (clipping, convex hull, triangulation, etc.)


// Utilities:

float edgeFunction(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& p);

bool isInsideEdge(const rsVector2DF& p, const rsVector2DF& e0, const rsVector2DF& e1);

rsVector2DF lineIntersection(const rsVector2DF& p0, const rsVector2DF& p1,
  const rsVector2DF& q0, const rsVector2DF& q1);











/** Returns a value proportional to how much the given point p is to the left of the directed edge
from a to b. The proportionlatity constant is the (signed) length of b-a. A point p to the left of 
the edge will return positive values, a point p to the right, negative values and a point on the 
edge zero. */
float edgeFunction(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& p);
// rename to leftOf (or rightTo) a->b is a directed edge or leftDistance

/** Intersection point between lines through p0,p1 and q0,q1 */
rsVector2DF lineIntersection(const rsVector2DF& p0, const rsVector2DF& p1,
  const rsVector2DF& q0, const rsVector2DF& q1);

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

/** Computes the area of the given polygon. */
float polygonArea(const std::vector<rsVector2DF>& p);

/** Given pixel coordinates x,y and the 3 vertices a,b,c of a triangle (in counterclockwise order), 
this function computes the area in which the triangle and pixel-square overlap. It's a value 
between 0 and 1, where 0 means the pixel and triangle do not intersect at all, 1 means the pixel is
fully covered by the triangle. */
float pixelCoverage(int x, int y, 
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c);


// new pixel coverage computations:
void unitSquareIntersections(const rsVector2DF& p, const rsVector2DF& q, 
  float& x0, float& x1, float& y0, float& y1); 
float unitSquareCut(const rsVector2DF& p, const rsVector2DF& q, 
  float& x0, float& x1, float& y0, float& y1, bool& quadCut); 
float pixelCoverage2(float x, float y, rsVector2DF a, rsVector2DF b, rsVector2DF c);
