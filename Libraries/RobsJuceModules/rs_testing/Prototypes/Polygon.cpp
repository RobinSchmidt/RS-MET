typedef rsVector2DF Vec2;
typedef std::vector<Vec2> ArrVec2;

//-------------------------------------------------------------------------------------------------
// Utilities

float edgeFunction(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& p) 
{
  return Vec2::crossProduct(b-a, p-a);
}

bool isInsideEdge(const rsVector2DF& p, const rsVector2DF& e0, const rsVector2DF& e1)
{
  return edgeFunction(e0, e1, p) <= 0.f;
}
// rename to leftOf (or rightTo)

rsVector2DF lineIntersection(const rsVector2DF& p0, const rsVector2DF& p1,
  const rsVector2DF& q0, const rsVector2DF& q1)
{
  // coeffs for the two implicit line equations:
  float a, b, c;
  rsLine2DF::twoPointToImplicit(p0.x, p0.y, p1.x, p1.y, &a, &b, &c, false);
  float A, B, C;
  rsLine2DF::twoPointToImplicit(q0.x, q0.y, q1.x, q1.y, &A, &B, &C, false);

  // solve 2x2 linear system M*v = r:
  // a*x + b*y = c
  // A*x + B*y = C
  // for x,y
  float M[2][2] = { {a,b}, {A,B} };
  float r[2]    = {  -c,    -C   };
  float v[2];
  rsLinearAlgebra::rsSolveLinearSystem2x2(M, v, r);
  return rsVector2DF(v[0], v[1]);
}
// move to rsLine2D
// actually, it computes an intersection point of the infinitely extended lines...hmm

float polygonSignedArea(const ArrVec2& p)
{
  if(p.size() < 3)
    return 0.f;
  float sum = Vec2::crossProduct(rsLast(p), p[0]);
  for(size_t i = 0; i < p.size()-1; i++)
    sum += Vec2::crossProduct(p[i], p[i+1]);
  return 0.5f * sum;
}

//-------------------------------------------------------------------------------------------------
// Clipping:

// internal function of Sutherland-Hodgman polygon clipper that clips the input polygon against a 
// given edge from e0 to e1, thereby adding zero, one or two vertices to the output polygon
void clipAgainstEdge(const std::vector<rsVector2DF>& in, std::vector<rsVector2DF>& out,
  const rsVector2DF& e0, const rsVector2DF& e1)
{
  //rsVector2DF s, p, i;  // don't use i for a vertex
  //s = in[in.size()-1];
  //for(size_t j = 0; j < in.size(); j++) {
  //  p = in[j];
  //  if(isInsideEdge(p, e0, e1))  { // cases 1,4    
  //    if(isInsideEdge(s, e0, e1))  // case 1 - add polygon vertex
  //      out.push_back(p);
  //    else {                       // case 4 - add intersection and polygon vertex
  //      i = lineIntersection(s, p, e0, e1);
  //      out.push_back(i);
  //      out.push_back(p);
  //    }
  //  }
  //  else  { // cases 2,3
  //    if(isInsideEdge(s, e0, e1)) { // case 2 - add intersection vertex
  //      i = lineIntersection(s, p, e0, e1);
  //      out.push_back(i);
  //    }
  //    else                         // case 3 - add no vertex
  //    {
  //      int dummy = 0;
  //    }
  //  }
  //  s = p;
  //}

  rsVector2DF s = in[in.size()-1];
  for(size_t j = 0; j < in.size(); j++) {
    rsVector2DF p = in[j];
    if(isInsideEdge(p, e0, e1))  {    // cases 1,4    
      if(!isInsideEdge(s, e0, e1)) 
        out.push_back(lineIntersection(s, p, e0, e1));
      out.push_back(p);
    }
    else  {                        // cases 2,3
      if(isInsideEdge(s, e0, e1))  // case 2 - add intersection vertex
        out.push_back(lineIntersection(s, p, e0, e1));
    }
    s = p;
  }
}
// optimization: both isInsideEdge and intersection need to compute the implicit line equation
// coeffs - can be done once (in production code)
// p: general polygon to be clipped, c: convex clipping polygon
std::vector<rsVector2DF> clipConvexPolygons2(const std::vector<rsVector2DF>& p, 
  const std::vector<rsVector2DF>& c)
{
  std::vector<rsVector2DF> r; // result
  clipAgainstEdge(p, r, c[c.size()-1], c[0]);
  for(size_t i = 0; i < c.size()-1; i++)
  {
    clipAgainstEdge(p, r, c[i], c[i+1]);
    //int dummy = 0;
  }
  return r;
}
// Sutherland-Hodgman algorithm (Foley, page 124ff)
// can p really be non-convex? in this case the output may have to be a set of polygons (one 
// non-convex could split into many polygons), see page 125 - or will the algorithm then 
// produce degenerate edges (page 929?)
// this version doesn't work


std::vector<rsVector2DF> clipAgainstEdge(const std::vector<rsVector2DF>& p,
  const rsVector2DF& e0, const rsVector2DF& e1)
{
  std::vector<rsVector2DF> r;
  if(p.size() == 0)
    return r;
  Vec2 S = p[p.size()-1];               // start of edge under consideration
  for(int i = 0; i < p.size(); i++) {   // loop over edges of polynomial
    Vec2 E = p[i];                      // end of edge under consideration
    if(isInsideEdge(E, e0, e1)) {
      if(!isInsideEdge(S, e0, e1))
        r.push_back(lineIntersection(S, E, e0, e1));
      r.push_back(E);
    }
    else {
      if(isInsideEdge(S, e0, e1))
        r.push_back(lineIntersection(S, E, e0, e1));
    }
    S = E;
  }
  return r;
}

std::vector<rsVector2DF> clipPolygon(const std::vector<rsVector2DF>& p, 
  const std::vector<rsVector2DF>& c)
{
  std::vector<rsVector2DF> r = p;
  Vec2 e0 = rsLast(c), e1;
  for(int i = 0; i < c.size(); i++)
  {
    e1 = c[i];
    r  = clipAgainstEdge(r, e0, e1);
    e0 = e1;
  }
  return r;
}
// production version should just be named clip and be a static function of class Polygon

/*
std::vector<rsVector2DF> clipConvexPolygons2(const std::vector<rsVector2DF>& p, 
const std::vector<rsVector2DF>& c)
{
ArrVec2 out = p;          // output polygon
Vec2 e0 = c[c.size()-1];  // current clip edge start
Vec2 e1 = c[0];           // current clip edge end
for(int i = 0; i < c.size(); i++)   {   // loop over clip polygon edges
ArrVec2 in = out;                     // holds the partially clipped polygon
if(in.size() == 0) 
continue;
out.clear();
Vec2 S = in[in.size()-1];             // start point of current edge in subject polygon
for(int j = 0; j < in.size(); j++) {  // loop over vertices in current input polygon
Vec2 E = in[j];                     // end point of current edge in subject polygon
if(isInsideEdge(E, e0, e1)) {
if(!isInsideEdge(S, e0, e1))
out.push_back(lineIntersection(S, E, e0, e1));  // add intersection vertex
out.push_back(E);                                 // add end vertex
}
else if(!isInsideEdge(S, e0, e1))
out.push_back(lineIntersection(S, E, e0, e1));
S = E;
}
e0 = e1;
e1 = c[i];
}
return out; // nope - this doesn't work yet
}
*/
// https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
// List outputList = subjectPolygon;   
// for (Edge clipEdge in clipPolygon) do
//    List inputList = outputList;
//    outputList.clear();
//    Point S = inputList.last;
//    for (Point E in inputList) do
//       if (E inside clipEdge) then
//          if (S not inside clipEdge) then
//             outputList.add(ComputeIntersection(S,E,clipEdge));
//          end if
//          outputList.add(E);
//       else if (S inside clipEdge) then
//          outputList.add(ComputeIntersection(S,E,clipEdge));
//       end if
//       S = E;
//    done
// done

// https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#C
// maybe use this as model for the production version - there's also a function to compute the
// windong of a polynomial (may be useful to use in assertions when counter/clockwise winding
// should be asserted)




// Experimental clipping:
void unitSquareIntersections(const Vec2& p, const Vec2& q, 
  float& x0, float& x1, float& y0, float& y1) // rename to bottom, top, left right
{
  // p, q and stand for a, b or c
  float t;
  t  = -p.y    /   (q.y - p.y);
  x0 =  p.x    + t*(q.x - p.x);
  t  = (1-p.y) /   (q.y - p.y);
  x1 =  p.x    + t*(q.x - p.x);
  t  = -p.x    /   (q.x - p.x);
  y0 =  p.y    + t*(q.y - p.y);
  t  = (1-p.x) /   (q.x - p.x);
  y1 =  p.y    + t*(q.y - p.y);
  // computations can be optimized (do this when the coverage function is finished and a unit test
  // is in place) - but maybe the way it currently looks is ideal for simd processing
  // get rid of pq prefix
}
float unitSquareCut(const Vec2& p, const Vec2& q, 
  float& x0, float& x1, float& y0, float& y1, bool& quadCut)
{
  // get rid of the pq prefix in the x0,x1,...names
  // maybe rename to unitSquareCutArea

  // maybe we need <=, >=

  //bool L = y0 > 0 && y0 < 1;  // left unit square edge crossed by triangle egde p,q
  //bool R = y1 > 0 && y1 < 1;  // right edge crossed 
  //bool B = x0 > 0 && x0 < 1;  // bottom edge crossed
  //bool T = x1 > 0 && x1 < 1;  // top edge crossed
  bool L = y0 >= 0 && y0 <= 1;  // left unit square edge crossed by triangle egde p,q
  bool R = y1 >= 0 && y1 <= 1;  // right edge crossed 
  bool B = x0 >= 0 && x0 <= 1;  // bottom edge crossed
  bool T = x1 >= 0 && x1 <= 1;  // top edge crossed
  // 4 booleans which are either all false or two are false and the other two are true. Below are
  // the enumerated combinations of true-values using parentheses for one that appeared before in 
  // the list (if it appeared, it did so in reverse order - which doesn't matter):
  // LR, LB, LT;   (BL), BR, BT;   (RL), (RB), RT;   (TL), (TR), (TB)
  // the combinations without parentheses are the non-redundant ones, there are 6 of them

  // it seems to be important to use strictly greater/less, otherwise, it returns zero when the
  // edge goes through a corner point, now in these cases, it may return a negative area (the 
  // absolute value is correct, though) -> more unit tests needed


  if(L) {
    if(R) {                  // LR, square is crossed horizontally
      quadCut = true;       // bite is quadrilateral
      if(p.x > q.x)        
        return 0.5f * ((1-y0) + (1-y1)); // edge cuts from right to left and cuts off top quad
      else 
        return 0.5f * (   y0  +    y1);  // edge cuts off bottom region
    }
    else {                   // square crossed diagonally
      quadCut = false;      // bite is triangular
      if(B) // LB
        return 0.5f * x0 * y0;     // bottom edge crossed - return bottom left triangular area
      else  // LT
        return 0.5f * x1 * (1-y0); // top edge crossed - return top-left triangular area
    }
  }
  // LR, LB, LT done (3 of 6)

  if(B) {
    if(T) {  // BT
      quadCut = true;
      if(p.y > q.y)  
        return 0.5f * (   x0  +    x1 );  // edge downward, cuts of left quad (right?)
      else           
        return 0.5f * ((1-x0) + (1-x1));  // verify
                                          // maybe they have to be exchanged?
    }
    else {   // BR
      quadCut = false;
      return 0.5f * (1-x0) * y1;
    }
  }
  // LR, LB, LT, BT, BR done (5 of 6)


  quadCut = false;

  if(T) // only RT remains
    return 0.5f * (1-x1) * (1-y1);

  // 8 return formulas corresponding to 4 possible triangles (top-left, top-right, bottom-left,
  // bottom-right) and 4 possible quadrilaterals (left, right, top, bottom)

  // todo: write unit testsm make sure, all branches are tested

  // production code does not need to compute L,R,T,B in advance - each can be computed as needed
  // -> saves a bit of logic, actually, the quadBite output variable is not used by outside
  // code...however, maybe keep it (at least in the prototype), maybe it's useful later

  return 0; // all booleans were false - nothing is cut off from the square
}
bool isInsideUnitSquare(Vec2 v)
{
  return v.x > 0 && v.x < 1 && v.y > 0 && v.y < 1;
}

float unitSquareCoverage(Vec2 a, Vec2 b, Vec2 c)
{
  // notation: abx0 denotes the x-coordinate of the intersection of the edge (a,b) with the line
  // in the x-direction for which y=0 (i.e. the x-axis), abx1 is the x coordinate of the 
  // intersection with the horizontal line at y=1, aby0: with vertical line at x=0, aby1: with
  // vertical line at x=1, similar definitions for triangle edges bc and ca
  float abB, abT, abL, abR;
  float bcB, bcT, bcL, bcR;
  float caB, caT, caL, caR;
  unitSquareIntersections(a, b, abB, abT, abL, abR);
  unitSquareIntersections(b, c, bcB, bcT, bcL, bcR);
  unitSquareIntersections(c, a, caB, caT, caL, caR);

  // mayb call them abT, abB for (a,b,top), etc
  // abx0 -> abB, abx1 -> abT, aby0 -> abL, aby1 ->abR

  // compute the areas that are cut off from the unit square by the 3 edges:
  bool  abQuad, bcQuad, caQuad;
  float abCut, bcCut, caCut;
  abCut = unitSquareCut(a, b, abB, abT, abL, abR, abQuad);
  bcCut = unitSquareCut(b, c, bcB, bcT, bcL, bcR, bcQuad);
  caCut = unitSquareCut(c, a, caB, caT, caL, caR, caQuad);

  // figure out which of the cut-areas have parts that were cut off twice:

  // condtions for overlap: 
  // -triangle vertex (for example: a) inside unit square
  //  -at least one but maybe both of the edges ab or ca are cutting off a quad


  float area = 1.f - abCut - bcCut - caCut;
  if(isInsideUnitSquare(a))
  {
    float overlap = 0; // amount of overlap in cut-areas resulting from ab and ca

    // ...

    area += overlap;
  }


  //int dummy = 0;



  // can - instead of computing the bite areas directly - compute the clip polygon from this info?
  // that would be better than just the coverage because we can use it to compute the center
  // of gravity and use that for linear deinterpolation
  //int nv = 0;   // number of vertices in clipped polygon (so far)
  Vec2 v[7];    // array for clipped polygon vertices

  return area; // preliminary
}
float pixelCoverage2(float x, float y, Vec2 a, Vec2 b, Vec2 c)
{
  Vec2 d(x, y);
  return unitSquareCoverage(a-d, b-d, c-d);
}
// Idea:
// -shift the original triangle, such the the pixel coordinates can be assumed to be 0,0 like:
//  a.x -= x; a.y -= y; b.x -= x; ...
// -for the 3 triangle edges (a,b), (b,c), (c,a), compute, how much each edge bites off from the
//  the unit square
// -the coverage is one minus the sum of these 3 values
// -each part that is biten off is either the empty polygon, a triangle or a quadrilateral
// -to compute the "bite", we need to compute the intersections of the respective triangle edge 
//  with the 4 lines that contain the edges of the unit square
// -these lines have simple equations: x = 0, x = 1, y = 0, y = 1, these formulas comform the 
//  implicit line equation format A*x + B*y + C = 0
// -....oooh - no this seems to work only, if a,b,c are outside the pixel - but if one of them is 
//  inside, the bite is not a simple triangle or quad anymore...hmmm
// -but maybe it can be made to work, if we consider additional bites..or wait...actually, in this
//  case, parts of the bites get biten off twice - maybe the areas of these parts can be computed 
//  and added back
// -if a is in the square then the bites ab and ca have an overlap that must be added back
//  -that overlap area may be a triangle or quadrilateral

// here has someone else the exact same problem:
// https://stackoverflow.com/questions/22634006/fast-calculation-of-the-intersection-area-of-a-triangle-and-the-unit-square

// maybe implement Liang/Barsky algorithm (clips a polygon against rectangle, see Foley, page 930)

// Foley, page 693 mentions Catmull's object-precision anitaliasing algorithm
// describein in detail here in the paper: A hidden-surface algorithm with anti-aliasing
// https://www.researchgate.net/publication/234810089_A_hidden-surface_algorithm_with_anti-aliasing

// general polygon clipping:
// https://en.wikipedia.org/wiki/Vatti_clipping_algorithm

// triangle/pixel coverage computation:
// http://www.cs.cmu.edu/afs/cs/academic/class/15462-s16/www/lec_slides/2.pdf

// book: Practical Algorithms for 3D Computer Graphics, Second Edition
// https://books.google.de/books?id=NKONAgAAQBAJ