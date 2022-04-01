#pragma once

// rename to Rasterization

// maybe get rid:
typedef RAPT::rsImage<float> rsImageF;
typedef RAPT::rsAlphaMask<float> rsAlphaMaskF;
typedef RAPT::rsImagePainter<float, float, float> rsImagePainterFFF;
typedef RAPT::rsImageDrawer<float, float, float> rsImageDrawerFFF;
typedef RAPT::rsLineDrawer<float, float, float> rsLineDrawerFFF;
typedef RAPT::rsPhaseScopeBuffer<float, float, double> rsPhaseScopeBufferFFD;

//-------------------------------------------------------------------------------------------------
// Utilities

/** Given pixel coordinates x,y and the 3 vertices a,b,c of a triangle (in counterclockwise order), 
this function computes the area in which the triangle and pixel-square overlap. It's a value 
between 0 and 1, where 0 means the pixel and triangle do not intersect at all, 1 means the pixel is
fully covered by the triangle. */
float pixelCoverage(int x, int y, 
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c);

//-------------------------------------------------------------------------------------------------
// Lines

void drawLineWuPrototype(rsImageF& img, float x0, float y0, float x1, float y1, float color);
void drawLineBresenham(rsImageF& img, int x0, int y0, int x1, int y1, float color);

void drawThickLine(rsImageF& img, float x0, float y0, float x1, float y1, float color,
  float thickness, bool roundCaps = false);
void drawThickLine2(rsImageF& img, float x0, float y0, float x1, float y1, float color,
  float thickness, int endCaps = 0);
void drawThickLineViaWu(rsImageF& img, float x0, float y0, float x1, float y1, float color, 
  float thickness);

void plotLineWidth(rsImageF& img, int x0, int y0, int x1, int y1, float wd); //?

// line profile functions
float lineIntensity1(float d, float t2);
float lineIntensity2(float d, float t2);
float lineIntensity3(float d, float t2);
float lineIntensity4(float d, float t2);


//-------------------------------------------------------------------------------------------------
// Triangles

/** Fills all pixels whose centers are inside the given triangle with the given color. If a pixel
center is on an edge, it will be considered inside, if it's a top or a left edge 
("top-left rule"). It corresponds to a convention where a pixel with coordinate i (i being the 
integer x or y coordinate) being defined as covering the range [i, i+1), i.e. an interval that 
is closed to the left and open to ther right, i.e. including i but excluding i+1.
The order of the vertices doesn't matter.
A top edge, is an edge that is exactly horizontal and is above the other edges.
A left edge, is an edge that is not exactly horizontal and is on the left side of the triangle. 
A triangle can have one or two left edges.
see:
https://docs.microsoft.com/en-us/windows/desktop/direct3d11/d3d10-graphics-programming-guide-rasterizer-stage-rules#triangle-rasterization-rules-without-multisampling
If this algorithms is used to render several triangles with shared edges (as in a triangle mesh),
each pixel will be visited only once. */
void drawTriangle(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color);

void drawTriangleAntiAliasedProto(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color);
// todo: make sure correct winding order - either make an assertion or re-order vertices here
// uses unweighted area sampling -> disadvantage: objects that are smaller than one pixel don't 
// seem to move at all when they move within one pixel such that no part ever sticks out
// for that, weighted area sampling (maybe with a conic filter) can be used - but exact weighted
// coverage computations are complicated ...maybe use a 3x3 pixel filter with piecewise constant
// weighting, like
// |1 2 1|
// |2 4 2| * 1/12,  ..maybe let the user set the filter kernel coeffs
// |1 2 1|
// default:
// |0 0 0|
// |0 1 0|
// |0 0 0|
// or maybe a 2x2 matrix:
// |1 1| * 1/4
// |1 1| 
// or maybe find the center of gravity of the polygon (triangle-clipped-to-square) and use the
// "paint-dot" operation (with bilinear deinterpolation)


// fast version:
// dispatcher between box-based and span-based:

/** An optimized version - it dispatches between drawTriangleAntiAliasedBoxBased and 
drawTriangleAntiAliasedSpanBased based on the size of the triangle (for small triangles, it may
not be worth it do to the span-based optimizations - the overhead pays off only for larger
triangles). */
void drawTriangleAntiAliased(rsImageDrawerFFF& drw, 
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color);

/** As a simple and obvious optimization, this function loops only over those pixels that are 
contained in the bounding box of the triangle. */
void drawTriangleAntiAliasedBoxBased(rsImageDrawerFFF& drw, 
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color);

/** As a further optimization, this functions computes 3 spans for each scanline where the first 
and last spans compute coverages and the middle section rund over the pixles that are fully 
covered, so no coverage needs to be computed (it's just 1 everywhere there). */
void drawTriangleAntiAliasedSpanBased(rsImageDrawerFFF& drw, 
  const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, float color);
// this does not work yet!!!


/** Class for representing an axis aligned rectangle. Mainly for simplifying computations for 
clipping aginst windows or pixels. */

class RectangleF  // just Rectangle doesn't compile (name conflict?)
{

public:

  float xMin = 0.f, xMax = 1.f, yMin = 0.f, yMax = 1.f;

  enum regions
  {
    INSIDE = 0,
    TOP_LEFT,
    CENTER_LEFT,
    BOTTOM_LEFT,
    BOTTOM_CENTER,
    BOTTOM_RIGHT,
    CENTER_RIGHT,
    TOP_RIGHT,
    TOP_CENTER
  };

  /** Returns in which of the region the given point falls. */
  int getRegion(const rsVector2DF& p);

  /** Returns whether or not the given point is inside the unit square. Points on the boundary are
  considered inside. */
  bool isInside(const rsVector2DF& p) { return getRegion(p) == INSIDE; }

  /** Given a line L defined by the parametric equation L(t) = a + t*(b-a), this function returns 
  the value of the parameter t for which the line L intersects the top edge of the unit square. If 
  the line is parallel, it will be plus or minus infinity. */
  float getIntersectionParameterTop(const rsVector2DF& a, const rsVector2DF& b)
  {
    return (yMax-a.y) / (b.y-a.y);
  }

  float getIntersectionParameterLeft(const rsVector2DF& a, const rsVector2DF& b)
  {
    return (xMin-a.x) / (b.x-a.x);
  }

  float getIntersectionParameterBottom(const rsVector2DF& a, const rsVector2DF& b)
  {
    return (yMin-a.y) / (b.y-a.y);
  }

  float getIntersectionParameterRight(const rsVector2DF& a, const rsVector2DF& b)
  {
    return (xMax-a.x) / (b.x-a.x);
  }

};

/** Given a triangle defined by the vertices a,b,c, this function returns the number of vertices of
the polygon that results from clipping the triangle to the unit square. It writes the polygon 
vertices into the array p which should be of length 7 (this is the maximum number of vertices that 
the resulting polygon can have). If 0 is returned, the triangle doesn't overlap the unit-square at 
all. */
int clipTriangleToUnitSquare(const rsVector2DF& a, const rsVector2DF& b, const rsVector2DF& c, 
  rsVector2DF* p);
// under construction


// todo: make class rsRenderer
// setMethod(int);   // rasterize or raytrace
// setAntiAlias(int); // off, unweighted area sampling, weighted area sampling, bilinear deinterpolation, supersampling
// setAreaSamplingFilter(rsImage* filterKernel); // or maybe pass a function-pointer f(float x, float y) definign the continuous kernel function?
// renderScene(rsScene*, rsImage*);
// renderTriangle(rsVertex* a, rsVertex* b, rsVertex* c);
// setVertexShader(rsVertexShader*);
// if i get serious about than, maybe i should call my company Media engineering tools - acronym 
// stays the same..and maybe add a paintbrush that paints a triangle to the logo :-)

// recommended book
// Physically Based Rendering: From Theory To Implementation
// https://www.amazon.de/Physically-Based-Rendering-Matt-Pharr/dp/0128006455/ref=dp_ob_title_bk
// has a website with source code (BSD licensed)
// http://www.pbrt.org/
// https://github.com/mmp/pbrt-v3

//=================================================================================================

/** Returns true if the pixel at coordinates x,y is an interior pixel of the given image. */
bool isInteriorPixel(int x, int y, const rsImageF& img);


/** Given an image "img", this function sets every entry in the "classes" pixel matrix to the
given "flatClassLabel", iff the corresponding pixel in img belongs to a flat region. The criterion
for a pixel to belong to a flat region is that it must have the same value/color as its 3x3 
neighborhood up to some given tolerance. The other entries of the "classes" matrix, i.e. those 
corresponding to pixels in non-flat regions are left as is. */
void classifyFlatPixels3x3(const rsImageF& img, rsImage<char>& classes, char flatClassLabel, 
  float tol = 0.f);

/** Under construction...

Finds the pixel closest to the given location xc, yc (the "center pixel") in the given image 
that satisfies the given predicate in the search direction given by dx, dy. For example, 
north-north-west would be given by dx = -1, dy = -2. It uses a Bresenham-style update equation for
scanning through the pixels and for each encountered pixel (where Bresenham would color the pixel),
it checks the predicate. If it's true, it returns the pixel's location, if it's false, it moves on
to the next pixel. The predicate should take the color of the pixel at xc, yc as first parameter 
and the color of the pixel along the search line as second parameter and return a bool. ..tbc...*/
template<class P>
rsVector2D<int> findClosestPixelWith(const rsImageF& img, P pred, int xc, int yc, int dx, int dy);


/** Returns an array of coordinates (x,y) for which C(x,y) == c. */
std::vector<rsVector2D<int>> findAll(const rsImage<char>& C, char c);

/** Turns flat regions in an image into smooth gradients. This algorithm is be used, for example, 
to make the typical outputs of fractal rendering algorithms look nicer. They typically assign a 
color to a pixel based on the number of iterations required to ensure divergence which leads to a 
lot of flat color regions liek in stepped color-gradients where we would very much prefer to see
a smooth gradient. This function aims to turn these steppy gradients onto smooth ones but without
blurring ther overall image. In some sense, it attempts to be some sort of opposite of the countour
filling algorithm. It's a bit like blurring but only applied to flat color regions. Well,
that's a simplification...tbc... */
int gradientifyFlatRegions(const rsImageF& in, rsImageF& out, int numPasses = 1);
// I think, it can't work in place: in and out must be distinct. Maybe the result should be passed
// as return value. We may need an internal temporary image with extended borders anyway...


int gradientifyFlatRegions2(const rsImageF& in, rsImageF& out, int numPasses = 1);


//=================================================================================================

struct rsPixelRGB // maybe rename to rsPixel24BitRGB
{
  // have constructors that take triples of floats, doubles, etc. they should optionally clip

  /** Default constructor. Sets r = g = b = 0 (actually, it just leaves them alone - but they are 
  intialized to 0 in the class, so they remain at 0).  */
  rsPixelRGB() {}


  rsPixelRGB(unsigned char gray) 
  { 
    r = g = b = gray;
  }

  rsPixelRGB(int gray) 
  { 
    r = g = b = (unsigned char) (gray);
  }


  /** Creates a pixel from 3 floating point numbers and optionally clips the values to the
  allowed range before converting to 8-bit integers. */
  rsPixelRGB(float R, float G, float B, bool clip = false)
  {
    if(clip) {
      R = rsClip(R, 0.f, 1.f);
      G = rsClip(G, 0.f, 1.f);
      B = rsClip(B, 0.f, 1.f);
    }
    r = (unsigned char)(255.f * R);
    g = (unsigned char)(255.f * G);
    b = (unsigned char)(255.f * B);
  }

  unsigned char r = 0, g = 0, b = 0;
};

inline bool rsAlmostEqual(const unsigned char& x, const unsigned char& y, const unsigned char& tol)
{
  if(x > y) {
    if(x - y > tol)
      return false; }
  else {
    if(y - x > tol)
      return false; }
  return true;
}

// explicit specialization:
inline bool rsAlmostEqual(const rsPixelRGB& x, const rsPixelRGB& y, const rsPixelRGB& tol)
{
  return rsAlmostEqual(x.r, y.r, tol.r)
    &&   rsAlmostEqual(x.g, y.g, tol.g) 
    &&   rsAlmostEqual(x.b, y.b, tol.b);
}

rsImage<rsPixelRGB> rsConvertImage(
  const rsImage<float>& R, const rsImage<float>& G, const rsImage<float>& B, bool clip);


void rsConvertImage(
  const rsImage<float>& R, const rsImage<float>& G, const rsImage<float>& B, bool clip, 
  rsImage<rsPixelRGB>& img);


rsImage<rsPixelRGB> rsConvertImage(const rsImage<rsFloat32x4>& in, bool clip = true);

void rsConvertImage(const rsImage<rsFloat32x4>& in, rsImage<rsPixelRGB>& out, bool clip = true);
// this function has a signature that is not consistent with those above..i think, this here makes
// more sense...update the functions above and all dependent code