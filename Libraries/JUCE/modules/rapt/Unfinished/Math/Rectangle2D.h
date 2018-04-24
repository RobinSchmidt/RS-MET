#ifndef RAPT_RECTANGLE2D_H
#define RAPT_RECTANGLE2D_H

  /**

  This is a class for representing 2 dimensional rectangles that are aligned with the coordinate
  axes.

  */

template<class CoordinateType>
class rsRectangle2D
{

public:


  /** \name Construction/Destruction */

  /** Constructor. */
  rsRectangle2D(CoordinateType xOrigin = 0, CoordinateType yOrigin = 0, CoordinateType width = 0, CoordinateType height = 0) :
    x(xOrigin), y(yOrigin), w(width), h(height)
  {

  }


  /** \name Setup */

  /** Copies the data members from another rectangle into this. */
  void setFromOtherRectangle(const rsRectangle2D& r)
  {
    x = r.x;
    y = r.y;
    w = r.w;
    h = r.h;
  }


  /** \name Inquiry */

  /** Returns the x-coordinate of the rectangle. */
  inline CoordinateType getX() const
  {
    return x;
  }

  /** Returns the y-coordinate of the rectangle. */
  inline CoordinateType getY() const
  {
    return y;
  }

  /** Returns the width of the rectangle. */
  inline CoordinateType getWidth() const
  {
    return w;
  }

  /** Returns the height of the rectangle. */
  inline CoordinateType getHeight() const
  {
    return h;
  }

  /** Returns the maximum x-coordinate of this rectangle. In most cases, this will be interpreted
  as the right edge. */
  inline int getMaxX() const
  {
    return x+w;
  }

  /** Returns the maximum y-coordinate of this rectangle. For upward y-axes (such as mostly in
  mathematics), this should be interpreted as top edge, for downward y-axes (such as mostly when
  drawing on screen), this should be interpreted as bottom edge. */
  inline int getMaxY() const
  {
    return y+h;
  }

  /** Returns the origin of the rectangle. */
  inline rsPoint2D<CoordinateType> getOrigin() const
  {
    return rsPoint2D<CoordinateType>(x, y);
  }

  /** Returns true when this rectangle contains the given point, false otherwise. Points on the
  border are considered inside. */
  bool containsPoint(CoordinateType x, CoordinateType y) const
  {
    if(x >= this->x && y >= this->y && x <= (this->x + w) && y <= (this->y + h))
      return true;
    else
      return false;
  }


  /** \name Operators */

  /** Compares two rectangles for equality. */
  bool operator==(const rsRectangle2D& r) const
  {
    if(x == r.x && y == r.y && w == r.w && h == r.h)
      return true;
    else
      return false;
  }

  /** Compares two triangles for inequality. */
  bool operator!=(const rsRectangle2D& r) const
  {
    return !(*this == r);
  }


  /** \name Misc */


  /** Returns the intersection between rectangles r1 and r2. If the intersection is empty, it
  will return the default rectangle with x = y = w = h = 0. */
  static rsRectangle2D<CoordinateType> intersect(const rsRectangle2D<CoordinateType> &r1,
    const rsRectangle2D<CoordinateType> &r2)
  {
    CoordinateType xMin, xMax, yMin, yMax, w, h;

    xMin = rsMax(r1.x, r2.x);
    yMin = rsMax(r1.y, r2.y);
    xMax = rsMin(r1.getMaxX(), r2.getMaxX());
    yMax = rsMin(r1.getMaxY(), r2.getMaxY());
    w    = xMax-xMin;
    h    = yMax-yMin;

    if(w <= 0 || h <= 0)
      return rsRectangle2D<CoordinateType>(0, 0, 0, 0);
    else
      return rsRectangle2D<CoordinateType>(xMin, yMin, w, h);
  }

  // \todo make a function for a union-rectangle (the smallest rectangle that contains both of 
  // the constituting rectangles to be united)


  /** \name Data */

  CoordinateType x, y, w, h;
    // \todo maybe use an rsPoint2D as "origin" and an rsSize2D as "size"

};

// \todo provide explicit intantiation for int

#endif
