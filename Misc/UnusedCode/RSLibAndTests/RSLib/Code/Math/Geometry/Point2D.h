#ifndef RS_POINT2D_H
#define RS_POINT2D_H

namespace RSLib
{

  template<class CoordinateType>
  class AffineTransform2D;
    // forward declaration, needed to make AffineTransform2D<CoordinateType> a friend of
    // Point2D<CoordinateType>

  /**

  This is a class for representing 2 dimensional points in the plane.

  \todo: rename into Vector2D - this is more appropriate because formally, it makes no sense to
         multiply numbers with points but with vectors, it does
  \todo: implement perp-dot product
  \todo: implement stereographic mapping to and from the Riemann sphere
         (see Visual Complex Analysis, Ch.2)

  */

  template<class CoordinateType>
  class rsPoint2D
  {

  public:

    //friend class AffineTransform2D<CoordinateType>;

    /** \name Construction/Destruction */

    /** Constructor. Initializes x and y with the passed parameters. */
    rsPoint2D(CoordinateType x = CoordinateType(0), CoordinateType y = CoordinateType(0))
    {
      this->x = x;
      this->y = y;
    }


    /** \name Inquiry */

    /** Returns the length of the vector from the origin to the point.
    \todo maybe rename to getCoordinateVectorLength - or rename the whole class into Vector2D -
    length dos not make sense for points
    */
    CoordinateType getLength() const
    {
      return rsSqrt(x*x + y*y); // try hypot instead
    }

    /** Returns the Euclidean distance of this point to the point in the argument. */
    CoordinateType distanceTo(const rsPoint2D<CoordinateType> &p) const
    {
      rsPoint2D<CoordinateType> deltaVector = *this-p;
      return deltaVector.getLength();
    }

    /** Returns true, when this point is inside a disc of radius 'margin' centered around the point
    in the argument. */
    bool isCloseTo(const rsPoint2D<CoordinateType> &p, CoordinateType margin) const
    {
      return distanceTo(p) <= margin;
    }

    /** \name Operators */

    /** Compares two points for equality. */
    bool operator==(const rsPoint2D& p) const
    {
      if( x == p.x && y == p.y )
        return true;
      else
        return false;
    }

    /** Compares two points for inequality. */
    bool operator!=(const rsPoint2D& p) const
    {
      if( x != p.x || y != p.y )
        return true;
      else
        return false;
    }

    /** Defines the negative of a point. */
    rsPoint2D operator-() const
    {
      return rsPoint2D(-x, -y);
    }

    /** Adds another point to this point and returns the result. */
    rsPoint2D& operator+=(const rsPoint2D &p)
    {
      x += p.x;
      y += p.y;
      return *this;
    }

    /** Subtracts another point from this one and returns the result. */
    rsPoint2D& operator-=(const rsPoint2D &p)
    {
      x -= p.x;
      y -= p.y;
      return *this;
    }

    /** Multiplies this point by a real number and returns the result. */
    rsPoint2D& operator*=(const CoordinateType &r)
    {
      x *= r; y *= r; return *this;
    }

    /** Divides this point by a real number and returns the result. */
    rsPoint2D& operator/=(const CoordinateType &r)
    {
      CoordinateType s = CoordinateType(1) / r;
      x *= s;
      y *= s;
      return *this;
    }

    /** Adds two points. */
    rsPoint2D<CoordinateType> operator+(const rsPoint2D<CoordinateType> &p) const
    {
      return rsPoint2D<CoordinateType>(x+p.x, y+p.y);
    }

    /** Subtracts two points. */
    rsPoint2D<CoordinateType> operator-(const rsPoint2D<CoordinateType> &p) const
    {
      return rsPoint2D<CoordinateType>(x-p.x, y-p.y);
    }

    /** Multiplies a point by a real number. */
    rsPoint2D<CoordinateType> operator*(const CoordinateType &r) const
    {
      return rsPoint2D<CoordinateType>(x*r, y*r);
    }

    /** Divides a point by a real number. */
    rsPoint2D<CoordinateType> operator/(const CoordinateType &r) const
    {
      CoordinateType s = CoordinateType(1) / r;
      return rsPoint2D<CoordinateType>(x*s, y*s);
    }


    /** \name Static Member Functions */

    /** Returns the cross-product of two points defined as x1*y2-x2*y1. The cross-product can be
    interpreted as the signed area of the parallelogram formed by the points (0,0),p1,p2,p1+p2
    where the sign is positive when p1 is clockwise from p2 with respect to the origin, negative
    when it's counterclockwise and zero when they are collinear (pointing in the same or opposite
    directions).  */
    static CoordinateType crossProduct(const rsPoint2D<CoordinateType> &p1,
                                       const rsPoint2D<CoordinateType> &p2)
    {
      return p1.x*p2.y - p2.x*p1.y;
    }

    /** Returns the scalar-product (aka inner product or dot product) of two points defined as
    x1*x2+y1*y2. When the two vectors are of unit length, it can be interpreted as the cosine of
    the angle between the two vectors. Because the cosine of a 90 degree angle is zero, this
    implies that a zero scalar-product indicates orthogonal vectors. If only p1=(x1,y1) is a unit
    vector, then the scalar-product can be interpreted as the length of the orthogonal projection
    of p2 onto p1. */
    static CoordinateType scalarProduct(const rsPoint2D<CoordinateType> &p1,
                                        const rsPoint2D<CoordinateType> &p2)
    {
      return p1.x*p2.x + p1.y*p2.y;
    }


    /** \name Data */

    CoordinateType x, y;

  };

  /** Multiplies a real number and a point. We need to define this operator outside the class
  because the left operand is not of class Point2D (but of some real number type). We use the fact
  that number*point == point*number to compute the result.  */
  template<class CoordinateType>
  inline rsPoint2D<CoordinateType> operator*(const CoordinateType &r,
                                             const rsPoint2D<CoordinateType> &p)
  {
    rsPoint2D<CoordinateType> tmp(p);
    tmp *= r;
    return tmp;
  }

}

#endif
