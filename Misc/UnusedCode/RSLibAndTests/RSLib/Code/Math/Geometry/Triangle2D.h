#ifndef RS_TRIANGLE2D_H
#define RS_TRIANGLE2D_H

namespace RSLib
{

  /**

  This is a class for representing 2 dimensional triangles in the plane.

  */

  template<class RealType>
  class rsTriangle2D
  {

  public:

    //friend class rsAffineTransform2D<RealType>;  // why?

    /** \name Construction/Destruction */

    /** Constructor. Initializes the 3 points. */
    rsTriangle2D(const rsPoint2D<RealType> &a, 
                 const rsPoint2D<RealType> &b, 
                 const rsPoint2D<RealType> &c)
    {
      this->a = a;
      this->b = b;
      this->c = c;
    }


    /** \name Inquiry */

    /** Returns true when this triangle contains the given point, false otherwise. */
    bool containsPoint(const rsPoint2D<RealType> &p) const
    {
      return isPointInsideTriangle(p, *this);
    }


    /** \name Operators */

    /** Compares two triangles for equality. */
    bool operator==(const rsTriangle2D& t) const
    {
      if( a == t.a && b == t.b && c == t.c )
        return true;
      else
        return false;
    }

    /** Compares two triangles for inequality. */
    bool operator!=(const rsTriangle2D& t) const
    {
      return !(*this == t);
    }


    /** \name Static Member Functions */

    /** Determines whether the two points p1 and p2 are on the same side of the line that connects 
    a and b.
    \todo move this function into a class 'Line' */
    static bool arePointsOnSameSideOfLine(const rsPoint2D<RealType> &p1, 
                                          const rsPoint2D<RealType> &p2, 
                                          const rsPoint2D<RealType> &a,
                                          const rsPoint2D<RealType> &b)
    {
      rsPoint2D<RealType> aa = a, bb = b, pp1 = p1, pp2 = p2; // needed because Point::operator- 
                                                              // is not const
      RealType cp1 = rsPoint2D<RealType>::crossProduct(bb-aa, pp1-aa);
      RealType cp2 = rsPoint2D<RealType>::crossProduct(bb-aa, pp2-aa);
      if( rsPoint2D<RealType>::scalarProduct(cp1, cp2) >= 0 )
        return true;
      else
        return false;
    }

    /** Determines whether the given point is inside the given triangle or not. */
    static bool isPointInsideTriangle(const rsPoint2D<RealType> &p, 
                                      const rsTriangle2D<RealType> &t)
    {
      if(  arePointsOnSameSideOfLine(p, t.a, t.b, t.c)
        && arePointsOnSameSideOfLine(p, t.b, t.a, t.c)
        && arePointsOnSameSideOfLine(p, t.c, t.a, t.b) )
      {
        return true;
      }
      else
        return false;
    }

  protected:

    rsPoint2D<RealType> a, b, c;

  };

  // a typedef'd explicit instantiation for coordinates of type double:
  typedef rsTriangle2D<double> rsDblTriangle2D; // needs the RSLib_API prefix?

}

#endif
