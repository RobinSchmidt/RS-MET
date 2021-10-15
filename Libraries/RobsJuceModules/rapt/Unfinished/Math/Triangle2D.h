#ifndef RAPT_TRIANGLE2D_H
#define RAPT_TRIANGLE2D_H

/** This is a class for representing 2 dimensional triangles in the plane.  */

template<class T>
class rsTriangle2D
{

public:

  //friend class rsAffineTransform2D<T>;  // why?

  /** \name Construction/Destruction */

  /** Constructor. Initializes the 3 points. */
  rsTriangle2D(const rsVector2D<T> &a, const rsVector2D<T> &b, const rsVector2D<T> &c)
  {
    this->a = a;
    this->b = b;
    this->c = c;
  }


  /** \name Inquiry */

  /** Returns true when this triangle contains the given point, false otherwise. */
  bool containsPoint(const rsVector2D<T> &p) const
  {
    return isPointInsideTriangle(p, *this);
  }


  /** \name Operators */

  /** Compares two triangles for equality. */
  bool operator==(const rsTriangle2D& t) const
  {
    if(a == t.a && b == t.b && c == t.c)
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
  static bool arePointsOnSameSideOfLine(const rsVector2D<T> &p1, const rsVector2D<T> &p2,
    const rsVector2D<T> &a, const rsVector2D<T> &b)
  {
    rsVector2D<T> aa = a, bb = b, pp1 = p1, pp2 = p2; 
    // needed because rsVector2D::operator- is not const -> fix this

    T cp1 = rsVector2D<T>::crossProduct(bb-aa, pp1-aa);
    T cp2 = rsVector2D<T>::crossProduct(bb-aa, pp2-aa);
    //if(rsVector2D<T>::scalarProduct(cp1, cp2) >= 0)  // old - wrong? doesn't seem to make sense
    if(cp1*cp2 >= 0)  // new - makes more sense but is a gues -> verify
      return true;
    else
      return false;
  }

  /** Determines whether the given point is inside the given triangle or not. */
  static bool isPointInsideTriangle(const rsVector2D<T> &p, const rsTriangle2D<T> &t)
  {
    if(arePointsOnSameSideOfLine(p, t.a, t.b, t.c)
      && arePointsOnSameSideOfLine(p, t.b, t.a, t.c)
      && arePointsOnSameSideOfLine(p, t.c, t.a, t.b))
    {
      return true;
    }
    else
      return false;
  }
  // can be simplified using the edge function,
  // see https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage
  // when doing so, check, if the unit test still passes

protected:

  rsVector2D<T> a, b, c; // vertices - maybe make public

};

// todo: getNormal (returns normal vector)

#endif
