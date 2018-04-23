#ifndef RAPT_AFFINETRANSFORM2D_H
#define RAPT_AFFINETRANSFORM2D_H

/**

This is a class for representing affine transformations of 2-dimensional points (of type
Point2D). An affine transformation is a linear transformation followed by a translation.
Mathematically, it can be expressed as:
\f[ \mathbf{y = Ax + b} \f]
with some Matrix \f$ \mathbf{A} \f$ to represent the linear transformation and a vector
\f$ \mathbf{b} \f$ for the translation. The class provides a bunch of static member functions
which can be used to construct primitive transforms. Then, general transfroms can be obtained
by chaining those primitive transforms, for example, a scaling followed by a translation
followed by a rotation for doubles can be constructed as:

\code
typedef AffineTransform2D<double> Trafo;  // defines an alias for the ugly template syntax
Trafo s = Trafo::scaling(2,3);            // scales x by 2 and y by 3
Trafo t = Trafo::translation(1,2);        // translates x by 1 and y by 2
Trafo r = Trafo::rotation(M_PI/4);        // rotates by 45 degrees
Trafo compoundTrafo = s.followedBy(t).followedBy(r);
\endcode

\todo check if the the code example works
\todo implement scaling around a center other than (0,0), chainBefore, rotated, sheared,
      translated, reflected, scaled
\todo setValues, setToIdentity
\todo maybe write a similar class for Moebius transforms of complex numbers (maybe with a
      baseclass ConformalMap)
\todo move the LaTeX comments away from the function declarations - they are cluttering them
too much */

template<class RealType>
class rsAffineTransform2D
{

public:

  /** \name Construction/Destruction */

  /** Constructor. Initializes to the identity transform with
  \f[
  \mathbf{A}
  =
  \begin{pmatrix}
  1 & 0 \\
  0 & 1
  \end{pmatrix} , \qquad
  \mathbf{b}
  =
  \begin{pmatrix}
  0 \\
  0
  \end{pmatrix}
  \f]
  */
  rsAffineTransform2D()
  {
    a11 = a22 = RealType(1);
    a12 = a21 = b1 = b2 = RealType(0);
  }

  /** Constructor. Initializes the transformation matrix with the given values, where:
  \f[
  \mathbf{A}
  =
  \begin{pmatrix}
  a_{11} & a_{12} \\
  a_{21} & a_{22}
  \end{pmatrix} , \qquad
  \mathbf{b}
  =
  \begin{pmatrix}
  b_1 \\
  b_2
  \end{pmatrix}
  \f]
  */
  rsAffineTransform2D(RealType a11, RealType a12, RealType a21, RealType a22,
    RealType b1, RealType b2)
  {
    this->a11 = a11; this->a12 = a12;
    this->a21 = a21; this->a22 = a22;
    this->b1  = b1;  this->b2  = b2;
  }


  /** \name Operators */

  /** Compares two transforms for equality. */
  bool operator==(const rsAffineTransform2D& t) const
  {
    if(a11==t.a11 && a12==t.a12 && a21==t.a21 && a22==t.a22 && b1==t.b1 && b2==t.b2)
      return true;
    else
      return false;
  }

  /** Compares two transformations for inequality. */
  bool operator!=(const rsAffineTransform2D& t) const
  {
    return !(*this == t);
  }


  /** \name Factory Functions */

  /** Returns an identity transform. */
  static const rsAffineTransform2D<RealType> identity()
  {
    return rsAffineTransform2D<RealType>(1, 0, 0, 1, 0, 0);
  }

  /** Returns a transformation that scales the \f$ x \f$ and \f$ y \f$ coordinates of a vector by
  the given scale factors. The transformation is of the form:
  \f[
  \mathbf{A}
  =
  \begin{pmatrix}
  s_x & 0 \\
  0   & s_y
  \end{pmatrix} , \qquad
  \mathbf{b}
  =
  \begin{pmatrix}
  0 \\
  0
  \end{pmatrix}
  \f]
  where \f$ s_x, sy \f$ are the scale factors for \f$ x \f$ and \f$ y \f$ respectively. */
  static const rsAffineTransform2D<RealType> scaling(RealType sx, RealType sy)
  {
    return rsAffineTransform2D<RealType>(sx, 0, 0, sy, 0, 0);
  }

  /** Returns a translation that translates x to x+b1 and y to y+b2. */
  static const rsAffineTransform2D<RealType> translation(RealType b1, RealType b2)
  {
    return rsAffineTransform2D<RealType>(1, 0, 0, 1, b1, b2);
  }

  /** Returns a (counterclockwise) rotation around the origin by the given angle. This
  transformation is of the form:
  \f[
  \mathbf{A}
  =
  \begin{pmatrix}
  \cos \phi & -\sin \phi \\
  \sin \phi &  \cos \phi
  \end{pmatrix} , \qquad
  \mathbf{b}
  =
  \begin{pmatrix}
  0 \\
  0
  \end{pmatrix}
  \f]
  where \f$ \phi \f$ is the angle in radians. */
  static const rsAffineTransform2D<RealType> rotation(RealType angleInRadians)
  {
    RealType c = (RealType)cos(angleInRadians);
    RealType s = (RealType)sin(angleInRadians);  // todo: use sinCos for doubles
    return rsAffineTransform2D<RealType>(c, -s, s, c, 0, 0);
  }

  /** Returns a (counterclockwise) rotation around an arbitrary center by the given angle. */
  static const rsAffineTransform2D<RealType> rotation(RealType angleInRadians,
    rsPoint2D<RealType> center)
  {
    rsAffineTransform2D<RealType> t =
      rsAffineTransform2D<RealType>::translation(-center.x, -center.y);
    t = t.followedBy(rsAffineTransform2D<RealType>::rotation(angleInRadians));
    t = t.followedBy(rsAffineTransform2D<RealType>::translation(center.x, center.y));
    return t;
  }

  /** Returns a (counterclockwise) rotation around the origin by \f$ \pi/2 \f$ (90 degrees) of
  the form
  \f[
  \mathbf{A}
  =
  \begin{pmatrix}
  0 & -1 \\
  1 &  0
  \end{pmatrix} , \qquad
  \mathbf{b}
  =
  \begin{pmatrix}
  0 \\
  0
  \end{pmatrix}
  \f] */
  static const rsAffineTransform2D<RealType> rotationByHalfPi()
  {
    return rsAffineTransform2D<RealType>(0, -1, 1, 0, 0, 0);
  }

  /** \todo: check, how exactly shear is defined. */
  static const rsAffineTransform2D<RealType> shear(RealType shearX, RealType shearY)
  {
    return  rsAffineTransform2D<RealType>(1, shearX, shearY, 1, 0, 0);
  }

  /** Creates a reflection about the line \f$ ax + by + c = 0 \f$. */
  static rsAffineTransform2D<RealType> reflection(RealType a, RealType b, RealType c)
  {
    rsAffineTransform2D<RealType> r;
    RealType tmp1 = -1 / (a*a + b*b);
    RealType tmp2 = (a*a - b*b);
    RealType tmp3 = 2*a*b;
    RealType tmp4 = tmp1*tmp2;
    RealType tmp5 = tmp1*tmp3;
    r.a11 =  tmp4;
    r.a12 =  tmp5;
    r.a21 =  tmp5;
    r.a22 = -tmp4;
    r.b1  =  tmp1 * 2*a*c;
    r.b2  =  tmp1 * 2*b*c;
    return r;
  }

  /** Creates a reflection about the line \f$ y = ax + b \f$. For the more general line equation
  \f$ ax + by + c = 0 \f$ (which can also represent vertical lines) use the three-parametric
  version of this function (but note how the meaning of 'a' and 'b' is different there. */
  static rsAffineTransform2D<RealType> reflection(RealType a, RealType b)
  {
    return rsAffineTransform2D<RealType>::reflection(-a, 1, -b);
  }

  /** Creates the unique affine transform that maps 3 given points \f$ p_1, p_2, p_3 \f$ into
  the 3 given image points \f$ q_1, q_2, q_3 \f$ such that
  \f[ \mathbf{q}_i = \mathbf{A p}_i + \mathbf{b}, \quad i = 1,2,3 \f]
  */
  static rsAffineTransform2D<RealType> from3PointsAndImages(
    const rsPoint2D<RealType> &p1, const rsPoint2D<RealType> &p2, const rsPoint2D<RealType> &p3,
    const rsPoint2D<RealType> &q1, const rsPoint2D<RealType> &q2, const rsPoint2D<RealType> &q3)
  {
    rsAffineTransform2D<RealType> t;
    RealType s = 1 / p1.x*(p3.y-p2.y)-p2.x*p3.y+p2.y*p3.x+p1.y*(p2.x-p3.x);
    t.a11 = -s * (p1.y*(q3.x-q2.x)-p2.y*q3.x+p3.y*q2.x+(p2.y-p3.y)*q1.x);
    t.a12 =  s * (p1.x*(q3.x-q2.x)-p2.x*q3.x+p3.x*q2.x+(p2.x-p3.x)*q1.x);
    t.a21 =  s * (p2.y*q3.y+p1.y*(q2.y-q3.y)-p3.y*q2.y+(p3.y-p2.y)*q1.y);
    t.a22 = -s * (p2.x*q3.y+p1.x*(q2.y-q3.y)-p3.x*q2.y+(p3.x-p2.x)*q1.y);
    t.b1  =  s * (p1.x*(p3.y*q2.x-p2.y*q3.x)+p1.y*(p2.x*q3.x-p3.x*q2.x)
      +(p2.y*p3.x-p2.x*p3.y)*q1.x);
    t.b2  =  s * (p1.x*(p3.y*q2.y-p2.y*q3.y)+p1.y*(p2.x*q3.y-p3.x*q2.y)
      +(p2.y*p3.x-p2.x*p3.y)*q1.y);
// \todo: precompute common subexpressions (maybe)
    return t;
  }


  /** \name Inquiry */

  /** Returns the determinant of this transform. */
  RealType getDeterminant() const
  {
    return a11*a22-a21*a12;
  }

  /** Returns true when this transform is singular - this essentially means that the plane is
  collapsed into a line or point. */
  bool isSingular() const
  {
    return getDeterminant() == 0;
  }


  /** \name Misc */

  /** Returns a transformation that represets this transformation followed by the transformation
  given in the argument. */
  rsAffineTransform2D<RealType> followedBy(const rsAffineTransform2D<RealType> &t)
  {
    rsAffineTransform2D<RealType> result;
    result.a11 = t.a11*a11 + t.a12*a21;
    result.a12 = t.a11*a12 + t.a12*a22;
    result.a21 = t.a21*a11 + t.a22*a21;
    result.a22 = t.a21*a12 + t.a22*a22;
    result.b1  = t.a11*b1  + t.a12*b2 + t.b1;
    result.b2  = t.a21*b1  + t.a22*b2 + t.b2;
    return result;
  }

  /** Returns the inverse transform of this transform. If the transform is given by:
  \f[ \mathbf{y = Ax + b} \f] then the inverse transform is:
  \f[ \mathbf{x = A^{-1} (y-b) = A^{-1}y - A^{-1}b} \f]
  \todo rename to getInverse
  */
  rsAffineTransform2D<RealType> inverted()
  {
    rsAssert(!isSingular()); // only non-singular transforms can be inverted

    rsAffineTransform2D<RealType> r;
    RealType scaler = 1 / getDeterminant();
    r.a11 =  scaler * a22;
    r.a12 = -scaler * a12;
    r.a21 = -scaler * a21;
    r.a22 =  scaler * a11;
    r.b1  = -(r.a11*b1 + r.a12*b2);
    r.b2  = -(r.a21*b1 + r.a22*b2);
    return r;
  }


  /** \name Application */

  /** Applies the transform to a point represented by coordinates x, y. */
  void applyTo(RealType &x, RealType &y);

  /** Applies the transform to an original point p and returns its image point.
  \todo maybe rename to applyTo
  */
  rsPoint2D<RealType> getTransformedPoint(const rsPoint2D<RealType> &originalPoint)
  {
    rsPoint2D<RealType> p = originalPoint;
    applyTo(p.x, p.y);
    return p;
  }


protected:

  /** \name Data */

  RealType a11, a12, a21, a22, b1, b2;

};

//-----------------------------------------------------------------------------------------------
// implementation:

/** Applies the transform to a point represented by coordinates x, y such that:
\f[
\begin{pmatrix}
x' \\
y'
\end{pmatrix}
=
\begin{pmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{pmatrix}
\begin{pmatrix}
x \\
y
\end{pmatrix}
+
\begin{pmatrix}
b_1 \\
b_2
\end{pmatrix}
\f] */
template<class RealType>
void rsAffineTransform2D<RealType>::applyTo(RealType &x, RealType &y)
{
  RealType xTmp = x;
  x = a11*xTmp + a12*y + b1;
  y = a21*xTmp + a22*y + b2;
}

#endif
