#ifndef RAPT_MOEBIUSTRANSFORM_H_INCLUDED
#define RAPT_MOEBIUSTRANSFORM_H_INCLUDED

/** This is a class for representing Moebius transformations in the complex plane. A Moebius
transformation is a defined as a transformation of the general form:

  \f[ w = M(z) = \frac{a z + b}{c z + d} \f]

where \f$ a, b, c, d \f$ are complex constants.



\todo implement normalization, implement fixpoint calculation
      (we need the rsSqrtC function for that), test it
\todo maybe make a baseclass rsConformalMap 

*/

template<class T>
class rsMoebiusTransform
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  /** Constructor. Initializes to the identity transform \f$ a = d = 1, b = c = 0 \f$. */
  rsMoebiusTransform() { setCoefficients(T(1), T(0), T(0), T(1)); }

  /** Constructor. Initializes the coefficients with the given values.  */
  rsMoebiusTransform(const rsComplex<T>& a, const rsComplex<T>& b,
    const rsComplex<T>& c, const rsComplex<T>& d)
  { setCoefficients(a, b, c, d); }

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setCoefficients(rsComplex<T> a, rsComplex<T> b,
    rsComplex<T> c, rsComplex<T> d)
  { this->a = a; this->b = b; this->c = c; this->d = d; }

  /* Normalizes this transform such that the determinant becomes unity. This will determine the
  coefficients up to an inversion of sign. not pinpoint the coefficients uniquely.  */
  void normalize()
  {
    rsComplex<T> s, phi;                // what is phi for?
    s = T(1) / rsSqrt(getDeterminant());
    a *= s;
    b *= s;
    c *= s;
    d *= s;
  }

  /*
  void forceRealPositiveA()
  {
    T r = a.getRadius();
    rsComplex<T> s = a.getConjugate() / r;
    a *= s;
    b *= s;
    c *= s;
    d *= s;
  }
  */

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  inline rsComplex<T> getA() const { return a; }
  inline rsComplex<T> getB() const { return b; }
  inline rsComplex<T> getC() const { return c; }
  inline rsComplex<T> getD() const { return d; }
  // remove, make a,b,c,d public

  /** Returns the determinant of this transform. */
  rsComplex<T> getDeterminant() const { return a*d - b*c; }

  /** Returns true when this transform is singular. This means that the whole complex plane will
  be collapsed the single point \f$ a/c \f$ by application of the trnasformation. */
  bool isSingular() const { return getDeterminant() == rsComplex<T>(0,0); }

  /** Returns true when this Moebius transform is an identity transformation. The condition for
  this to be the case is \f$ a = d \neq 0, b = c = 0 \f$. */
  bool isIdentity() const
  {
    rsComplex<T> O(0, 0); // zero, for comparison
    if(b != O || c != O || a != d)
      return false;
    if(a == O)
      return false;
    return true;
  }

  /** Returns the 2 fixpoints of this transform which are points that are mapped to
  themselves such that z = M(z) = (a*z + b) / (c*z + d). The fixpoints are solutions of the
  quadratic equation c*z^2 + (d-a)*z - b = 0 */
  void getFixPoints(rsComplex<T> &fp1, rsComplex<T> &fp2) const
  {
    rsSolveQuadraticEquation((d-a)/c, -b/c, fp1, fp2);
  }
  // rename to getFixedPoints

  /* Returns the multiplier associated with the 1st fixpoint (the one returned in fp1 in
  getFixPoints). The mulitplier associated with the other fixpoint fp2 is the reciprocal
  of this value. */
  /*
  rsComplex<T> getMultiplier()
  {
    rsComplex<T> fp1, fp2;
    getFixPoints(fp1, fp2);
    return (a - c*fp1) / (a - c*fp2);  // Eq.42
    // ckeck, if denominator may become zero - maybe use Eq.43 instead
  }
  */

  /*
  bool isElliptic()   const { return (a+d).isReal() && (a+d).getRadius() < 2; }
  bool isHyperbolic() const { return (a+d).isReal() && (a+d).getRadius() > 2; }
  bool isParabolic()  const { return (a+d) == 2     || (a+d) == -2; }
  bool isLoxodromic() const { return !(a+d).isReal(); }
  wrong: these equations hold only for a normalized Moebius transform, maybe compute the
  "multiplier" of the transform and check its real and imaginary parts or something
  */

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Compares two transforms for equality. Equality is defined here as the condition, that all
  coefficients should be equal, even though a Moebius transform is unique only up to an arbitrary
  (non-zero) scaling of all coefficients. */
  bool operator==(const rsMoebiusTransform& t) const
  {
    return a == t.a && b == t.b && c == t.c && d == t.d;
    /*
    if(a == t.a && b == t.b && c == t.c && d == t.d)
      return true;
    else
      return false;
      */
  }

  /** Compares two transformations for inequality. */
  bool operator!=(const rsMoebiusTransform& t) const
  {
    return !(*this == t);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Factory Functions (move to Lifetime) */

  /** Returns an identity transform. */
  static const rsMoebiusTransform<T> identity()
  {
    return rsMoebiusTransform<T>(1, 0, 0, 1);
  }

  /** Returns a transform that scales a complex number \f$ z \f$ by some complex constant 
  \f$ k \f$. */
  static const rsMoebiusTransform<T> complexScaling(rsComplex<T> k)
  {
    return rsMoebiusTransform<T>(k, 0, 0, 1);
  }

  // \todo:
  // scaling(rsComplex<T> k) { a = Complex(k, 0); b = c = 0; d = 1; }
  // rotation: a = exp(i*phi); b = c = 0; d = 1; -> complexScaling(exp(i*phi))
  // inversion: a = d = 0; b = c = 1;
  // translation: a = 1; b = u; c = 0; d = 1; // u is the translation vector/number
  // rotationByHalfPi: complexScaling(rsComplex<T>(0, 1))

  /** Creates the unique (up to scaling of the coefficients) Moebius transform that maps 3 given
  points \f$ z_1, z_2, z_3 \f$ to the 3 given image points \f$ w_1, w_2, w_3 \f$ such that
  \f[ w_i = M(z_i), \quad i = 1,2,3 \f]  */
  static rsMoebiusTransform<T> from3PointsAndImages(
    const rsComplex<T> &z1, const rsComplex<T> &z2, const rsComplex<T> &z3,
    const rsComplex<T> &w1, const rsComplex<T> &w2, const rsComplex<T> &w3)
  {
    rsMoebiusTransform<T> r;
    r.a = (w2-w1)*w3*z3+(w1*w2-w2*w3)*z2+(w1*w3-w1*w2)*z1;
    r.b = ((w1*w3-w1*w2)*z2+(w1*w2-w2*w3)*z1)*z3+(w2-w1)*w3*z1*z2;
    r.c = (w2-w1)*z3+(w1-w3)*z2+(w3-w2)*z1;
    r.d = ((w3-w2)*z2+(w1-w3)*z1)*z3+(w2-w1)*z1*z2;
    return r;
    // \todo: optimize (by factoring out and precomputation of common subexpressions)
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Related Transformations */

  /** Returns a transformation that represets this transformation followed by the transformation
  given in the argument. */
  rsMoebiusTransform<T> followedBy(const rsMoebiusTransform<T> &t)
  {
    rsMoebiusTransform<T> r;
    r.a = t.a*a + t.b*c;
    r.b = t.a*b + t.b*d;
    r.c = t.c*a + t.d*c;
    r.d = t.c*b + t.d*d;
    return r;
  }
  // ToDo: 
  // -Verify, if this formula is still correct when the multiplications t.a*a, ..etc. are
  //  non-commutative. If not, maybe some of the arguments should be swapped. Then try this with a 
  //  type T that has a non-commutative multiplication such as rsMatrix2D. I have no idea, if 
  //  Moebius transformations with non-commutative underlying types are a thing, though - but 
  //  better be safe than sorry. But I think, it's actually the correct order, if "t" is the left
  //  operand and "this" is the right operand which is the case in function composition.

  /** Returns the inverse transform of this transform. */
  rsMoebiusTransform<T> getInverse()
  {
    rsAssert(!isSingular()); // only non-singular transforms can be inverted
    rsMoebiusTransform<T> r;
    r.a =  d;
    r.b = -b;
    r.c = -c;
    r.d =  a;
    return r;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Application */

  /** Applies the transform to the given complex number. */
  void applyTo(rsComplex<T> &z) const;
  // todo: return the value instead of manipulating the input argument, rename to apply - or just
  // implement the () operator

  /** Applies the transform to the given complex number and returns the mapped number w. */
  rsComplex<T> getMappedNumber(const rsComplex<T> &z) const
  {
    rsComplex<T> w = z;
    applyTo(w);
    return w;
  }
  // get rid

protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  rsComplex<T> a, b, c, d;  // use rsComplex

};

//-------------------------------------------------------------------------------------------------
// implementation:

template<class T>
void rsMoebiusTransform<T>::applyTo(rsComplex<T> &z) const
{
  z = (a*z + b) / (c*z + d);
}
// move into class

#endif
