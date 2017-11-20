#ifndef RS_MOEBIUSTRANSFORM_H
#define RS_MOEBIUSTRANSFORM_H

namespace RSLib
{

  /**

  This is a class for representing Moebius transformations in the complex plane. A Moebius 
  transformation is a defined as a transformation of the general form: 

  \f[ w = M(z) = \frac{a z + b}{c z + d} \f]

  where \f$ a, b, c, d \f$ are complex constants.

  \todo implement normalization, implement fixpoint calculation 
        (we need the rsSqrtC function for that), test it
  \todo maybe make a baseclass ConformalMap

  */

  template<class RealType>
  class rsMoebiusTransform
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Initializes to the identity transform \f$ a = d = 1, b = c = 0 \f$. */
    rsMoebiusTransform()
    {
      setCoefficients(1, 0, 0, 1);
    }

    /** Constructor. Initializes the coefficients with the given values.  */
    rsMoebiusTransform(rsComplex<RealType> a, rsComplex<RealType> b, 
                       rsComplex<RealType> c, rsComplex<RealType> d)
    {
      setCoefficients(a, b, c, d);
    }


    /** \name Setup */

    void setCoefficients(rsComplex<RealType> a, rsComplex<RealType> b, 
                         rsComplex<RealType> c, rsComplex<RealType> d)
    {
      this->a = a; 
      this->b = b;
      this->c = c; 
      this->d = d;
    }

    /* Normalizes this transform such that the determinant becomes unity. This will determine the
    coefficients up to an inversion of sign.
    not pinpoint the coefficients uniquely.  */
    void normalize()
    {
      rsComplex<RealType> s, phi;
      s = 1.0 / rsSqrtC(getDeterminant());
      a *= s;
      b *= s;
      c *= s;
      d *= s;
    }

    /*
    void forceRealPositiveA()
    {
      RealType r = a.getRadius();
      rsComplex<RealType> s = a.getConjugate() / r;
      a *= s;
      b *= s;
      c *= s;
      d *= s;
    }
    */


    /** \name Inquiry */

    inline rsComplex<RealType> getA() const { return a; }
    inline rsComplex<RealType> getB() const { return b; }
    inline rsComplex<RealType> getC() const { return c; }
    inline rsComplex<RealType> getD() const { return d; }

    /** Returns the determinant of this transform. */
    rsComplex<RealType> getDeterminant() const
    {
      return a*d - b*c;
    }

    /** Returns true when this transform is singular. This means that the whole complex plane will
    be collapsed the single point \f$ a/c \f$ by application of the trnasformation. */
    bool isSingular() const
    {
      return getDeterminant() == 0;
    }

    /** Returns true when this Moebius transform is an identity transformation. The condition for 
    this to be the case is \f$ a = d \neq 0, b = c = 0 \f$. */
    bool isIdentity() const
    {
      if( b != 0 || c != 0 || a != d )
        return false;
      if( a == 0 )
        return false;
      return true;
    }

    /** Returns the 2 fixpoints of this transform which are points that are mapped to 
    themselves such that z = M(z) = (a*z + b) / (c*z + d). The fixpoints are solutions of the 
    quadratic equation c*z^2 + (d-a)*z - b = 0 */
    void getFixPoints(rsComplex<RealType> &fp1, rsComplex<RealType> &fp2) const
    {
      solveQuadraticEquation((d-a)/c, -b/c, fp1, fp2);
    }

    /* Returns the multiplier associated with the 1st fixpoint (the one returned in fp1 in 
    getFixPoints). The mulitplier associated with the other fixpoint fp2 is the reciprocal 
    of this value. */
    /*
    rsComplex<RealType> getMultiplier()
    {
      rsComplex<RealType> fp1, fp2;
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


    /** \name Operators */

    /** Compares two transforms for equality. Equality is defined here as the condition, that all 
    coefficients should be equal, even though a Moebius transform is unique only up to an arbitrary
    (non-zero) scaling of all coefficients. */
    bool operator==(const rsMoebiusTransform& t) const
    {
      if( a == t.a && b == t.b && c == t.c && d == t.d )
        return true;
      else
        return false;
    }

    /** Compares two transformations for inequality. */
    bool operator!=(const rsMoebiusTransform& t) const
    {
      return !(*this == t);
    }


    /** \name Factory Functions */

    /** Returns an identity transform. */
    static const rsMoebiusTransform<RealType> identity()
    {
      return rsMoebiusTransform<RealType>(1, 0, 0, 1);
    }

    /** Returns a transform that sacles a complex number \f$ z \f$ by some complex constant 
    \f$ k \f$. */
    static const rsMoebiusTransform<RealType> complexScaling(rsComplex<RealType> k)
    {
      return rsMoebiusTransform<RealType>(k, 0, 0, 1);
    }

    // \todo:
    // scaling(rsComplex<RealType> k) { a = Complex(k, 0); b = c = 0; d = 1; }
    // rotation: a = exp(i*phi); b = c = 0; d = 1; -> complexScaling(exp(i*phi))
    // inversion: a = d = 0; b = c = 1;
    // translation: a = 1; b = u; c = 0; d = 1; // u is the translation vector/number
    // rotationByHalfPi: complexScaling(rsComplex<RealType>(0, 1))

    /** Creates the unique (up to scaling of the coefficients) Moebius transform that maps 3 given 
    points \f$ z_1, z_2, z_3 \f$ to the 3 given image points \f$ w_1, w_2, w_3 \f$ such that
    \f[ w_i = M(z_i), \quad i = 1,2,3 \f]  */
    static rsMoebiusTransform<RealType> from3PointsAndImages(
      const rsComplex<RealType> &z1, const rsComplex<RealType> &z2, const rsComplex<RealType> &z3, 
      const rsComplex<RealType> &w1, const rsComplex<RealType> &w2, const rsComplex<RealType> &w3)
    {
      rsMoebiusTransform<RealType> r;
      r.a = (w2-w1)*w3*z3+(w1*w2-w2*w3)*z2+(w1*w3-w1*w2)*z1;         
      r.b = ((w1*w3-w1*w2)*z2+(w1*w2-w2*w3)*z1)*z3+(w2-w1)*w3*z1*z2;
      r.c = (w2-w1)*z3+(w1-w3)*z2+(w3-w2)*z1;
      r.d = ((w3-w2)*z2+(w1-w3)*z1)*z3+(w2-w1)*z1*z2;
      return r;
      // \todo: optimize (by factoring out and precomputation of common subexpressions)
    }


    /** \name Related Transformations */

    /** Returns a transformation that represets this transformation followed by the transformation 
    given in the argument. */
    rsMoebiusTransform<RealType> followedBy(const rsMoebiusTransform<RealType> &t)
    {
      rsMoebiusTransform<RealType> r;
      r.a = t.a*a + t.b*c;
      r.b = t.a*b + t.b*d;
      r.c = t.c*a + t.d*c;
      r.d = t.c*b + t.d*d;
      return r;
    }

    /** Returns the inverse transform of this transform. */
    rsMoebiusTransform<RealType> getInverse()
    {
      rsAssert( !isSingular() ); // only non-singular transforms can be inverted
      rsMoebiusTransform<RealType> r;
      r.a =  d;
      r.b = -b;
      r.c = -c;
      r.d =  a;
      return r;
    }


    /** \name Application */

    /** Applies the transform to the given complex number. */
    void applyTo(rsComplex<RealType> &z) const;

    /** Applies the transform to the given complex number and returns the mapped number w. */
    rsComplex<RealType> getMappedNumber(const rsComplex<RealType> &z) const
    {
      rsComplex<RealType> w = z;
      applyTo(w);
      return w;
    }

  protected:

    /** \name Data */

    rsComplex<RealType> a, b, c, d;

  };

  //-----------------------------------------------------------------------------------------------
  // implementation:

  template<class RealType>
  void rsMoebiusTransform<RealType>::applyTo(rsComplex<RealType> &z) const
  {
    z = (a*z + b) / (c*z + d);
  }

}

#endif
