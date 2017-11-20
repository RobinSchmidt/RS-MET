#ifndef RS_COMPLEX_H
#define RS_COMPLEX_H

namespace RSLib
{

  /**

  This is a class for complex numbers. It defines the basic arithmetic operations between complex
  numbers as well as the special cases when one of the operands is real.

  ATTENTION: do not define any further member variables, nor switch the ordering of re and im
  because the FourierTransformer classes rely on the fact, that an rsComplex<double> number
  consists of two doubles re, im (in that order) and nothing else (the low-level algorithms
  actually run on buffers of doubles). It's generally a bit dangerous to make assumptions about
  the memory layout of member variables, but for performance reasons we do it anyway. The
  unit-test suite has a function for checking the sizes and memory-layout of Complex<double> and
  Complex<float> and will complain when these assumptions fail to hold.

  \todo: try hypot() to calculate the absolute value (might be faster and/or more accurate than
         current implementation?)
  \todo: write LaTeX compatible comments with math notation
  \todo: use a sinCos functions for cartesian-to-polar conversion, maybe write a function for this
         conversion independent of class Complex
  \todo: replace the division, and absolute-value calculation by the more accurate algorithms
         from Numerical Recipies
  \todo: in the division operators, we first compute a scale-factor = 1 / (stuff) in order to 
         replace a division in subsequent operations by multiplication. that, however, rules out to 
         use the template for integers (if division is to be used). maybe get rid of the 
         precomputation and leave that optimization to the compiler. but we should do performance 
         tests, if this really is inconsequential (i hope so)

  */

  template<class RealType>
  class rsComplex
  {

  public:

    /** \name Public Data Members */

    /** Real part */
    RealType re;

    /** Imaginary part */
    RealType im;


    /** \name Construction/Destruction */

    /** Constructor. Initializes real and imaginary part to zero. */
    rsComplex() { re = im = RealType(0); }

    /** Constructor. Initializes real part to the argument "reInit" and imaginary part to zero. */
    rsComplex(RealType reInit)
    {
      re = RealType(reInit);
      im = RealType(0);
    }

    /** Constructor. Initializes real and imaginary parts with the parameters. */
    rsComplex(RealType reInit, RealType imInit)
    {
      re = RealType(reInit);
      im = RealType(imInit);
    }

    /** Constructor for implicit type-conversion of the complex number with RealType2 into a
    complex number with RealType. */
    template<class RealType2>
    rsComplex(rsComplex<RealType2> zInit)
    {
      re = (RealType) zInit.re;
      im = (RealType) zInit.im;
    }


    /** \name Operators */

    /** Compares two complex numbers of equality. */
    bool operator==(const rsComplex& z) const
    {
      if( re == z.re && im == z.im )
        return true;
      else
        return false;
    }

    /** Compares two complex numbers of inequality. */
    bool operator!=(const rsComplex& z) const
    {
      if( re != z.re || im != z.im )
        return true;
      else
        return false;
    }

    /** Defines the negative of a complex number. */
    rsComplex operator-() const
    { return rsComplex(-re, -im); }

    /** Adds another complex number to this complex and returns the result. */
    rsComplex& operator+=(const rsComplex &z)
    {
      re += z.re;
      im += z.im;
      return *this;
    }

    /** Adds a real number to this complex and returns the result. */
    rsComplex& operator+=(const RealType &r)
    {
      re += r;
      return *this;
    }

    /** Subtracts another complex number from this complex and returns the result. */
    rsComplex& operator-=(const rsComplex &z)
    {
      re -= z.re;
      im -= z.im;
      return *this;
    }

    /** Subtracts a real number from this complex and returns the result. */
    rsComplex& operator-=(const RealType &r)
    {
      re -= r;
      return *this;
    }

    /** Multiplies this complex number by another complex number and returns the result. */
    rsComplex& operator*=(const rsComplex &z)
    {
      RealType reNew = re*z.re - im*z.im;
      RealType imNew = re*z.im + im*z.re; // ...we could assign this directly to im
      re             = reNew;
      im             = imNew;
      return *this;
    }

    /** Multiplies this complex number by a real number and returns the result. */
    rsComplex& operator*=(const RealType &r)
    {
      re *= r;
      im *= r;
      return *this;
    }

    /** Divides this complex number by another complex number and returns the result. */
    rsComplex& operator/=(const rsComplex &z)
    {
      RealType scale = 1.0 / (z.re*z.re + z.im*z.im);
      RealType reNew = scale*( re*z.re  + im*z.im  );
      RealType imNew = scale*( im*z.re  - re*z.im  );
      re             = reNew;
      im             = imNew;
      return *this;
    }

    /** Divides this complex number by a real number and returns the result. */
    rsComplex& operator/=(const RealType &r)
    {
      RealType scale = 1.0 / r;
      re *= scale;
      im *= scale;
      return *this;
    }


    /** \name Setup */

    /** Adjusts the radius of this complex number leaving the angle unchanged. */
    void setRadius(RealType newRadius)
    {
      RealType phi = getAngle();
      im  = sin(phi);
      re  = cos(phi);     // \todo: use: sinCos(phi, &im, &re);
      re *= newRadius;    // re = newRadius * cos(phi);
      im *= newRadius;    // im = newRadius * sin(phi);
    }

    /** Adjusts the angle of this complex number leaving the magnitude unchanged. */
    void setAngle(RealType newAngle)
    {
      RealType r = getRadius();
      im = sin(newAngle);
      re = cos(newAngle);    // \todo: use sinCos(newAngle, &im, &re);
      re *= r;               // re = r * cos(newAngle);
      im *= r;               // im = r * sin(newAngle);
    }

    /** Sets the radius and angle of this complex number. */
    void setRadiusAndAngle(RealType newRadius, RealType newAngle)
    {
      im  = sin(newAngle);
      re  = cos(newAngle);       // \todo: use sinCos(newAngle, &im, &re);
      re *= newRadius;           // re = newRadius * cos(newAngle);
      im *= newRadius;           // im = newRadius * sin(newAngle);
    }


    /** \name Inquiry */

    /** Returns the radius of this complex number. */
    RealType getRadius() const
    {
      return rsSqrt(re*re + im*im); // try hypot, measure performance
      //return hypot(re, im);     // actually slower
      // \todo: use algorithm from numerical recipies (more precise) - maybe use a separate
      // function getRadiusPrecise or something
    }

    /** Returns the angle of this complex number. */
    RealType getAngle() const
    {
      if((re==0.0) && (im==0))
        return 0.0;
      else
        return atan2(im, re);
    }

    /** Returns the complex conjugate of this complex number. */
    rsComplex getConjugate() const
    {
      return rsComplex(re, -im);
    }

    /** Returns the reciprocal of this complex number. */
    rsComplex getReciprocal() const
    {
      RealType scaler = 1.0 / (re*re + im*im);
      return Complex(scaler*re, -scaler*im);
    }

    /** Returns true, if this complex number is purely real. */
    bool isReal() const
    {
      return im == RealType(0);
    }

    /** Returns true, if this complex number is purely imaginary. */
    bool isImaginary() const
    {
      return re == RealType(0);
    }

    /** Returns true if real or imaginary part (or both) are plus or minus infinity, false
    otherwise. */
    bool isInfinite() const
    {
      if(  re == RS_INF(RealType) || re == -RS_INF(RealType)
        || im == RS_INF(RealType) || im == -RS_INF(RealType) )
        return true;
      else
        return false;
    }

  };

  // some binary operators are defined outside the class, such that the left hand operand does
  // not necessarily need to be of class rsComplex:

  /** Adds two complex numbers. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator+(const rsComplex<RealType> &z, const rsComplex<RealType> &w)
  {
    return rsComplex<RealType>(z.re+w.re, z.im+w.im);
  }

  /** Adds a complex and a real number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator+(const rsComplex<RealType> &z, const RealType &r)
  {
    return rsComplex<RealType>(z.re+r, z.im);
  }

  /** Adds a real and a complex number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator+(const RealType &r, const rsComplex<RealType> &z)
  {
    return rsComplex<RealType>(z.re+r, z.im);
  }

  /** Subtracts two complex numbers. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator-(const rsComplex<RealType> &z, const rsComplex<RealType> &w)
  {
    return rsComplex<RealType>(z.re-w.re, z.im-w.im);
  }

  /** Subtracts a real number from a complex number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator-(const rsComplex<RealType> &z, const RealType &r)
  {
    return rsComplex<RealType>(z.re-r, z.im);
  }

  /** Subtracts a complex number from a real number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator-(const RealType &r, const rsComplex<RealType> &z)
  {
    return rsComplex<RealType>(r-z.re, -z.im);
  }

  /** Multiplies two complex numbers. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator*(const rsComplex<RealType> &z, const rsComplex<RealType> &w)
  {
    return rsComplex<RealType>(z.re*w.re-z.im*w.im, z.re*w.im+z.im*w.re);
  }

  /** Multiplies a complex number and a real number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator*(const rsComplex<RealType> &z, const RealType &r)
  {
    return rsComplex<RealType>(z.re*r, z.im*r);
  }

  /** Multiplies a real number and a complex number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator*(const RealType &r, const rsComplex<RealType> &z)
  {
    return rsComplex<RealType>(z.re*r, z.im*r);
  }

  /** Divides two complex numbers. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator/(const rsComplex<RealType> &z, const rsComplex<RealType> &w)
  {
    return rsComplex<RealType>( ( z.re*w.re + z.im*w.im) / (w.re*w.re + w.im*w.im),
                                ( z.im*w.re - z.re*w.im) / (w.re*w.re + w.im*w.im)  );

    // the division-optimization below causes problems with nested complex numbers (1 nesting level
    // works, but 2 don't):
    //RealType scale = rsComplex<RealType>(1.0) / (w.re*w.re + w.im*w.im);
    //return rsComplex<RealType>( scale*( z.re*w.re + z.im*w.im),     // real part
    //                            scale*( z.im*w.re - z.re*w.im)  );  // imaginary part
  }

  /** Divides a complex number by a real number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator/(const rsComplex<RealType> &z, const RealType &r)
  {
    RealType scale = 1.0 / r;
    return rsComplex<RealType>(scale*z.re, scale*z.im);
  }

  /** Divides a real number by a complex number. */
  template<class RealType>
  RS_INLINE rsComplex<RealType> operator/(const RealType &r, const rsComplex<RealType> &z)
  {
    RealType scale = r / (z.re*z.re + z.im*z.im);
    return rsComplex<RealType>(scale*z.re, -scale*z.im);
  }

  /* Returns the (primitive) square-root of a complex number (algorithm from Numerical Recipies).
  \todo: move to a sepaprate file ComplexFunctions.h
  */
  template<class RealType>
  rsComplex<RealType> rsSqrtC(rsComplex<RealType> z)
  {
    rsComplex<RealType> c;
    RealType x, y, w, r;
    if( (z.re == 0.0) && (z.im == 0.0) )
    {
      c.re = 0.0;
      c.im = 0.0;
      return c;    // return rsComplex<RealType>(0, 0);
    }
    else
    {
      x = fabs(z.re);  // use rsAbs
      y = fabs(z.im);
      if( x >= y )
      {
        r = y/x;
        w = rsSqrt(x) * rsSqrt(0.5*(1.0+rsSqrt(1.0+r*r)));
      }
      else
      {
        r = x/y;
        w = rsSqrt(y) * rsSqrt(0.5*(r+rsSqrt(1.0+r*r)));
      }
      if( z.re >= 0.0 )
      {
        c.re = w;
        c.im = z.im/(2.0*w);
      }
      else
      {
        c.im = (z.im >= 0) ? w : -w; // if(z.im >= 0) c.im = w; else c.im = -w;
        c.re = z.im/(2.0*c.im);
      }
      return c;
    }
  }

  /** Implements the generic rsAbs function for complex numbers. */
  template<class RealType>
  RealType rsAbs(rsComplex<RealType> z)
  {
    return z.getRadius();
  }

  /* Computes the 2 roots (zeros) of the quadratic equation x^2 + p*x + q = 0. */
  template<class RealType>
  void solveQuadraticEquation(rsComplex<RealType> p, rsComplex<RealType> q,
                              rsComplex<RealType> &root1, rsComplex<RealType> &root2)
  {
    rsComplex<RealType> tmp = rsSqrtC(0.25*p*p - q);
    root1 = -0.5*p + tmp;
    root2 = -0.5*p - tmp;
  }

  // explicit instantiations:
  typedef rsComplex<double> RSLib_API rsComplexDbl;
}

#endif
