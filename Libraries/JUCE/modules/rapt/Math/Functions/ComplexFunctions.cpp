#ifndef RS_COMPLEXFUNCTIONS_INL
#define RS_COMPLEXFUNCTIONS_INL

#include "ComplexFunctions.h"

namespace RSLib
{

  template<class T>
  rsComplex<T> rsAcdC(rsComplex<T> w, double k)
  {
    if( k == 1.0 )
    {
      RS_DEBUG_BREAK;    // k may not be equal to 1
      return rsComplex<T>(0.0, 0.0);
    }

    int    M = 7;     // fixed number of Landen iterations
    double v[7];      // array to store the vector of descending elliptic moduli
    rsLanden(k, M, v);

    int    n;
    double v1;
    for(n = 0; n < M; n++)
    {
      if( n == 0 )
        v1 = k;
      else
        v1 = v[n-1];
      w = w / (1.0 + rsSqrtC(1.0 - w*w * v1*v1)) * 2.0/(1.0+v[n]);
    }

    rsComplex<T> u;
    if( w == 1.0 )
      u = 0.0;
    else
      u = (2.0/PI) * rsAcosC(w);

    double K, Kprime;
    rsEllipticIntegral(k, &K, &Kprime);
    double R = Kprime/K;

    u.re = rsSrem(u.re, 4.0);
    u.im = rsSrem(u.im, 2.0*R);

    return u;
  }

  template<class T>
  rsComplex<T> rsAcosC(rsComplex<T> z)
  {
    rsComplex<T> tmp = z + rsComplex<T>(0.0, 1.0) * rsSqrtC( 1.0 - z*z );
    tmp = rsLogC(tmp);
    return tmp * rsComplex<T>(0.0, -1.0);
  }

  template<class T>
  rsComplex<T> rsAcoshC(rsComplex<T> z)
  {
    rsComplex<T> tmp = z + rsSqrtC( z - 1.0 ) * rsSqrtC( z + 1.0 );
    return rsLogC(tmp);
  }

  template<class T>
  rsComplex<T> rsAsinC(rsComplex<T> z)
  {
    rsComplex<T> tmp = z * rsComplex<T>(0.0, 1.0) + rsSqrtC( 1.0 - z*z );
    tmp = rsLogC(tmp);
    return tmp * rsComplex<T>(0.0, -1.0);
  }

  template<class T>
  rsComplex<T> rsAsinhC(rsComplex<T> z)
  {
    rsComplex<T> tmp = z + rsSqrtC( z*z + 1.0 ) ;
    return rsLogC(tmp);
  }

  template<class T>
  rsComplex<T> rsAsnC(rsComplex<T> w, double k)
  {
    return 1.0 - rsAcdC(w, k);
  }

  template<class T>
  rsComplex<T> rsAtanC(rsComplex<T> z)
  {
    rsComplex<T> tmp = ( rsComplex<T>(0.0, 1.0 ) + z) / ( rsComplex<T>(0.0, 1.0) - z );
    tmp = rsLogC(tmp);
    return tmp * rsComplex<T>(0.0, -0.50);
  }

  template<class T>
  rsComplex<T> rsAtanhC(rsComplex<T> z)
  {
    rsComplex<T> tmp = ( rsComplex<T>(1.0, 0.0 ) + z ) / ( rsComplex<T>(1.0, 0.0) - z );
    tmp = rsLogC(tmp);
    return 0.5 * tmp;
  }

  template<class T>
  rsComplex<T> rsCdC(rsComplex<T> u, double k)
  {
    int    M = 7;     // fixed number of Landen iterations
    double v[7];      // array to store the vector of descending elliptic moduli
    rsLanden(k, M, v);

    // initialization:
    rsComplex<T> w = rsCosC(u * PI/2.0);

    // ascending Landen/Gauss transformation:
    for(int n=M-1; n>=0; n--)
      w = (1.0+v[n])*w / (1.0+v[n]*w*w);

    return w;
  }

  template<class T>
  rsComplex<T> rsCosC(rsComplex<T> z)
  {
    rsComplex<T> tmp = rsExpC(rsComplex<T>(0.0, 1.0) * z); // tmp = exp(i*z);
    tmp        += rsExpC(rsComplex<T>(0.0, -1.0)* z);      // tmp = exp(i*z) + exp(-i*z)
    return 0.5 * tmp;
  }

  template<class T>
  rsComplex<T> rsCoshC(rsComplex<T> z)
  {
    return 0.5 * ( expC(z) + expC(-z) );
  }

  template<class T>
  rsComplex<T> rsExpC(rsComplex<T> z)
  {
    rsComplex<T> result;
    T tmp = exp(z.re);
    rsSinCos(z.im, &result.im, &result.re);
    result.re *= tmp;
    result.im *= tmp;
    return result;
  }

  template<class T>
  rsComplex<T> rsLogC(rsComplex<T> z)
  {
    rsComplex<T> tmp;
    tmp.re = log(z.getRadius());
    tmp.im = z.getAngle();
    return tmp;
  }

  template<class T>
  rsComplex<T> rsPowC(rsComplex<T> basis, rsComplex<T> exponent)
  {
    if( basis != rsComplex<T>(0.0, 0.0) )
      return rsExpC( exponent * rsLogC((rsComplex<T>)basis) );
    else // basis == 0
    {
      if( exponent.re == 0.0 && exponent.im == 0.0 )
        return rsComplex<T>(1.0, 0.0);
      else if( exponent.re < 0.0 && exponent.im == 0.0 )
        return rsComplex<T>(RS_INF(T), 0.0);
      else if( exponent.re > 0.0 && exponent.im == 0.0 )
        return rsComplex<T>(0.0, 0.0);
      else
        return rsExpC( exponent * rsLogC((rsComplex<T>)basis) );
    }
  }

  template<class T>
  rsComplex<T> rsSinC(rsComplex<T> z)
  {
    rsComplex<T> tmp = rsExpC(rsComplex<T>(0.0,  1.0) * z); // tmp = exp(i*z);
    tmp        -= rsExpC(rsComplex<T>(0.0, -1.0) * z);      // tmp = exp(i*z) - exp(-i*z)
    return tmp * rsComplex<T>(0.0, -0.5);
  }

  template<class T>
  rsComplex<T> rsSinhC(rsComplex<T> z)
  {
    return 0.5 * ( expC(z) - expC(-z) );
  }

  template<class T>
  rsComplex<T> rsSnC(rsComplex<T> u, double k)
  {
    int    M = 7;     // fixed number of Landen iterations
    double v[7];      // array to store the vector of descending elliptic moduli
    rsLanden(k, M, v);

    // initialization:
    rsComplex<T> w = rsSinC(u * PI/2.0);

    // ascending Landen/Gauss transformation:
    for(int n = M-1; n >= 0; n--)
      w = (1.0+v[n])*w / (1.0+v[n]*w*w);

    return w;
  }

  /*
  template<class T>
  rsComplex<T> RSLib::rsSqrtC(rsComplex<T> z)
  {
  double r = rsSqrt(z.getRadius());
  double p = 0.5*(z.getAngle());
  double s, c;
  sinCos(p, &s, &c);
  return rsComplex<T>(r*c, r*s);            // re = r*cos(p), im = r*sin(p)
  }
  // drag the implementation for Complex.h over here
  */

  template<class T>
  rsComplex<T> rsTanC(rsComplex<T> z)
  {
    return sinC(z) / cosC(z);
  }

  template<class T>
  rsComplex<T> rsTanhC(rsComplex<T> z)
  {
    return sinhC(z) / coshC(z);
  }

  template<class T>
  int rsGetNumFiniteValues(rsComplex<T> *a, int N)
  {
    int result = 0;
    for(int n = 0; n < N; n++)
    {
      if( !a[n].isInfinite() )
        result++;
    }
    return result;
  }

  template<class T>
  int rsCopyFiniteValues(rsComplex<T> *z, rsComplex<T> *zF, int N)
  {
    int m = 0;
    for(int n = 0; n < N; n++)
    {
      if( !z[n].isInfinite() )
      {
        zF[m] = z[n];
        m++;
      }
    }
    return m;
  }

  template<class T>
  rsComplex<T> rsProductOfFiniteFactors(rsComplex<T> *a, int N)
  {
    rsComplex<T> result = rsComplex<T>(1.0, 0.0);
    for(int n = 0; n < N; n++)
    {
      if( !a[n].isInfinite() )
        result *= a[n];
    }
    return result;
  }

  template<class T>
  int rsOnlyLeftHalfPlane(rsComplex<T> *z, rsComplex<T> *zL, int N)
  {
    int m = 0;
    for(int n = 0; n < N; n++)
    {
      if( z[n].re <= 0.0 )
      {
        zL[m] = z[n];
        m++;
      }
    }
    return m;
  }

  template<class T>
  int rsOnlyUpperHalfPlane(rsComplex<T> *z, rsComplex<T> *zU, int N)
  {
    int m = 0;
    for(int n = 0; n < N; n++)
    {
      if( z[n].im >= 0.0 )
      {
        zU[m] = z[n];
        m++;
      }
    }
    return m;
  }

  template<class T>
  void rsZeroNegligibleImaginaryParts(rsComplex<T> *z, int length, double threshold)
  {
    for(int n = 0; n < length; n++)
    {
      if( fabs(z[n].im) < threshold )
        z[n].im = 0.0;
    }
  }

  template<class T>
  void rsConjugate(rsComplex<T> *z, int length)
  {
    for(int n = 0; n < length; n++)
      z[n].im = -z[n].im;
  }

  template<class T>
  bool rsComplexLessByReIm(const rsComplex<T> left, const rsComplex<T> right)
  {
    if( left.re < right.re )
      return true;
    else if( right.re < left.re )
      return false;
    else
    {
      // real parts are equal - compare by imaginary parts:
      if( left.im < right.im )
        return true;
      else if( right.re < left.re )
        return false;
      else
        return false; // both complex numbers are equal
    }
  }

  template<class T>
  void rsSortComplexArrayByReIm(rsComplex<T> *z, int length)
  {
    int i;
    std::vector<rsComplex<T> > zv;
    zv.reserve(length);
    for(i = 0; i < length; i++)
      zv.push_back(z[i]);
    std::sort(zv.begin(), zv.end(), rsComplexLessByReIm<T> );
    for(i = 0; i < length; i++)
      z[i] = zv[i];
  }

  template<class T>
  bool rsAreAllValuesReal(rsComplex<T> *z, int length, T relativeTolerance)
  {
    for(int i = 0; i < length; i++)
    {
      if( rsAbs(z[i].im) > rsAbs(relativeTolerance*z[i].re) )
        return false;
    }
    return true;
  }

}

#endif
