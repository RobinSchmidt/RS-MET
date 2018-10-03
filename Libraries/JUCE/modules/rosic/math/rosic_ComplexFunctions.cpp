//#include <vector>
//
//#include "rosic_ComplexFunctions.h"
//using namespace rosic;

Complex rosic::acdC(Complex w, double k)
{
  if( k == 1.0 )
  {
    DEBUG_BREAK;    // k may not be equal to 1
    return Complex(0.0, 0.0);
  }

  int    M = 7;     // fixed number of Landen iterations
  double v[7];      // array to store the vector of descending elliptic moduli
  landen(k, M, v);

  int    n;
  double v1;
  for(n=0; n<M; n++)
  {
    if( n == 0 )
      v1 = k;
    else
      v1 = v[n-1];
    w = w / (1.0 + sqrtC(1.0 - w*w * v1*v1)) * 2.0/(1.0+v[n]);
  }

  Complex u;
  if( w == 1.0 )
    u = 0.0;
  else
    u = (2.0/PI) * acosC(w);

  double K, Kprime;
  ellipticIntegral(k, &K, &Kprime);
  double R = Kprime/K;

  u.re = srem(u.re, 4.0);
  u.im = srem(u.im, 2.0*R);

  return u;
}

Complex rosic::acosC(Complex z)
{
  Complex tmp = z + Complex(0.0, 1.0) * sqrtC( 1.0 - z*z );
  tmp = logC(tmp);
  return tmp * Complex(0.0, -1.0);
}

Complex rosic::acoshC(Complex z)
{
  Complex tmp = z + sqrtC( z - 1.0 ) * sqrtC( z + 1.0 );
  return logC(tmp);
}

Complex rosic::asinC(Complex z)
{
  Complex tmp = z * Complex(0.0, 1.0) + sqrtC( 1.0 - z*z );
  tmp = logC(tmp);
  return tmp * Complex(0.0, -1.0);
}

Complex rosic::asinhC(Complex z)
{
  Complex tmp = z + sqrtC( z*z + 1.0 ) ;
  return logC(tmp);
}

Complex rosic::asnC(Complex w, double k)
{
  return 1.0 - acdC(w, k);
}

Complex rosic::atanC(Complex z)
{
  Complex tmp = ( Complex(0.0, 1.0 ) + z) / ( Complex(0.0, 1.0) - z );
  tmp = logC(tmp);
  return tmp * Complex(0.0, -0.50);
}

Complex rosic::atanhC(Complex z)
{
  Complex tmp = ( Complex(1.0, 0.0 ) + z ) / ( Complex(1.0, 0.0) - z );
  tmp = logC(tmp);
  return 0.5 * tmp;
}

Complex rosic::cdC(Complex u, double k)
{
  int    M = 7;     // fixed number of Landen iterations
  double v[7];      // array to store the vector of descending elliptic moduli
  landen(k, M, v);

  // initialization:
  Complex w = cosC(u * PI/2.0);

  // ascending Landen/Gauss transformation:
  for(int n=M-1; n>=0; n--)
    w = (1.0+v[n])*w / (1+v[n]*w*w);

  return w;
}

Complex rosic::cosC(Complex z)
{
  Complex tmp = expC(Complex(0.0, 1.0) * z); // tmp = exp(i*z);
  tmp        += expC(Complex(0.0, -1.0)* z); // tmp = exp(i*z) + exp(-i*z)     
  return 0.5 * tmp; 
}

Complex rosic::coshC(Complex z)
{  
  return 0.5 * ( expC(z) + expC(-z) ); 
}

Complex rosic::expC(Complex z)
{
  Complex result;
  double tmp = exp(z.re);
  RAPT::rsSinCos(z.im, &result.im, &result.re);
  result.re *= tmp;
  result.im *= tmp;
  return result;
  //return Complex(tmp*cos(z.im), tmp*sin(z.im));
}

Complex rosic::logC(Complex z)
{
  Complex tmp; 
  tmp.re = log(z.getRadius());
  tmp.im = z.getAngle();
  return tmp;
}

Complex rosic::powC(Complex basis, Complex exponent)
{
  if( basis != Complex(0.0, 0.0) )
    return expC( exponent * logC((Complex)basis) );
  else // basis == 0
  {
    if( exponent.re == 0.0 && exponent.im == 0.0 )
      return Complex(1.0, 0.0);
    else if( exponent.re < 0.0 && exponent.im == 0.0 )
      return Complex(INF, 0.0);
    else if( exponent.re > 0.0 && exponent.im == 0.0 )
      return Complex(0.0, 0.0);
    else
      return expC( exponent * logC((Complex)basis) );
  }
}

Complex rosic::sinC(Complex z)
{
  Complex tmp = expC(Complex(0.0,  1.0) * z); // tmp = exp(i*z);
  tmp        -= expC(Complex(0.0, -1.0) * z); // tmp = exp(i*z) - exp(-i*z)
  return tmp * Complex(0.0, -0.5); 
}

Complex rosic::sinhC(Complex z)
{  
  return 0.5 * ( expC(z) - expC(-z) ); 
}

Complex rosic::snC(Complex u, double k)
{
  int    M = 7;     // fixed number of Landen iterations
  double v[7];      // array to store the vector of descending elliptic moduli
  landen(k, M, v);

  // initialization:
  Complex w = sinC(u * PI/2.0);

  // ascending Landen/Gauss transformation:
  for(int n=M-1; n>=0; n--)
    w = (1.0+v[n])*w / (1+v[n]*w*w);

  return w;
}

Complex rosic::sqrtC(Complex z)
{
  double r = sqrt(z.getRadius());
  double p = 0.5*(z.getAngle());  
  double s, c;
  RAPT::rsSinCos(p, &s, &c);
  return Complex(r*c, r*s);            // re = r*cos(p), im = r*sin(p)
}

Complex rosic::tanC(Complex z)
{
  return sinC(z) / cosC(z);
}

Complex rosic::tanhC(Complex z)
{
  return sinhC(z) / coshC(z);
}

int rosic::getNumFiniteValues(Complex *a, int N)
{
  int result = 0;
  for(int n = 0; n < N; n++)
  {
    if( !a[n].isInfinite() )
      result++;
  }
  return result;
}

int rosic::copyFiniteValues(Complex *z, Complex *zF, int N)
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

Complex rosic::productOfFiniteFactors(Complex *a, int N)
{
  Complex result = Complex(1.0, 0.0);
  for(int n = 0; n < N; n++)
  {
    if( !a[n].isInfinite() )
      result *= a[n];
  }
  return result;
}

int rosic::onlyLeftHalfPlane(Complex *z, Complex *zL, int N)
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

int rosic::onlyUpperHalfPlane(Complex *z, Complex *zU, int N)
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

void rosic::zeroNegligibleImaginaryParts(Complex *z, int length, double threshold)
{
  for(int n = 0; n < length; n++)
  {
    if( fabs(z[n].im) < threshold )
      z[n].im = 0.0;
  }
}

bool rosic::complexLessByReIm(const Complex left, const Complex right)
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

void rosic::sortComplexArrayByReIm(Complex *z, int length)
{
  int i;
  std::vector<Complex> zv;
  zv.reserve(length);
  for(i = 0; i < length; i++)
    zv.push_back(z[i]);
  std::sort(zv.begin(), zv.end(), complexLessByReIm);
  for(i = 0; i < length; i++)
    z[i] = zv[i];
}
