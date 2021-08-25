
template<class T>
std::complex<T> rsAcdC(std::complex<T> w, T k)
{
  if(k == 1.0)
  {
    RS_DEBUG_BREAK;    // k may not be equal to 1
    return std::complex<T>(0.0, 0.0);
  }

  int M = 7;     // fixed number of Landen iterations
  T v[7];      // array to store the vector of descending elliptic moduli
  rsLanden(k, M, v);

  int    n;
  T v1;
  for(n = 0; n < M; n++)
  {
    if(n == 0)
      v1 = k;
    else
      v1 = v[n-1];
    //w = w / (T(1) + rsSqrtC(T(1) - w*w * v1*v1)) * 2.0/(1.0+v[n]);
    w = w / (T(1) + sqrt(T(1) - w*w * v1*v1)) * T(2)/(T(1)+v[n]);
  }

  std::complex<T> u;
  if(w == T(1))
    u = 0.0;
  else
    u = T(2.0/PI) * rsAcosC(w);

  T K, Kprime;
  rsEllipticIntegral(k, &K, &Kprime);
  T R = Kprime/K;

  u.real(rsSrem(u.real(), T(4)  ));
  u.imag(rsSrem(u.imag(), T(2*R)));

  return u;
}

template<class T>
std::complex<T> rsAcosC(std::complex<T> z)
{
  std::complex<T> tmp = z + std::complex<T>(T(0), T(1)) * sqrt(T(1) - z*z);
  tmp = log(tmp);
  return tmp * std::complex<T>(0.0, -1.0);
}

template<class T>
std::complex<T> rsAcoshC(std::complex<T> z)
{
  std::complex<T> tmp = z + rsSqrtC(z - 1.0) * rsSqrtC(z + 1.0);
  return log(tmp);
}

template<class T>
std::complex<T> rsAsinC(std::complex<T> z)
{
  std::complex<T> tmp = z * std::complex<T>(0.0, 1.0) + rsSqrtC(1.0 - z*z);
  tmp = log(tmp);
  return tmp * std::complex<T>(0.0, -1.0);
}

template<class T>
std::complex<T> rsAsinhC(std::complex<T> z)
{
  std::complex<T> tmp = z + rsSqrtC(z*z + 1.0);
  return log(tmp);
}

template<class T>
std::complex<T> rsAsnC(std::complex<T> w, T k)
{
  return T(1) - rsAcdC(w, k);
}

template<class T>
std::complex<T> rsAtanC(std::complex<T> z)
{
  std::complex<T> tmp = (std::complex<T>(0.0, 1.0) + z) / (std::complex<T>(0.0, 1.0) - z);
  tmp = log(tmp);
  return tmp * std::complex<T>(0.0, -0.50);
}

template<class T>
std::complex<T> rsAtanhC(std::complex<T> z)
{
  std::complex<T> tmp = (std::complex<T>(1.0, 0.0) + z) / (std::complex<T>(1.0, 0.0) - z);
  tmp = rsLogC(tmp);
  return 0.5 * tmp;
}

template<class T>
std::complex<T> rsCdC(std::complex<T> u, T k)
{
  int M = 7;     // fixed number of Landen iterations
  T v[7];      // array to store the vector of descending elliptic moduli
  rsLanden(k, M, v);

  // initialization:
  std::complex<T> w = rsCosC(u * T(PI/2.0));

  // ascending Landen/Gauss transformation:
  for(int n=M-1; n>=0; n--)
    w = (T(1)+v[n])*w / (T(1)+v[n]*w*w);

  return w;
}

template<class T>
std::complex<T> rsCosC(std::complex<T> z)
{
  std::complex<T> tmp = exp(std::complex<T>(0.0, 1.0) * z); // tmp = exp(i*z);
  tmp        += exp(std::complex<T>(0.0, -1.0)* z);      // tmp = exp(i*z) + exp(-i*z)
  return T(0.5) * tmp;
}

template<class T>
std::complex<T> rsCoshC(std::complex<T> z)
{
  return T(0.5) * (expC(z) + expC(-z));
}

/*
template<class T>
std::complex<T> rsExpC(std::complex<T> z)
{
  std::complex<T> result;
  T tmp = exp(z.real());
  rsSinCos(z.imag(), &result.imag(), &result.real());
  result.real() *= tmp;
  result.imag() *= tmp;
  return result;
}

template<class T>
std::complex<T> rsLogC(std::complex<T> z)
{
  std::complex<T> tmp;
  tmp.real() = log(z.getRadius());
  tmp.imag() = z.getAngle();
  return tmp;
}
*/


template<class T>
std::complex<T> rsPowC(std::complex<T> basis, std::complex<T> exponent)
{
  if(basis != std::complex<T>(0.0, 0.0))
    return exp(exponent * rsLogC((std::complex<T>)basis));
  else // basis == 0
  {
    if(exponent.real() == 0.0 && exponent.imag() == 0.0)
      return std::complex<T>(1.0, 0.0);
    else if(exponent.real() < 0.0 && exponent.imag() == 0.0)
      return std::complex<T>(RS_INF(T), 0.0);
    else if(exponent.real() > 0.0 && exponent.imag() == 0.0)
      return std::complex<T>(0.0, 0.0);
    else
      return exp(exponent * rsLogC((std::complex<T>)basis));
  }
}

template<class T>
std::complex<T> rsSinC(std::complex<T> z)
{
  std::complex<T> tmp = exp(std::complex<T>(0.0, 1.0) * z); // tmp = exp(i*z);
  tmp        -= exp(std::complex<T>(0.0, -1.0) * z);      // tmp = exp(i*z) - exp(-i*z)
  return tmp * std::complex<T>(0.0, -0.5);
}

template<class T>
std::complex<T> rsSinhC(std::complex<T> z)
{
  return T(0.5) * (expC(z) - expC(-z));
}

template<class T>
std::complex<T> rsSnC(std::complex<T> u, T k)
{
  int M = 7;     // fixed number of Landen iterations
  T v[7];        // array to store the vector of descending elliptic moduli
  rsLanden(k, M, v);

  // initialization:
  std::complex<T> w = rsSinC(u * T(PI/2.0));

  // ascending Landen/Gauss transformation:
  for(int n = M-1; n >= 0; n--)
    w = (T(1)+v[n])*w / (T(1)+v[n]*w*w);

  return w;
}

/*
template<class T>
std::complex<T> RSLib::rsSqrtC(std::complex<T> z)
{
double r = rsSqrt(z.getRadius());
double p = 0.5*(z.getAngle());
double s, c;
sinCos(p, &s, &c);
return std::complex<T>(r*c, r*s);            // re = r*cos(p), im = r*sin(p)
}
// drag the implementation for Complex.h over here
*/

template<class T>
std::complex<T> rsTanC(std::complex<T> z)
{
  return sinC(z) / cosC(z);
}

template<class T>
std::complex<T> rsTanhC(std::complex<T> z)
{
  return sinhC(z) / coshC(z);
}

template<class T>
int rsGetNumFiniteValues(std::complex<T> *a, int N)
{
  int result = 0;
  for(int n = 0; n < N; n++)
  {
    if(!rsIsInfinite(a[n]))
      result++;
  }
  return result;
}

template<class T>
int rsCopyFiniteValues(const std::complex<T> *z, std::complex<T> *zF, int N)
{
  int m = 0;
  for(int n = 0; n < N; n++) {
    if(!rsIsInfinite(z[n])) {
      zF[m] = z[n];
      m++; }}
  return m;
}

template<class T>
std::complex<T> rsProductOfFiniteFactors(std::complex<T> *a, int N)
{
  std::complex<T> result = std::complex<T>(1.0, 0.0);
  for(int n = 0; n < N; n++)
  {
    if(!rsIsInfinite(a[n]))
      result *= a[n];
  }
  return result;
}

template<class T>
int rsOnlyLeftHalfPlane(std::complex<T> *z, std::complex<T> *zL, int N)
{
  int m = 0;
  for(int n = 0; n < N; n++)
  {
    if(z[n].real() <= 0.0)
    {
      zL[m] = z[n];
      m++;
    }
  }
  return m;
}

template<class T>
int rsOnlyUpperHalfPlane(std::complex<T> *z, std::complex<T> *zU, int N)
{
  int m = 0;
  for(int n = 0; n < N; n++)
  {
    if(z[n].imag() >= 0.0)
    {
      zU[m] = z[n];
      m++;
    }
  }
  return m;
}

template<class T>
void rsZeroNegligibleImaginaryParts(std::complex<T> *z, int length, T threshold)
{
  for(int n = 0; n < length; n++)
  {
    if(fabs(z[n].imag()) < threshold)
      z[n].imag(0);
  }
}

template<class T>
void rsConjugate(std::complex<T> *z, int length)
{
  for(int n = 0; n < length; n++)
    z[n].imag(-z[n].imag());
}

template<class T>
bool rsComplexLessByReIm(const std::complex<T>& left, const std::complex<T>& right)
{
  if(left.real() < right.real())
    return true;
  else if(right.real() < left.real())
    return false;
  else
  {
    // real parts are equal - compare by imaginary parts:
    if(left.imag() < right.imag())
      return true;
    else if(right.imag() < left.imag())
      return false;
    else
      return false; // both complex numbers are equal
  }
}

template<class T>
bool rsComplexLessByImRe(const std::complex<T>& left, const std::complex<T>& right)
{
  if(left.imag() < right.imag())
    return true;
  else if(right.imag() < left.imag())
    return false;
  else
  {
    // imaginary parts are equal - compare by real parts:
    if(left.real() < right.real())
      return true;
    else if(right.real() < left.real())
      return false;
    else
      return false; // both complex numbers are equal
  }
}

template<class T>
void rsSortComplexArrayByReIm(std::complex<T> *z, int length)
{
  rsHeapSort(z, length, rsComplexLessByReIm<T>);
}

template<class T>
bool rsAreAllValuesReal(std::complex<T> *z, int length, T relativeTolerance)
{
  for(int i = 0; i < length; i++)
  {
    if(rsAbs(z[i].imag()) > rsAbs(relativeTolerance*z[i].real()))
      return false;
  }
  return true;
}

template<class T>
bool rsAreNeighborsConjugates(std::complex<T> *z, int length, T tol)
{
  for(int i = 0; i < length-1; i+=2)
  {
    if(z[i].real() == RS_INF(T) && z[i+1].real() == RS_INF(T))
      continue;
    if( abs(z[i].real() - z[i+1].real()) > tol )
      return false;
    if( abs(z[i].imag() + z[i+1].imag()) > tol )
      return false;
  }
  return true;
}
