template<class T>
void rsLanden(T k, int M, T* v)
{
  // performs M iterations of the Landen transformation of an elliptic modulus k and returns the 
  // results in the array v which must be of length M+1
  int n;
  if(k == 0.0 || k == 1.0)
  {
    v[0] = k;
    for(n=1; n<M; n++)
      v[n] = 0.0;
    return;
  }
  else
  {
    for(n=0; n<M; n++)
    {
      k    = (k/(T(1) + sqrt(T(1)-k*k)));
      k   *= k;
      v[n] = k;
    }
  }
}

//template<class T>
//T rsIdentity(T x)
//{
//  return x;
//}

template<class T>
void rsEllipticIntegral(T k, T *K, T *Kprime, int M)
{
  T kmin = T(1e-6);
  T kmax = sqrt(1 - kmin*kmin);
  T kp, L;
  T* v  = new T[M];
  T* vp = new T[M];
  int n;

  if(k == 1.0)
    *K = RS_INF(T);
  else if(k > kmax)
  {
    kp = sqrt(1 - k*k);
    L  = -log(kp / 4);
    *K = L + (L-1)*kp*kp/4;
  }
  else
  {
    rsLanden(k, M, v);
    for(n=0; n<M; n++)
      v[n] += 1.0;
    *K = rsArrayTools::product(v, M) * T(0.5*PI);
  }

  if(k == 0.0)
    *Kprime = RS_INF(T);
  else if(k < kmin)
  {
    L       = -log(k/4);
    *Kprime = L + (L-1)*k*k/4;
  }
  else
  {
    kp = sqrt(1 - k*k);
    rsLanden(kp, M, vp);
    for(n=0; n<M; n++)
      vp[n] += 1.0;
    *Kprime = rsArrayTools::product(vp, M) * T(0.5*PI);
  }

  delete[] v;
  delete[] vp;
}
// ToDo: get rid of new/delete - use static sized stack allocated memory or, if possible, try to 
// convert everything to on-the-fly computations or let the caller pass a workspace. Also, we don't
// need different arrays for v and vp. We can just reuse v. Maybe we should warp the block from 
// rsLanden(..) to *K = into a samll lambda function, also rsLanden itself could be an internal
// lambda


// see also:
// https://github.com/jgaeddert/liquid-dsp/blob/master/src/filter/src/ellip.c

// Computes auxiliary functions f(x) and g(x) needed in the evaluation of Si(x) and Ci(x).
// This approximation follows "Handbook of mathematical functions" by Abramowitz/Stegun, 
// page 232.
template<class T>
void rsSinCosIntegralAux(T x, T* f, T* g) // usable fo x >= 1
{
  T x2 = x*x;    // x^2
  T x4 = x2*x2;  // x^4
  T x6 = x4*x2;  // x^6
  T x8 = x4*x4;  // x^8

  // Eq. 5.2.38:
  static const T
    fa1 =  38.027264, fb1 =  40.021433,
    fa2 = 265.187033, fb2 = 322.624911,
    fa3 = 335.677320, fb3 = 570.236280,
    fa4 =  38.102495, fb4 = 157.105423;
  *f = (x8 + fa1*x6 + fa2*x4 + fa3*x2 + fa4) /
    (x*(x8 + fb1*x6 + fb2*x4 + fb3*x2 + fb4));

  // Eq. 5.2.39:
  static const T
    ga1 =  42.242855, gb1 =   48.196927,
    ga2 = 302.757865, gb2 =  482.485984,
    ga3 = 352.018498, gb3 = 1114.978885,
    ga4 =  21.821899, gb4 =  449.690326;
  *g =  (x8 + ga1*x6 + ga2*x4 + ga3*x2 + ga4) /
    (x2*(x8 + gb1*x6 + gb2*x4 + gb3*x2 + gb4));
}
template<class T>
T rsSineIntegralViaAuxFunctions(T x) // usable for |x| > 1
{
  T xAbs = fabs(x);
  T s, c, f, g;
  rsSinCos(xAbs, &s, &c);
  rsSinCosIntegralAux(xAbs, &f, &g);
  return rsSign(x) * (PI/2 - f*c - g*s); // Eq. 5.2.8
}
template<class T>
T rsSineIntegralViaPowerSeries(T x)  // recommended for |x| <= 3
{
  T xAbs = fabs(x);  // |x|
  T sum  = 0.0;      // accumulator for result
  T mx2  = -x*x;     // -x^2
  int    k    = 1;        // k := 2*n+1
  T xk   = x;        // accumulator for x^k = x^(2*n+1)
  for(int n = 0; n < 15; n++)
  {
    T old = sum;
    sum += rsInverseFactorials[k] * xk / k;
    if(sum == old)
      break; // converged - maybe measure, if this really leads to efficiency improvements on 
             // the average -> throw arguments between 0 and 3 at the function, see if it performs 
             // better with or without the conditional early exit
             // also, if we know that |x| < 3, we might get away with fewer than 15 terms, 
             // in experiments, it never exceeded 13, when |x| < 3 ...but maybe we should shift
             // the crossover-point for the algorithm higher than 3
    xk *= mx2;
    k  += 2;
  }
  return sum;
}
template<class T>
T rsSineIntegralViaContinuedFractions(T x)  // recommended for |x| > 3
{
  static const int maxIterations = 100;
  T xAbs = fabs(x);
  T a;
  std::complex<T> b(1.0, xAbs);
  //std::complex<T> c = 1.0 / RS_TINY;    // gives compiler warning
  std::complex<T> c = 1.0 / 1.175494e-38; // == float minimum, fixes warning
  std::complex<T> d = 1.0 / b;
  std::complex<T> h = d;
  std::complex<T> cd;
  for(int i = 1; i <= maxIterations; i++)
  {
    a   = -i*i;
    b  += 2.0;
    d   = 1.0 / (a*d + b);
    c   = b + a/c;
    cd  = c*d;             // factor
    h  *= cd;              // product
    if(fabs(cd.real()-1.0) + fabs(cd.imag()) < RS_EPS(T))
      break;
  }
  T si, co;
  rsSinCos(xAbs, &si, &co);
  h *= std::complex<T>(co, -si);
  return rsSign(x) * (0.5*PI + h.imag());
}
template<class T>
T rsSineIntegral(T x)
{
  T xAbs = fabs(x);
  if(xAbs <= 3.0)
    return rsSineIntegralViaPowerSeries(x);
  else
  {
    return rsSineIntegralViaContinuedFractions(x); // exact but slow (really? -> measure)
    //return rsSineIntegralViaAuxFunctions(x);       // approximation
  }
}
// compare with this implementation:
// http://www.mymathlib.com/functions/sin_cos_integrals.html

template<class T>
T rsSinc(T x)
{
  if(rsAbs(x) < RS_EPS(T))
    return 1.0;
  return sin(x) / x;
}

template<class T>
T rsNormalizedSinc(T x)
{
  return rsSinc(PI*x);
}


/*

Ideas:

-rsSoftAbs(x) = x * tanh(a*x) where a controls the softness, see: 
   https://www.desmos.com/calculator/7at9bzeszu
 maybe also try some other sigmoid functions instead of tanh. But be aware that this will make 
 small signals even smaller which may be undesirable in amplitude envelope estimation. But in this
 context, the softening would be pointless anyway because the abs is typically lowpassed 
 afterwards. Maybe it's useful for a soft rectifier-like distortion. ...maybe make a FuncShaper 
 preset "SoftRectifier". With a = 1, near zero, it will look like x^2 and far from zero, it will 
 look like |x|



*/