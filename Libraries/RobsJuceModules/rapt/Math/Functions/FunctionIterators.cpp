template<class T>
rsComplexExponentialIterator<T>::rsComplexExponentialIterator(std::complex<T> a, std::complex<T> z)
{
  this->w = a;
  this->z = z;
}

template<class T>
void rsComplexExponentialIterator<T>::resetValue(std::complex<T> initialValue)
{ 
  w = initialValue;
}

template<class T>
void rsComplexExponentialIterator<T>::setZ(std::complex<T> newZ)
{
  z = newZ;
}

//=================================================================================================

/*
template<class T>
rsSineIterator<T>::rsSineIterator()
{
  a1 =  1.0806046117362795;
  s1 = -0.84147098480789650;
  s2 = -0.90929742682568171;
    // calling setup(1, 0, 1) would compute these values, but that would be more costly.
}

template<class T>
rsSineIterator<T>::rsSineIterator(T w, T p, T a)
{
  setup(w, p, a);
}
*/

template<class T>
void rsSineIterator<T>::setup(T w, T p, T a)
{
  a1 = 2.0*cos(w);
  s1 = a*sin(p-    w);
  s2 = a*sin(p-2.0*w);
  // Try to optimize computations of s1,s2 using addition theorems - i think, we may save one call
  // to sin - but maybe that should be done in a separate function "setupFast"
}

//=================================================================================================

template<class T, int N>
void rsPolynomialIterator<T, N>::setup(const T* aIn, T h, T x0)
{
  using Poly = rsPolynomial<T>;

  // Given the coefficients of a degree N polynomial in a(x) and a stepsize h, this function 
  // computes the coefficients for the degree N-1 correction polynomial c(x) such that
  // a(x) + c(x) = a(x+h). It may be used in place, i.e. c may point to the same array as a.
  auto corrector = [](const T* a, T h, int n, T* c)
  {
    for(int i = 0; i <= n-1; i++)
    {
      T bi = 0;
      for(int j = i; j <= n; j++)
        bi += a[j] * rsPow(h, j-i) * T(rsBinomialCoefficient(j, j-i));
      c[i] = bi - a[i];
    }
  };
  // maybe make this available as public static member function - oh, but maybe it should go into 
  // class rsPolynomial so it doesn't get instantiated again and again for every N
  // maybe it should be called getDifferencePolynomial or finiteDifference, forwardDifference or 
  // something ...well, there is already rsPolynomial::finiteDifference which actually already does
  // the same thing...but the implementation there sucks (allocates heap memory, etc.)...i think, 
  // it should be replaced by (a better version of) the code above

  T a[N+1];
  rsArrayTools::copy(aIn, a, N+1);
  for(int m = N; m >= 0; m--)
  {
    y[m] = Poly::evaluate(x0, a, m); 
    corrector(a, h, m, a);
  }

  // ToDo: 
  // -optimize the call to rsBinomialCoefficient and maybe also avoid pow...maybe provide a variant
  //  of the function that lets the caller pass a table of binomial coefficients and keep this 
  //  function as convenience function for prototyping (maybe add a qualilfier "Slow") to make it
  //  apparent at the call site
  //  ...or can we compute them at compile time by template trickery? see
  //  https://en.wikipedia.org/wiki/Compile-time_function_execution
  // -Maybe make a special implementation for cubic polynomials that hardcodes the required 
  //  binomial coeffs
  // -Maybe have two template parameters, one for the state y and one for the type used for a,h,x0 
  //  here. The idea is to use double precision for a,h,x0 and single for y, such that the initial
  //  state does not already have too much roundoff error. 
  //  (See Hamming: Numerical Methods..., 2nd ed., pg. 156)
}

//=================================================================================================

template<class T>
void rsSineSweepIterator<T>::setup(const rsSweepParameters<T>& p)
{
  T coeffsPhs[4];
  fitCubicWithDerivative(p.t0, p.t1, p.p0, p.p1, p.w0, p.w1, 
    &coeffsPhs[3], &coeffsPhs[2], &coeffsPhs[1], &coeffsPhs[0]);
    // ToDo: this API sucks! change it, so we can just pass the pointer to coeffsPhs

  T coeffsLogAmp[4];
  fitCubicWithDerivative(p.t0, p.t1, p.g0, p.g1, p.f0, p.f1, 
    &coeffsLogAmp[3], &coeffsLogAmp[2], &coeffsLogAmp[1], &coeffsLogAmp[0]);

  rsComplex<T> coeffs[4];
  for(int i = 0; i < 4; i++)
    coeffs[i] = rsComplex<T>(coeffsLogAmp[i], coeffsPhs[i]);

  core.setup(coeffs, T(1), p.t0);

  // Could it be numerically better, if we normalize the time interval to 0..1 before computing the
  // polynomial coeffs? -> experiments needed
}



/*

ToDo:

-For rsSineSweepIterator, use a self-written rsComplex class. std::complex doesn't seem to admit 
 usage of simd vector types for T, or in fact, any type other than float, double, long double. See
 https://en.cppreference.com/w/cpp/numeric/complex
   "The specializations std::complex<float>, std::complex<double>, and std::complex<long double> 
    are LiteralTypes for representing and manipulating complex numbers. The effect of instantiating 
    the template complex for any other type is unspecified." 
 Seems like std::complex is not as flexible i i'd like it to be. ...that applies also to 
 rsComplexExponentialIterator

-rsSineCosineIterator : public rsComplexExponentialIterator
-rsExponentialIterator, rsLinearIterator, rsQuadraticIterator, rsCubicIterator, rsCubicExpIterator
-for the polynomial iterators, see Salomon - Computer Graphics, page 275ff ("Fast Calculation of 
 the Curve" and page 698 ("Forward Differences")
-other functions: sinh, cosh, tanh, 1/cosh, 1/cosh^2 (gaussian-like?)

Ideas:

We can derive iterations for approximating certian functions y = f(x) by constructing an ODE for
the function in question and then use a numeric ODE solver to obtain successive solution values. 
Maybe a multistep method such as the Adams-Bashforth methods could be used: 
  https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods
  https://www.sciencedirect.com/science/article/pii/S0898122110007777
where the initial few values are computed directly by f(x). For example, consider
  y = x^p
Taking the derivative gives:
  y' = p*x^(p-1) = p*x^p / x = p*y/x
so the ODE is y' = p*y/x. Or let's take
  y  = 1/x = x^(-1)
  y' = -x^(-2) = - x^(-1) * x^(-1) = -y^2
so the ODE is: y' = -y^2. To actually compute y = x^p via an ODE solver, we could define z = 1/x, 
then use the z' = -z^2 equation to compute 1/x and then use that result in y' = p*y/x, so as to
avoid to compute the division by x for each iterate (we use the iteratively computed approximation
to 1/x instead of dividing by x. This is similar to what rsPowerIterator does.



Notes:
here's an interesting thread about a recursive sine oscillator:
  https://dsp.stackexchange.com/questions/124/how-to-implement-a-digital-oscillator
especially the amplitude drift compensation approach with a taylor expansion of
  1 / (sqrt(re^2 + im^2)) ~= (1/2) * (3 - (re^2 + im^2))
every 1000 (or something) samples

Here also: https://www.vicanek.de/articles/QuadOsc.pdf ..Martin says that this recursion has better
numeric stability

*/
