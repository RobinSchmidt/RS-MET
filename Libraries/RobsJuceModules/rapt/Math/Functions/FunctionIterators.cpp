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
void rsPolynomialIterator<T, N>::setup(T* a, T h, T x0)
{
  using Poly = rsPolynomial<T>;

  y[N] = Poly::evaluate(x0, a, N); 
  T c[N+1];
  for(int m = N-1; m >= 0; m--)
  {


    for(int i = 0; i <= N-1; i++)
    {
      T bi(0);
      for(int j = i; j <= N; j++)
        bi += a[j] * pow(h, j-i) * rsBinomialCoefficient(j, j-i);
      c[i] = bi - a[i]; 
      // this is wrong! in each iteration, the a-array needs to be the c-array from the previous 
      // iteration, not the original a-array all the time
    }

    y[m] = Poly::evaluate(x0, c, m); 


  }




  int dummy = 0;

  // ToDo: optimize the call to rsBinomialCoefficient and pow
}





/*

ToDo:
-rsSineCosineIterator : public rsComplexExponentialIterator
-rsExponentialIterator, rsLinearIterator, rsQuadraticIterator, rsCubicIterator, 
 rsPolynomialIterator, rsPolyExpIterator, rsCubicExpIterator
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

*/