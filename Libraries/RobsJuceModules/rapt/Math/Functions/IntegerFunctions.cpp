template<class TUInt>
TUInt rsBinomialCoefficient(TUInt n, TUInt k)
{
  if(k == 0 || k == n)
    return 1;
  else if(2*k > n)
    return rsBinomialCoefficient(n, n-k);
  else
  {
    int result = n-k+1;
    for(TUInt i = 2; i <= k; i++)
    {
      result *= n-k+i;
      result /= i;
    }
    return result;
  }
}
// Maybe it could be optimized by accumulating a numerator and denominator and doing one division
// at the end. But such an algorithm would be more prone to internal overflow, I think.
// See also: https://www.youtube.com/watch?v=FpeqDFmBc2A at around 33:44


template<class TUInt>
TUInt rsBinomialCoefficientUpTo20(TUInt n, TUInt k)
{
  rsAssert(n <= 20);   // can only be used for n <= 20, otherwise internal overflow occurs
  unsigned long long nL  = n;
  unsigned long long kL  = k;
  return (TUInt)(rsProduct(kL+1ULL, nL) / rsProduct(1ULL, nL-kL));
}

template<class TUInt>
void rsBinomialDistribution(double *P, TUInt n, double p)
{
  double q = 1-p;
  std::vector<TUInt> bnk(n+1);
  rsGetLineOfPascalTriangle(&bnk[0], n);  // binomial coefficients
  for(int k = 0; k <= n; k++)
    P[k] = bnk[k] * pow(p, k) * pow(q, n-k);
}

template<class TInt>
TInt rsDelta(TInt i, TInt j)
{
  return (TInt)(i == j); // \todo: maybe inline this and templatize it on the return type
}
// todo; inline, maybe rename to rsKroneckerDelta  ..."delta" is too generic - it means different 
// things in different contexts

template<class TUInt>
TUInt rsFactorial(TUInt n)
{
  TUInt result = 1;
  for(TUInt i=1; i<=n; i++)
    result *= i;
  return result;
}

template<class TUInt>
TUInt rsGcd(TUInt m, TUInt n)
{
  TUInt lo  = rsMin(n, m);
  TUInt hi  = rsMax(n, m);
  TUInt tmp = hi;
  while(lo != 0)
  {
    tmp = hi;
    hi  = lo;
    lo  = tmp % lo;
  }
  return hi;
}
// todo: implement the gcd of an array of numbers ...is it actually true that 
// gcd(a,b,c) == gcd(gcd(a,b),c)? seems plausible because the gcd of a,b,c cannot be greater than 
// the gcd of a,b. if it's true, that relation could be used to implement the array version 
// efficiently like: g = a[0]; for(int i = 1; i < N; i++) g = gcd(g, a[i]); return g;
// -> figure out -> yes, seems to be the case:
// https://www.geeksforgeeks.org/gcd-two-array-numbers/
// https://stackoverflow.com/questions/21128981/finding-gcd-of-array-code-c-language
// Is it possible to also output the quotients of both inputs after dividing by the gcd? Often, we
// compute the gcd in order to reduce fractions in which case we then also need to divide both
// numbers by their gcd which is another two divisions. Could the algo directly output
// g = gcd(m, n) and m/g, n/g to avoid having to compute m/g, n/g by two extra divisions? Note that
// this is *not* what the extended Euclidean algo produces.

template<class TInt>
TInt rsGeneralizedDelta(TInt superscripts[], TInt subscripts[], TInt N)
{
  return rsLeviCivita(superscripts, N) * rsLeviCivita(subscripts, N);
  // is this formula right? what if leviCivita returns 0 but the subscripts
}

template<class TUInt>
TUInt rsMultinomialCoefficient(TUInt *k, TUInt kSize)
{
  // algorithm uses recursion 26.4.2 from http://dlmf.nist.gov/26.4
  TUInt sum    = k[kSize-1];
  TUInt result = 1;
  for(int j = kSize-2; j >= 0; j--)
  {
    sum    += k[j];
    result *= rsBinomialCoefficient(sum, k[j]);
  }
  return result;

  // maybe we could use a recursion for the factor rsBinomialCoefficient(sum, k[j]) as well?
  // that would make the whole algorithm O(n)
  // B(n,h) * B(n-h,k) == B(n,k) * B(n-k,h) may be useful? (formula from wikipedia)
  // ...hmm - i can't see, how this could be possible - it's probably not
}

template<class TUInt>
TUInt rsMultinomialCoefficientUpTo12(TUInt *k, TUInt kSize)
{
  // naive algorithm based directly on the definition:
  TUInt n      = rsArrayTools::sum(k, kSize);
  TUInt result = rsFactorial(n); // maybe let rsFactorial return an rsUint64, so we may have a
                                    // larger range until overflow occurs
  rsAssert(n <= 12);   // can only be used for n <= 12, otherwise internal overflow occurs
  for(TUInt i = 0; i < kSize; i++)
    result /= rsFactorial(k[i]);
  return result;
}

template<class TUInt>
void rsCreatePascalTriangle(TUInt *pt, TUInt numLines)
{
  pt[0] = 1;
  int k = 1;
  for(TUInt i = 1; i < numLines; i++)
  {
    pt[k]   = 1;
    pt[k+i] = 1;
    for(TUInt j = 1; j <= i-1; j++)
      pt[k+j] = pt[k+j-i-1] + pt[k+j-i];
    k += i+1;
  }
}

template<class TUInt>
void rsGetLineOfPascalTriangle(TUInt *c, TUInt n)
{
  rsArrayTools::fillWithZeros(c, n+1);
  c[0] = 1;
  for(TUInt i = 0; i <= n; i++)
  {
    for(TUInt j = i; j >= 1; j--)
      c[j] = c[j-1] + c[j];
  }
  // maybe this can be optimized by using the symmetry of the Pascal triangle
}

template<class T>
void rsNextPascalTriangleLine(const T* x, T* y, int N)
{
  rsAssert(N >= 0);
  T xL = T(1);
  y[0] = T(1);
  for(int i = 1; i < N; i++) 
  { 
    T xR = x[i]; 
    y[i] = xL + xR;
    xL   = xR;  
  }
  y[N] = T(1);
}
// -maybe this can be optimized using symmetry by doing something like
//  y[i] = y[i+k] = xL + xR where k depends on i and N - or maybe y[i] = y[N-i] = xL + xR?
// -maybe generalize to compute coefficients: y[i] = a*xL + b*xR...i think, this computes the coeff
//  for x^k * y^(n-k) in (a*x + b*y)^n

template<class T>
void rsPascalTriangleLine(T* y, int N)
{
  for(int n = 0; n <= N; n++) 
    rsNextPascalTriangleLine(y, y, n);
}

// see: https://en.wikipedia.org/wiki/Pascal%27s_triangle
// maybe implement also:
// https://en.wikipedia.org/wiki/Trinomial_triangle
// https://en.wikipedia.org/wiki/(2,1)-Pascal_triangle
// https://en.wikipedia.org/wiki/Bell_triangle
// https://en.wikipedia.org/wiki/Bernoulli%27s_triangle
// https://en.wikipedia.org/wiki/Leibniz_harmonic_triangle
// https://en.wikipedia.org/wiki/Eulerian_number#Basic_properties


template<class TInt>
TInt rsLeviCivita(TInt indices[], TInt N)
{
  TInt result = 1;
  for(int i = 0; i < N-1; i++)
  {
    for(int j = i+1; j < N; j++)
    {
      int d = indices[j] - indices[i];
      if(d == 0)
        return 0;
      else if(d < 0)
        result *= -1;
    }
  }
  return result;
}

template<class TUInt>
TUInt rsLcm(TUInt m, TUInt n)
{
  return n*(m/rsGcd(n, m)); // == m*(n/rsGcd(n, m)) == (n*m)/rsGcd(n, m) but the last is more
                            // prone to internal overflow
}

template<class TInt>
void rsStirlingNumbersFirstKind(TInt **s, TInt nMax)
{
  int n, k;
  for(n = 0; n <= nMax; n++)
    s[n][n] = 1;
  for(n = 1; n <= nMax; n++)
    s[n][0] = 0;
  for(n = 1; n <= nMax; n++)
  {
    for(k = 1; k < n; k++)
      s[n][k] = s[n-1][k-1] - (n-1)*s[n-1][k];
  }
}

template<class TInt>
void rsStirlingNumbersSecondKind(TInt **S, TInt nMax)
{
  int n, k;
  for(n = 0; n <= nMax; n++)
    S[n][n] = 1;
  for(n = 1; n <= nMax; n++)
    S[n][0] = 0;
  for(n = 1; n <= nMax; n++)
  {
    for(k = 1; k < n; k++)
      S[n][k] = S[n-1][k-1] + k*S[n-1][k];
  }
}

template<class TInt>
TInt rsSum(TInt min, TInt max)
{
  if(max < min)
    return 0;
  else
    return (max*(max+1) - min*(min-1)) / 2;
}
// ToDo: inline this, i.e. move to .h file
// the formula comes from using this https://de.wikipedia.org/wiki/Gau%C3%9Fsche_Summenformel for
// min and max separately, subtracting the results and simplifying
// maybe the range for which this function doesn't overflow can be extended by using one of the 
// equivalent  forms:
// H*(H+1) - L*(L-1) = H^2 + H - L^2 + L = -L^2 + H^2 + H + L = (H-L+1)*(H+L)
// ...maybe this could also make it more efficient to evaluate? -> benchmark
// https://www.wolframalpha.com/input/?i=h%5E2+%2B+h+-+l%5E2+%2B+l

/*
template<class TInt>
TInt rsWrapAround(TInt numberToWrap, TInt length)
{
  while(numberToWrap >= length)
    numberToWrap -= length;
  while(numberToWrap < 0)
    numberToWrap += length;
  return numberToWrap;
}
*/




/*

ToDo:

-Document the growth and overflow behavior of the functions to compute binomial coefficients.
 Mathematically, the asymptotic growth of binomial coeffs as function of n can be computed using
 the Stirling formula for the factorial. That might be a useful thing to know for an instantiation
 of the functions for floating point numbers. For integers, we care about the exact overflow limit 
 but with floats, we may be interested in the asymptotic growth as well. Here is some Sage code for
 this:
 # Approximates asymptotic growth of the binomial coefficient n-choose-(n/2) which is the largest
 # in the n-th line of Pascal's triangle (the one in the center) using Stirling's formula. 
 # See: http://www.hananayad.com/teaching/syde423/binomialCoefficient.pdf
 var("n")
 num = sqrt(2*pi*n)*(n/e)^n                  # approximates n!
 den = (sqrt(2*pi*(n/2))*((n/2)/e)^(n/2))^2  # approximates ((n/2)!)^2
 C  = num/den                                # approximates binomial coeff n-choose-(n/2)
 C2 = 2^n * sqrt(2*pi*n)/(pi*n)              # = C, manually simplified
 C3 = 2^n * 2 / sqrt(2*pi*n)                 # = C, according to paper
 plot([C,C2,C3], 0, 8)                       # yep, plots look the same
 C, C.simplify_full()                        # Sage doesn't really simplify it at all - WTF?
 

-Maybe implement a sqrt function for integers with semantics similar to integer division: if the 
 real result is non-integer, just round down. For example sqrt(15) = 3.8729..., so we would return
 3. 3*3 = 9, so we would get an analog of a remainder of 15-9 = 6. I don't know, if such a function
 makes sense mathematically, but we may need it to satisfy the compiler for certain class template 
 instantiations where the template was initially conceived for float/double/complex/... and uses 
 sqrt calls because (for example, in class rsPolynomial, sqrt gets called in the root-finder algo).
 ...but then, these algos may return nonsensical results, so it would actually be better, if the
 respective functions wouldn't even be compiled for an <int> instantiation. Maybe, we should drag 
 such functions out of the class rsPolynomial and make them free functions instead...or maybe 
 assemble them as static functions in some other class.


*/