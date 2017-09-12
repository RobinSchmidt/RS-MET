#ifndef RS_INTEGERFUNCTIONS_H
#define RS_INTEGERFUNCTIONS_H

namespace RSLib
{

  /** Calculates the binomial coefficient n over k. */
  RSLib_API unsigned int rsBinomialCoefficient(unsigned int n, unsigned int k);

  /** Similar to binomialCoefficient, but uses a simpler algorithm that works only for n <= 20
  because of internal overflow. The algorithm has O(n) time-complexity and O(1)
  memory-compexity. */
  RSLib_API unsigned int rsBinomialCoefficientUpTo20(unsigned int n, unsigned int k);

  /** Fills the array P of length n+1 with the values of the binomial distribution for a given n 
  and p. Here, n is the number of experiments ("coin tosses"), p is the probability to see a 
  particular result ("heads") in a single experiment (toss) and the array entry P[k], k = 0...n is 
  the probability to see k heads in n tosses. The probability is given by 
  P[k] = B(n,k) * p^k * (1-p)^(n-k) where B(n,k) is the binomial coefficient "n-choose-k". Works 
  only up to n = 34 - above that, internal overflow occurs. */
  RSLib_API void rsBinomialDistribution(double *P, unsigned int n, double p);

  /** Kronecker delta function. Returns 1, if i == j and 0 otherwise. */
  RSLib_API int rsDelta(int i, int j);

  /** Calculates the factorial of some integer n >= 0.
  \todo: un-inline, return unsigned long long (rsInt64), determine maximum for "n" and state it in
  the comment, catch errros when it's called with values above maxN. */
  RSLib_API unsigned int rsFactorial(unsigned int n);

  template <class T>
  RS_INLINE T rsFactorial(T n);

  /** Calculates the greatest common divisor of n and m. */
  RSLib_API unsigned int rsGcd(unsigned int m, unsigned int n);

  /** Implements the generalized Kronecker delta. This is defined as +1 if the
  superscripts-array is an even permutation of the subscripts-array, -1 if the
  superscripts-array is an odd permutation of the subscripts-array and 0 otherwise.
  @see leviCivita  */
  RSLib_API int rsGeneralizedDelta(int superscripts[], int subscripts[], int N);

  /** Computes the multinomial coefficient defined as:
  c = (      n     ) = n! / (k1! * k2! * ... * km!)
      (k1,k2,...,km)                          
  where m = kSize and n = sum(k). Scales like O(n^2) in computation and O(1) in memory usage. */
  RSLib_API rsUint32 rsMultinomialCoefficient(rsUint32 *k, rsUint32 kSize);

  /** Returns the same result as rsMultinomialCoefficient for n = sum(k) <= 12. Uses a naive 
  algorithm based directly on the definition of the multinomial coefficient. I don't really know,
  which algorithm is more efficient. \todo: measure performance. */
  RSLib_API rsUint32 rsMultinomialCoefficientUpTo12(rsUint32 *k, rsUint32 kSize);

  /** Fills the passed array with all lines of a pascal-triangle up to numLines. The minimum 
  required length for the array is given by length >= (numLines*(numLines+1))/2. To access values 
  from this array, you can use rsPascalTriangle().  */
  RSLib_API void rsCreatePascalTriangle(rsUint32 *pascalTriangle, rsUint32 numLines);

  /** Returns the k-th value of the n-th line in a pascal triangle that is supposed to have been
  created previously via rsCreatePascalTriangle. In this creation, you must make sure to create at
  least n lines. */
  RS_INLINE unsigned int rsPascalTriangle(rsUint32 *pascalTriangle, rsUint32 n, rsUint32 k)
  {
    return pascalTriangle[((n*(n+1))>>1)+k]; 
  }
  // \todo provide general functions triangleArrayRead/Write/Allocate/Free

  /** Returns one line of the Pascal triangle in the c-array which must be of length n+1, 
  beginning with index 0. The array returned is:
  n
  0:     1
  1:    1 1
  2:   1 2 1
  3:  1 3 3 1
  4: 1 4 6 4 1
  etc.
  so the c-array should be of length n+1. The k-th element of this array is the n-over-k binomial
  coefficient, so the function may be used to compute the binomial coefficients for all k at once,
  given some value of n. This is more efficient then calling binomialCoefficient for each k. 
  (\todo: is it really more efficient? measure this) */
  RSLib_API void rsGetLineOfPascalTriangle(unsigned int *c, unsigned int n);
    // maybe rename to rsBinomialCoeffs

  /** Maps an integer index in the range 0...numIndices-1 into a normalized floating point number
  in the range 0...1. */
  RS_INLINE float rsIndexToNormalizedValue(int index, int numIndices);

  /** Calculates the least common multiple of n and m. */
  RSLib_API unsigned int rsLcm(unsigned int m, unsigned int n);

  /** Levi-Civita symbol. The function assumes that the passed "indices" array has N elements and
  each element is some number from the range 1...N. The function returns 1, if the array is an even
  permutation of the numbers 1,2,3...N in their natural order where "even permutation" means that
  it can be obtained from the naturally ordered array by an even number (including zero) of swaps.
  It returns -1, if the array is an odd permutation (with likewise definition). If it is no
  permutation at all (which is the case when at least one value appears more than once in the
  array), it returns 0.
  The function may also be used the same way, when the indices array contains numbers from
  0...N-1 instead of 1...N. */
  RSLib_API int rsLeviCivita(int indices[], int N);

  /** Maps a normalized floating point number in the range 0...1 into an integer index in the range
  0...numIndices-1. */
  RS_INLINE int rsNormalizedValueToIndex(float normalizedValue, int numIndices);

  /** Power function for unsigned integers. */
  RSLib_API rsUint32 rsPowInt(rsUint32 base, rsUint32 exponent);

  /** Product of all integers between min and max (min and max inclusive). It is templated so it
  may be used for different types of integers (signed/unsigned, long/short, etc.).
  \todo: move to FunctionTemplates
  */
  template <class T>
  RS_INLINE T rsProduct(T min, T max);

  /** Fills the 2-dimensional array s (which is supposed to be of size (nMax+1) times (nMax+1)) with 
  Stirling numbers of the first kind. The array is actually triangular and looks like this:
  n |k=  0   1   2   3   4   5
  -----------------------------
  0 |    1     
  1 |    0   1    
  2 |    0  -1   1
  3 |    0   2  -3   1   
  4 |    0  -6  11  -6   1
  5 |    0  24 -50  35 -10   1
  so, for k > n, the array will not be touched. When you access the array, use n as the 1st index and 
  k as the 2nd index. The Stirling numbers of the first kind s(n,k) are the coefficients in the 
  expansion x_(n) = sum_k=0^n s(n,k) x^n, where x_(n) is the falling factorial, defined as 
  x_(n) := x*(x-1)*(x-2)*...*(x-n+1).  */
  RSLib_API void rsStirlingNumbersFirstKind(int **s, int nMax);
   // \todo use "triangular" array for s

  /** Fills the 2-dimensional array s (which is supposed to be of size (nMax+1) times (nMax+1)) with 
  Stirling numbers of the first kind. The array is actually triangular and looks like this:
  n |k=  0   1   2   3   4   5
  -----------------------------
  0 |    1     
  1 |    0   1    
  2 |    0   1   1
  3 |    0   1   3   1   
  4 |    0   1   7   6   1
  5 |    0   1  15  25  10   1
  so, for k > n, the array will not be touched. When you access the array, use n as the 1st index and 
  k as the 2nd index. The Stirling numbers of the second kind S(n,k) are the coefficients in the 
  expansion x^n = sum_k=0^n S(n,k) x_(n). @see rsStirlingNumbersFirstKind */
  RSLib_API void rsStirlingNumbersSecondKind(int **S, int nMax);

  /** Sum of all integers between min and max (min and max inclusive). */
  RSLib_API int rsSum(int min, int max);
    // \todo: if possible, make a generalized version that returns the sum of all k^n for some
    // fixed n, k = min...max (see, if wolfram alpha can solve the sum) - it says:
    // http://www.wolframalpha.com/input/?i=sum&a=*C.sum-_*Calculator.dflt-&f2=k%5En&f=Sum.sumfunction_k%5En&f3=p&f=Sum.sumlowerlimit%5Cu005fp&f4=q&f=Sum.sumupperlimit2%5Cu005fq&a=*FVarOpt.1-_**-.***Sum.sumvariable---.*--
    // sum_{k=p}^q k^n = Z(-n,p) - Z(-n,q+1) where Z is the Hurwitz-Zeta function which, for 
    // negative integer -n can be expressed in terms of Bernoulli polynomials:
    // Z(-n,x) = -B_{n+1}(x) / (n+1), according to wikipedia. 
    // http://en.wikipedia.org/wiki/Hurwitz_zeta_function So, it seems, we need to be able to
    // (construct and) evaluate Bernoulli polynomials for this.
    // http://en.wikipedia.org/wiki/Bernoulli_polynomials#Sums_of_pth_powers
    // http://en.wikipedia.org/wiki/Faulhaber%27s_formula

  /** Wraps an integer number into the permitted range (0...length-1). */
  RSLib_API int rsWrapAround(int numberToWrap, int length);


    
  // more to come: primeFactors (should return an array of pairs of BigInts with the factor itself 
  // and its exponent), etc..

  //===============================================================================================
  // implementation of inlined functions::

  template <class T>
  RS_INLINE T rsFactorial(T n)
  {
    return rsProduct(T(1), n);
  }

  RS_INLINE float rsIndexToNormalizedValue(int index, int numIndices)
  {
    return (float) (2 * index + 1) / (float) (2 * numIndices);
  }

  RS_INLINE int rsNormalizedValueToIndex(float normalizedValue, int numIndices)
  {
    return (int) floor(normalizedValue * numIndices);
  }

  template <class T>
  RS_INLINE T rsProduct(T min, T max)
  {
    if( max < min )
      return T(1);
    else
    {
      T accu = T(1);
      for(T i = min; i <= max; i++)
        accu *= i;
      return accu;
    }
  }

}

#endif
