#ifndef rosic_IntegerFunctions_h
#define rosic_IntegerFunctions_h

//#include "../basics/rosic_FunctionTemplates.h"

namespace rosic
{

  // todo: move some of them into the cpp file

  /** Returns one line of the Pascal triangle in the c-array, beginning with index 0. The array returned is:
  n 
  0: 1
  1: 1 1
  2: 1 2 1
  3: 1 3 3 1
  4: 1 4 6 4 1 
  etc., so the c-array should be of length n+1. The k-th element of this array is the n-over-k binomial coefficient, so the function may be 
  used to compute the binomial coefficients for all k at once, given some value of n. This is more efficient then calling 
  binomialCoefficient for each k: as a matter of fact, binomialCoefficient calls getLineOfPascalTriangle and then returns the k-th 
  element. */
  void getLineOfPascalTriangle(unsigned int *c, unsigned int n);

  /** Calculates the binomial coefficient n over k by constructing the Pascal triangle. The time-compexity of the algorithm is O(n^2) and 
  the memory-complexity is O(n). The algorithm doesn't have internal overflow issues, i.e. it returns the correct result as long as the 
  result itself does not have overflow. */
  unsigned int binomialCoefficient(unsigned int n, unsigned int k);

  /** Similar to binomialCoefficient, but uses a simpler algorithm that works only for n <= 20 because of internal overflow. The algorithm
  has O(n) time-complexity and O(1) memory-compexity. */
  unsigned int binomialCoefficientUpTo20(unsigned int n, unsigned int k);

  /** Calculates the factorial of some integer n >= 0. 
  \todo: un-inline, return unsigned long long (int64), determine maximum for "n" and state it in the comment, catch errros when it's called
  with values above maxN. */
  INLINE unsigned int factorial(unsigned int n);

  /** Calculates the greatest common divisor of n and m. */
  INLINE unsigned int gcd(unsigned int m, unsigned int n);

  /** Calculates the least common multiple of n and m. */
  INLINE unsigned int lcm(unsigned int m, unsigned int n);

  /** Power function for integers. Multiplies the base exponent-1 time with itself and returns the 
  result. */
  INLINE int powInt(int base, int exponent);

  /** Product of all integers between min and max (min and max inclusive). */
  template <class T>
  INLINE T product(T min, T max);

  /** Sum of all integers between min and max (min and max inclusive). */
  INLINE int sum(int min, int max);

  /** Wraps an integer number into the permitted range (0...length-1). */  
  INLINE int wrapAround(int numberToWrap, int length);

  //===============================================================================================
  // implementation:



  INLINE unsigned int factorial(unsigned int n)
  {
    unsigned int result = 1;
    for(unsigned int i=1; i<=n; i++)
      result *= i;
    return result;
  }

  INLINE unsigned int gcd(unsigned int m, unsigned int n)
  {
    unsigned int lo  = RAPT::rsMin(n, m);
    unsigned int hi  = RAPT::rsMax(n, m);
    unsigned int tmp = hi;
    while( lo != 0 )
    {
      tmp = hi;
      hi  = lo;
      lo  = tmp % lo;
    }
    return hi;
  }

  INLINE unsigned int lcm(unsigned int m, unsigned int n)
  {
    return n*m / gcd(n, m);
  }

  INLINE int powInt(int base, int exponent)
  {
    int result = base;
    for(int p=1; p<exponent; p++)
      result *= base;
    return result;
  }

  template <class T>
  INLINE T product(T min, T max)
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

  INLINE int sum(int min, int max)
  {
    if( max < min )
      return 0;
    else
      return ( max*(max+1) - min*(min-1) ) / 2;
  }

  INLINE int wrapAround(int numberToWrap, int length)
  {
    while( numberToWrap >= length )
      numberToWrap -= length;
    while( numberToWrap < 0 )
      numberToWrap += length;   
    return numberToWrap;
  }

} // end namespace rosic

#endif // #ifndef rosic_IntegerFunctions_h