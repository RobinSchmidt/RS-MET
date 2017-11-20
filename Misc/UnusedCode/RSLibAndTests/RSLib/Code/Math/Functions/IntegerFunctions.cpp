#ifndef RSLib_IntegerFunctions_cpp
#define RSLib_IntegerFunctions_cpp

namespace RSLib
{

  unsigned int rsBinomialCoefficient(unsigned int n, unsigned int k)
  {
    if( k == 0 || k == n )
      return 1;
    else if( 2*k > n )
      return rsBinomialCoefficient(n, n-k);
    else
    {
      int result = n-k+1;
      for(unsigned int i = 2; i <= k; i++)
      {
        result *= n-k+i;
        result /= i;
      }
      return result;
    }
  }

  unsigned int rsBinomialCoefficientUpTo20(unsigned int n, unsigned int k)
  {
    rsAssert( n <= 20 );   // can only be used for n <= 20, otherwise internal overflow occurs
    unsigned long long nL  = n;
    unsigned long long kL  = k;
    return (unsigned int) (rsProduct(kL+1ULL, nL) / rsProduct(1ULL, nL-kL));
  }

  void rsBinomialDistribution(double *P, unsigned int n, double p)
  {
    double q = 1-p;
    std::vector<unsigned int> bnk(n+1);
    rsGetLineOfPascalTriangle(&bnk[0], n);  // binomial coefficients
    for(int k = 0; k <= n; k++)
      P[k] = bnk[k] * pow(p, k) * pow(q, n-k);
  }

  int rsDelta(int i, int j)
  {
    return (int) (i == j); // \todo: maybe inline this and templatize it on the return type
  }

  unsigned int rsFactorial(unsigned int n)
  {
    unsigned int result = 1;
    for(unsigned int i=1; i<=n; i++)
      result *= i;
    return result;
  }

  unsigned int rsGcd(unsigned int m, unsigned int n)
  {
    unsigned int lo  = rsMin(n, m);
    unsigned int hi  = rsMax(n, m);
    unsigned int tmp = hi;
    while( lo != 0 )
    {
      tmp = hi;
      hi  = lo;
      lo  = tmp % lo;
    }
    return hi;
  }

  int rsGeneralizedDelta(int superscripts[], int subscripts[], int N)
  {
    return rsLeviCivita(superscripts, N) * rsLeviCivita(subscripts, N);
    // is this formula right? what if leviCivita returns 0 but the subscripts
  }

  rsUint32 rsMultinomialCoefficient(rsUint32 *k, rsUint32 kSize)
  {
    // algorithm uses recursion 26.4.2 from http://dlmf.nist.gov/26.4
    rsUint32 sum    = k[kSize-1];
    rsUint32 result = 1;
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

  rsUint32 rsMultinomialCoefficientUpTo12(rsUint32 *k, rsUint32 kSize)
  {
    // naive algorithm based directly on the definition:
    rsUint32 n      = rsSum(k, kSize);
    rsUint32 result = rsFactorial(n); // maybe let rsFactorial return an rsUint64, so we may have a
                                      // larger range until overflow occurs
    rsAssert( n <= 12 );   // can only be used for n <= 12, otherwise internal overflow occurs
    for(rsUint32 i = 0; i < kSize; i++)
      result /= rsFactorial(k[i]);
    return result;
  }

  void rsCreatePascalTriangle(rsUint32 *pt, rsUint32 numLines)
  {
    pt[0] = 1;
    int k = 1;
    for(rsUint32 i = 1; i < numLines; i++)
    {
      pt[k]   = 1;
      pt[k+i] = 1;
      for(rsUint32 j = 1; j <= i-1; j++) 
        pt[k+j] = pt[k+j-i-1] + pt[k+j-i];
      k += i+1;
    }
  }

  void rsGetLineOfPascalTriangle(unsigned int *c, unsigned int n)
  {
    rsFillWithZeros(c, n+1);
    c[0] = 1;
    for(unsigned int i = 0; i <= n; i++)
    {
      for(unsigned int j = i; j >= 1; j--)
        c[j] = c[j-1] + c[j];
    }
    // maybe this can be optimized by using the symmetry of the Pascal triangle
  }

  int rsLeviCivita(int indices[], int N)
  {
    int result = 1;
    for(int i = 0; i < N-1; i++)
    {
      for(int j = i+1; j < N; j++)
      {
        int d = indices[j] - indices[i];
        if( d == 0 )
          return 0;
        else if( d < 0 )
          result *= -1;
      }
    }
    return result;
  }

  unsigned int rsLcm(unsigned int m, unsigned int n)
  {
    return n*m / rsGcd(n, m);
  }

  rsUint32 rsPowInt(rsUint32 base, rsUint32 exponent)
  {
    rsUint32 result = 1;
    for(rsUint32 p = 1; p <= exponent; p++)
      result *= base;
    return result;
  }

  void rsStirlingNumbersFirstKind(int **s, int nMax)
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

  void rsStirlingNumbersSecondKind(int **S, int nMax)
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

  int rsSum(int min, int max)
  {
    if( max < min )
      return 0;
    else
      return ( max*(max+1) - min*(min-1) ) / 2;
  }

  int rsWrapAround(int numberToWrap, int length)
  {
    while( numberToWrap >= length )
      numberToWrap -= length;
    while( numberToWrap < 0 )
      numberToWrap += length;
    return numberToWrap;
  }

}

#endif
