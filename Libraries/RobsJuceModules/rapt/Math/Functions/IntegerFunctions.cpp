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
  return n*m / rsGcd(n, m);
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