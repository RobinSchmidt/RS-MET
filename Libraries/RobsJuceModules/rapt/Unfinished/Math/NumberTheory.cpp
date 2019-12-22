template <class T>
T rsIntSqrt(T x)
{
  T y(1);
  for(T k = x; k != 0; k /= 4)
    y *= 2;
  T z = (y + x/y)/2;
  while( z < y )
  {
    y = z;
    z = (y + x/y)/2;
  }
  return y;
}

template<class T>
void rsFindPrimesUpTo(std::vector<T> &primes, T upperLimit)
{
  // sieve of Eratosthenes for odd numbers of the form 2*i+1:
  rsFlagArray primeFlags(upperLimit/2);
  primeFlags.setAllTrue();
  T r  = rsIntSqrt(upperLimit);
  T i  = 3;
  T i2 = i/2;
  primeFlags.setFlagFalse(0);           // "1" is not considered a prime
  while( i <= r )
  {
    if( primeFlags.isFlagTrue(i2) )
    {
      //printf("%d %s", 2*i2+1, " ");
      for(T j = (i*i)/2; j < primeFlags.getNumFlags(); j += i)
        primeFlags.setFlagFalse(j);
    }
    i2 = (T) primeFlags.getNextTrueFlag(i2+1);
    i  = 2*i2 + 1;
  }

  // collect the primes:
  T numPrimes = (T) primeFlags.getNumTrueFlags()+1; // +1 because "2" isn't included
  primes.clear();
  primes.reserve(numPrimes);
  primes.push_back(2);
  for(i = 1; i < primeFlags.getNumFlags(); i++)
  {
    if( primeFlags.isFlagTrue(i) )
      primes.push_back(2*i+1);
  }
}

template<class T>
inline rsUint32 rsFirstIndexWithMultipleOf(T p, T start)
{
  rsUint32 i = start % p;
  if( i != 0 )
    i = ((p-i) + p*(~i & 1)) / 2; // i = ((p-i) + p*isEven(i)) / 2;
  return i;
}
template<class T>
void rsGetPrimeFlags(bool *primeFlags, T start, rsUint32 numFlags, T *primeTable)
{
  memset(primeFlags, 1, numFlags*sizeof(bool));
  rsUint32 flagIndex;
  rsUint32 primeIndex = 1;
  T prime = primeTable[primeIndex];
  T end   = start + 2*(numFlags-1);
  T limit = rsIntSqrt(end);
  while( prime <= limit )
  {
    prime     = primeTable[primeIndex++];
    flagIndex = rsFirstIndexWithMultipleOf(prime, start);
    while( flagIndex < numFlags )
    {
      primeFlags[flagIndex] = 0;
      flagIndex += prime;
    }
  }
  // Remark: i tried optimizing it by pre-sieving multiples of 3 before the outer loop and then
  // run the outer loop over the remaining primes with two passes of the inner loop with a
  // stepsize of 3*prime, to avoid crossing out multiples of 3 again. This seems not to be worth
  // the trouble: it blows up the code-amount by more than a factor 2 and makes it messy and gives
  // only about 10% increase in performance.
}

template<class T>
rsUint32 rsAppendFlaggedPrimes(bool *primeFlags, T start, rsUint32 numFlags, T *primeTable,
                               rsUint32 writeStart, rsUint32 tableLength)
{
  rsUint32 np = writeStart;
  for(rsUint32 i = 0; i < numFlags && np < tableLength; i++)
  {
    if( primeFlags[i] != 0 )
      primeTable[np++] = start + 2*i;
  }
  return np - writeStart;
}

template<class T>
void rsFillPrimeTable(T *primes, rsUint32 numPrimes, rsUint32 bufferSize)
{
  if( bufferSize == 0 )
    bufferSize = 32768;
  bool *flags = new bool[bufferSize];
  primes[0] = 2;
  primes[1] = 3;
  rsUint32 numPrimesFound = 2;  // rsMax(2, partialFill)
  T        start          = 5;  // primes[numPrimesFound-1] + 2
  //rsUint32 start = 5;           // primes[numPrimesFound-1] + 2 (old, error with gcc on windows)
  rsUint32 numFlags;
  while( numPrimesFound < numPrimes )
  {
    numFlags = rsMin(bufferSize, numPrimesFound);
    rsGetPrimeFlags(flags, start, numFlags, primes);
    numPrimesFound += rsAppendFlaggedPrimes(flags, start, numFlags,
                                            primes, numPrimesFound, numPrimes);
    start += 2*numFlags;
  }
  delete[] flags;
}

// todo: make it work also for negative numbers - in this case, the first factor should be -1
template<class T>
void rsPrimeFactors(T x, std::vector<T>& factors, std::vector<T>& exponents, 
  std::vector<T> *primeTable)
{
  factors.clear();
  exponents.clear();
  if(x == 0 || x == 1) {
    factors.push_back(x);
    exponents.push_back(1);
    return;
  }

  T limit = rsIntSqrt(x);
  bool tableIsTemporary = (primeTable == nullptr);
  if( tableIsTemporary ) {
    primeTable = new std::vector<T>;

    //rsFindPrimesUpTo(*primeTable, limit);
    rsFindPrimesUpTo(*primeTable, x); 
    // this function takes the rsIntSqrt itself internally and we dont want sqrt(sqrt(x))...
    // take a closer look, make unit test - the unit tests now take MUCH longer such that it feels
    // like a hang - but i think, it takes legitimately a long time (the unit test code for these
    // factorings is currently commented out)
  }

  T i  = 0;
  T np = 0;
  T p = (*primeTable)[0];
  while( p <= limit && i < primeTable->size() ) { // 2nd condition needed to avoid access violation
    p = (*primeTable)[i];                         // for temporary tables
    if( x % p == 0 ) {
      factors.push_back(p);
      exponents.push_back(1);
      x /= p;
      while( x % p == 0 ) {
        exponents[np]++;
        x /= p;
      }
      limit = rsIntSqrt(x);
      np++;
    }
    i++;
  }
  if( x != T(1) ) {
    factors.push_back(x);
    exponents.push_back(1);
  }

  if( tableIsTemporary )
    delete primeTable;

  // \todo: maybe, when no table is passed, we should not create one but instead just use 2 and
  // every odd number as trial divisors - when a composite trial-divisor is tried, it will not be a
  // factor of the remaining part anymore anyway (for example, when 15 is tried, it can't be a
  // factor, because 3 and 5 have already been factored out)
}

template<class T>
void rsEGCD(T x, T y, T& a, T& b, T& g)
{
  T q, u, v, w, a1, b1, g1;
  bool swap = false;
  if( y > x )
  {
    rsSwap(x, y);
    swap = true;
  }
  a = T(1); b = T(0); g = x;
  u = T(0); v = T(1); w = y;
  while( w > T(0) )
  {
    q  = g/w;
    a1 = a; b1 = b; g1 = g;
    a  = u; b  = v; g  = w;
    u = a1-q*u; v = b1-q*v; w = g1-q*w;
  }
  if( swap )
    rsSwap(a, b);
}

template <class T>
T rsModularPow(const T& base, const T& modulus, rsUint64 exponent)
{
  T result = rsUnityValue(base);
  T square(base);
  while( true )
  {
    if( exponent & 1 )
      result = (result * square) % modulus;
    exponent /= 2;
    if( exponent == 0 )
      break;
    square = (square * square) % modulus;
  }
  return result;
}

template <class T>
T rsModularInverse(const T& x, const T& m)
{
  T a, b, g;
  rsEGCD(x, m, a, b, g);
  if( g == 1 )
  {
    a = a % m;    // is this required?
    if( a < 0 )
      a += m;
    return a;
  }
  else
    return T(0);
}

template <class T>
T rsPrimeModularInverse(const T& x, const T& p)
{
  return rsModularPow(x, p, p-2);
  // To make it applicable also for nonprime moduli, we could use totient(p)-1 instead of p-2 for
  // the last parameter for rsModularPow and then check, if the result times x (mod p) is equal to
  // 1 and if so, return the result, otherwise return 0. But the totient is expensive to evaluate,
  // so it seems better to use the EGCD-based algorithm in cases of nonprime moduli.
}

template <class T>
T rsPrimeModularInverse2(const T& x, const T& m)
{
  T z = x % m;
  T a(1);
  T q;
  while( z != T(1) )
  {
    if( z == T(0) )
      return z;
    q = -m/z;
    z = m + q*z;
    a = (q*a) % m;
  }
  if( a < 0 )
    a += m;
  return a;
}

template <class T>
T rsChineseRemainderWeights(T* m, T* w, rsUint32 count)
{
  //T M = rsProduct(m, count);
  T M = rsArrayTools::product(m, count);
  T Mi;
  for(rsUint32 i = 0; i < count; i++)
  {
    Mi   = M  / m[i];
    w[i] = Mi * rsModularInverse(Mi, m[i]);
  }
  return M;
}

template <class T>
T rsApplyChineseRemainderTheorem(T* r, T* w, T M, rsUint32 count)
{
  return rsArrayTools::weightedSum(w, r, count) % M;
}

template <class T>
T rsChineseRemainderTheorem(T* r, T* m, rsUint32 count)
{
  T *w = new T[count];
  T  M = rsChineseRemainderWeights(m, w, count);
  T  R = rsApplyChineseRemainderTheorem(r, w, M, count);
  delete[] w;
  return R;
}
