


template<class T>
bool rangeStartLess(const rsRange<T>& r1, const rsRange<T>& r2)
{
  return r1.getMin() < r2.getMin();
}

template<class T>
void splitRange(const rsRange<T>& r, rsRange<T>& rl, rsRange<T>& ru, T ratio)
{
  rl = rsRange<T>(r.getMin(),  r.getMin() + ratio*r.getSize()); // lower part of range r
  ru = rsRange<T>(rl.getMax(), r.getMax());                     // upper part of range r
}

template<class T>
void rsRatioGenerator<T>::rangeSplits(T* splitPoints, int numSplitPoints, T ratio, int strategy)
{
  // prototype implementation - may not be optimal in terms of efficiency, but is algorithmically
  // easier to understand

  int N = numSplitPoints-1;  // number of segments
  T r = ratio;
  typedef rsRange<T> Range;
  std::vector<Range> s(N);  // array/set of segments/intervals/ranges (initially empty)
  s[0] = Range(0.0, 1.0);   // seed set-of-ranges is { [0,1) }
  int n = 1;                // current number of ranges in the set
  while(n < N) 
  {
    int k = rsArray::maxIndex(&s[0], n);    // index of largest range in the current set

    // split s[k] into two ranges:
    Range rl, ru;
    if(strategy == 0) splitRange(s[k], rl, ru, r);       // always use r
    else if(strategy == 1)                               // alternating, ...
      if(rsIsOdd(n))  splitRange(s[k], rl, ru, r);       // ...odd n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...even n uses 1-r
    else if(strategy == 2)                               // alternating, ...
      if(rsIsEven(n)) splitRange(s[k], rl, ru, r);       // ...even n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...odd n uses 1-r
    // is it somehow possible to continuously morph between the odd and even version of this
    // strategy? ...maybe just "crossfade" the result split-point arrays? that would add another
    // potentially interesting dimension for tweaking ...maybe even vector-crossfade between
    // skew-right/skew-left/alternate-odd/alternate-even?

    // the lower part of s[k] replaces sk, the upper part gets appended to the array:
    s[k] = rl;
    s[n] = ru;
    n++;
  }
  // the cost of this algorithm is O(N^2) because in each iteration, we search for the maximum 
  // which is itself an O(N) operation - can we avoid this search by somehow keeping the array of
  // ranges sorted (by size), thereby turning it into an O(N) algorithm?

  // sort the array ranges by their start-point:
  rsHeapSort(&s[0], N, &rangeStartLess);
  // this is an O(N*log(N)) operation - so if we can turn the above into an O(N) operation, the 
  // overall comlexity of range splitting would still be O(N*log(N)) - certainly much better than
  // O(N^2) - but to achieve O(N), we would have to avoid the final sorting too - maybe by always
  // keeping a version sorted by size and another version sorted by start around?

  // copy range start-points into an array split-points:
  for(n = 0; n < N; n++)
    splitPoints[n] = s[n].getMin();
  splitPoints[n] = 1.0;  // ...and add the end of the last range
  // hmm..the last is always 1 and the first is always 0 - an these actually cannot be properly 
  // called split-point - maybe return only the inner, true split-points - then N = numSplitPoins+1
}


template<class T>
void rsRatioGenerator<T>::fillRatioTable(T* r, int N)
{
  typedef RatioKind RK;
  switch(kind)
  {
  case RK::metallic:      for(int i = 0; i < N; i++) r[i] = metallic(T(i));   break;
  case RK::primeSqrt:     for(int i = 0; i < N; i++) r[i] = primeSqrt(i);     break;
  case RK::primeSqrtDiff: for(int i = 0; i < N; i++) r[i] = primeSqrtDiff(i); break;

  case RK::rangeSplitSkewed: rangeSplits(r, N, p1, 0); break;
  case RK::rangeSplitOdd:    rangeSplits(r, N, p1, 1); break;
  case RK::rangeSplitEven:   rangeSplits(r, N, p1, 2); break;
  }
}


/*
Ideas:
what about ratios whose continued fraction expansion (CFE) starts with a single n and then has only
1s thereafter - they should be also "very irrational" in the sense of needing a long CFE - i think, 
more so than the other metallic ratios since what matters for the convergence of the CFE are those 
coeffs that come late in the sequence
http://www.peacefromharmony.org/docs/7-27_Stakhov_Math_of_Harmony_EN.pdf
here's a generalization witha 2nd parameter
"Further Generalization of Golden Mean in Relation to Euler's 'Divine' Equation":
https://arxiv.org/ftp/math/papers/0611/0611095.pdf



*/