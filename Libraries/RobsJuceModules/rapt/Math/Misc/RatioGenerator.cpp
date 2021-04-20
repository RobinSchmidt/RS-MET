


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
    int k = rsArrayTools::maxIndex(&s[0], n);    // index of largest range in the current set

    // split s[k] into two ranges:
    Range rl, ru;
    if(strategy == 0) splitRange(s[k], rl, ru, r);        // always use r
    else if(strategy == 1) {                              // alternating, ...
      if(rsIsOdd(n)) { splitRange(s[k], rl, ru, r);    }  // ...odd n uses r
      else           { splitRange(s[k], rl, ru, 1-r);  }} // ...even n uses 1-r
    else if(strategy == 2) {                              // alternating, ...
      if(rsIsEven(n)) { splitRange(s[k], rl, ru, r);   }  // ...even n uses r
      else            { splitRange(s[k], rl, ru, 1-r); }} // ...odd n uses 1-r
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
  rsHeapSort(&s[0], N, &rangeStartLess); // maybe make sorting optional
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
  bool sorted = false;
  typedef RatioKind RK;
  switch(kind)
  {
  case RK::metallic:       { for(int i=0; i<N; i++) r[i] = metallic(T(i), p1); sorted = true; } break;
  case RK::primePower:     { for(int i=0; i<N; i++) r[i] = primePower(i);  sorted = true; } break;
  case RK::primePowerDiff: { for(int i=0; i<N; i++) r[i] = primePowerDiff(i);             } break;
  case RK::rangeSplitSkewed: { rangeSplits(r, N, p1, 0); sorted = true; } break;
  case RK::rangeSplitOdd:    { rangeSplits(r, N, p1, 1); sorted = true; } break;
  case RK::rangeSplitEven:   { rangeSplits(r, N, p1, 2); sorted = true; } break;

  case RK::linToExp:   
  { 
    for(int i = 0; i < N; i++)
    {
      T linVal = T(1) + T(i) / T(N-1);
      T expVal = exp(linVal) / exp(2);      // optimize!
      r[i]     = (T(1)-p1)*linVal + p1*expVal;
      sorted = true;
    }
    //rangeSplits(r, N, p1, 2); sorted = true; 
  } break;


  }
  if(!sorted)
    rsHeapSort(r, N); 
  rsAssert(rsIsFiniteNumbers(r, N));           // catch numeric singularities in the algos
  //rsAssert(rsIsSortedStrictlyAscending(r, N)); // catch algorithmic flaws

  // have boolean options like invert (map to 1..2 and take reciprocal), reverse 
  // (map to 0...1 and take 1-..), extract frac part
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

apply nonlinear function to contract toward the middle or spread toward the outsides, also apply 
skew - maybe like this: pre-skew -> contract -> post-skew (maybe figure out a formula that does it 
all at once - maybe some things simplify)

//rsArrayTools::transformRange(&incs[0], &incs[0], numOscs, T(0), T(1));  // ratios in 0...1



// this has been used in an early version of rsBlepOscArray:

    // linear progression of increments - this produces really bad beating:
    T minRatio = T(1) - T(0.5) * detune;
    T maxRatio = T(1) + T(0.5) * detune;
    T ratioInc = (maxRatio-minRatio) / T(numOscs-1);
    for(int i = 0; i < numOscs; i++)
    {
      incs[i] = minRatio + i * ratioInc;

      // test:
      //incs[i] = rsMetallicRatio(T(i+1));

      //incs[i] = RG::metallic(T(i));
      //incs[i] = ratioGenerator->primeSqrt(i);
      incs[i] = ratioGenerator->primeSqrtDiff(i);

      //incs[i] = sqrt( T(i+1) );

      // use 1.x where x is the fractional part:
      incs[i] = 1 + (incs[i] - floor(incs[i])); // sorted = false;
    }
    // todo: keep track of whether or not the incs-array is sorted (some algorithms produce sorted
    // arrays, others don't) and if they are not sorted, sort them afterwards

// or this:

    // experiment random incs:
    rsNoiseGenerator<T> prng;
    prng.setSeed(1);
    T randomness = 2.0;
    T minInc = T(1) - T(0.5) * randomness;
    T maxInc = T(1) + T(0.5) * randomness;
    prng.setRange(minInc, maxInc);
    for(int i = 0; i < numOscs; i++)
      incs[i] = prng.getSample();
    rsArrayTools::cumulativeSum(&incs[0], &incs[0], numOscs);
    // can give okayish results

todo: other options: linear progression of frequencies (increment reciprocals), geometric
spacing...hmm - can we somehow "invert" the generalized mean formula to obtain a generalized
spreading function - we can do it for arithmetic mean, harmonic mean and geometric mean - but
what should be in between?

see rosic::SuperOscillator::setFreq for other spacing strategies - factor out the code from 
there to make it usable here too

todo: maybe compute the mean ratio and use it to center the mean frequency/increment - but 
what mean should we use? probably the arithmetic, but the geometric or harmonic seems 
plausible as well - maybe make it a user option (use a generalized mean and let the user set 
the exponent). actually i think, the harmonic mean of the increments (corresponding to the
arithmetic mean of the frequencies) makes more sense - experimentation needed....
    
*/