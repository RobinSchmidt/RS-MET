

template<class T, class TTol>
void rsSparseRationalFunction<T, TTol>::_reduce()
{
  SparsePoly gcd = SparsePoly::greatestCommonDivisor(num, den, false);
  num._divideBy(gcd);
  den._divideBy(gcd);

  // ToDo:
  //
  // - Maybe automatically make the denominator monic here. May the _divideBy call destroy the
  //   canonical representation of num and den? ...Figure out!
  //
  // - Does it make sense to pass false to greatestCommonDivisor as last argument? This avoids the
  //   make-it-monic step. Could it be advantageous to use a monic gcd? Maybe it would 
  //   automatically ensure that our den is again monic when it was monic before? In this case, I 
  //   think, we should pass true instead.
  //
  // - In rsFraction, the implementation of reduce() looks like this:
  //
  //     T gcd = rsGcd(num, den); num /= gcd; den /= gcd;
  //
  //   Maybe we can make it look the same here. Maybe then, the duplicated code can even be 
  //   factored out into a free function template rsReduceToLowestTerms(T& num, T& den).
  //
  // - Maybe make a version of this function that doesn't always allocate, i.e. uses temporary 
  //   workspace  params. It should avoid allocations as long as the temp params have enough 
  //   capacity.
}


