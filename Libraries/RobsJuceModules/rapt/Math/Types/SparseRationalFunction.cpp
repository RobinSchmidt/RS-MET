
template<class T, class TTol>
void rsSparseRationalFunction<T, TTol>::_reduce()
{
  SparsePoly gcd = SparsePoly::greatestCommonDivisor(num, den, false);
  //num._divideBy(gcd); // Oh! This works only if gcd would be a momomial!
  //den._divideBy(gcd);
  // We need to call a more general division function

  rsError("Not yet correctly implemented!");

  // ToDo:
  //
  // - Maybe automatically make the denominator monic here. May the _divideBy call destroy the
  //   canonical representation of num and den? ...Figure out!
  //
  // - Does it make sense to pass false to greatestCommonDivisor as last argument? This avoids the
  //   make-it-monic step in the gcd algo. Could it be advantageous to use a monic gcd? Maybe it 
  //   would automatically ensure that our den is again monic when it was monic before? In this 
  //   case, I think, we should pass true instead.
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

template<class T, class TTol>
bool rsSparseRationalFunction<T, TTol>::_isCanonical()
{
  bool can = true;
  can &= num.isCoprimeTo(den);
  can &= num._isCanonical();
  can &= den._isCanonical();
  can &= den.isMonic();
  return can;
}

template<class T, class TTol>
void rsSparseRationalFunction<T, TTol>::weightedSum(
  const rsSparseRationalFunction<T, TTol>& p, T wp,
  const rsSparseRationalFunction<T, TTol>& q, T wq,
  rsSparseRationalFunction<T, TTol>* r, T tol)
{
  r->den = p.den * q.den;
  r->num = wp * p.num * q.den  +  wq * q.num * p.den;

  //r->_canonicalize();  // Uncomment this!

  // This can probably be optimized with respect to avoid unnecessary temporary objects and heap
  // allocations. We may also use the gcd instead of just cross-mutiplying the denominators. But 
  // whether that will be an optimization or a pessimization....well - it may depend on the inputs
  // but I guess that most of the time, it will be a pessimization.
}

template<class T, class TTol>
void rsSparseRationalFunction<T, TTol>::weightedSumDestructive(
  rsSparseRationalFunction<T, TTol>* p, T wp,
  rsSparseRationalFunction<T, TTol>* q, T wq,
  rsSparseRationalFunction<T, TTol>* r, T tol)  // Get rid of tol param!

{
  rsAssert(rsAreAddressesDistinct(*p, *q));
  //rsAssert(rsAreAddressesDistinct(*r, *p));  // We may actually allow this!
  rsAssert(rsAreAddressesDistinct(*r, *q));
  // Maybe we can relax this? It would be really nice if r could be equal to at least one of p or
  // q. Requiring p and q to be distinct is not such a big problem. I think, we can allow this, if
  // SP::weightedSum can work in place. But at the moment. I think, it can't. But maybe it can be
  // made so. Looking at the code, it seems like it could work when the result aliases to the 1st 
  // argument. Test and document this! A test indicates that this may indeed work out. Investigate
  // this further and document! We could perhaps make it work to also allow r == q by swapping p 
  // and q (and wp and wq) in this case. But what if r == p == q? ...well...in that case, we could 
  // leave the denominator of r (and p and q) alone and just multiply the numerator by the scaler 
  // (wp+wq), I think. I think, the p == q != r case could possibly also be handled - just copy p
  // or q into r and then scale by (wp+wq)

  using SP = rsSparsePolynomial<T, TTol>;

  //              arg1        arg2        result
  SP::multiply(   p->num,     q->den,     &p->num);  // Replace p->num by p->num * q->den
  SP::multiply(   q->num,     p->den,     &q->num);  // Replace q->num by q->num * p->den
  SP::weightedSum(p->num, wp, q->num, wq, &r->num);  // Establish r->num
  SP::multiply(   p->den,     q->den,     &r->den);  // Establish r->den

  //r->_canonicalize();  // Uncomment this!

  // ToDo:
  //
  // - Document exactly, how it can be used with respect to which pointers must be distinct and 
  //   which ones may alias (and to what). Document why it's called "destructive". It is because
  //   it may destroy the input parameters in the process of computing the output. It's meant to
  //   be used in place when memory usage should be optimized and the inputs become irrelevant
  //   after the computation. The regular weightedSum function allocates temporary memory due to 
  //   usage of the =,*,+ operators. This function should not allocate iff all paremeters have 
  //   allocated enough capacity to hold the (intermediate) results.
}
// Needs tests
