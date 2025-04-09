
template<class T, class TTol>
void rsSparseRationalFunction<T, TTol>::_reduce()
{
  SparsePoly gcd = SparsePoly::greatestCommonDivisor(num, den, false);
  num = num / gcd;
  den = den / gcd;

  //rsError("Not yet correctly implemented!");

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
void rsSparseRationalFunction<T, TTol>::_canonicalize()
{
  // Maybe we should do this first because it may affect the behavior of _reduce() etc.?
  TTol tol = rsMax(num.getRoundoffTolerance(), den.getRoundoffTolerance());
  setRoundoffTolerance(tol);
  // ..and maybe we should not call setRoundoffTolerance but rather just do:
  // num.tol = tol; den.tol = tol;
  // because the setter may trigger a redundant action of removing zeros. But i think, we can't
  // access num.tol, den.tol here They are protected. But maybe the (possibly redundant) action
  // isn't a big deal.

  _reduce(); 
  _canonicalizeNumAndDen();
  _makeDenominatorMonic();
  // Maybe we should first do _canon.. and then _reduce...
}

template<class T, class TTol>
bool rsSparseRationalFunction<T, TTol>::_isCanonical() const
{
  bool can = true;
  can &= num.isCoprimeTo(den);
  can &= num._isCanonical();
  can &= den._isCanonical();
  can &= den.isMonic();
  can &= num.getRoundoffTolerance() == den.getRoundoffTolerance();
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


//-------------------------------------------------------------------------------------------------
// rsSparseDigitalTransferFunction

template<class T, class TTol>
void rsSparseDigitalTransferFunction<T, TTol>::addPreDelay(int amountInSamples)
{
  if(amountInSamples < 0)
  {
    rsError("Negative predelay is not allowed");
    return;
    // We could allow it if we already have some predelay and the given amount would just 
    // reduce it. Maybe we can relax the restriction to amountInSamples >= -num.getMinPower() or
    // something. If the minimum power is 3, we could allow a predelay amount of -3. This would 
    // then just reduce the predelay to zero.
  }

  num.shiftPowers(amountInSamples);

  // Maybe do num.shiftPowers(rsMax(amountInSamples, -getPreDelay()) ); and write unit tests for
  // this
}

template<class T, class TTol>
bool rsSparseDigitalTransferFunction<T, TTol>::_isCanonical() const
{
  bool ok = true;

  // Numerator and denominator polynomials should not be empty:
  ok &= num.getNumTerms() > 0;    // Maybe empty numerator is admissible to represent H(z) = 0? 
  ok &= den.getNumTerms() > 0;

  // We assume the num, den polynomials to be in canonical representation:
  ok &= num._isCanonical();
  ok &= den._isCanonical();

  // Filter should satisfy the a0 == 1 normalization property:
  ok &= den.getPower(0) == 0;
  ok &= den.getCoeff(0) == T(1);  // Should we use a tolerance? ...but maybe not.

  // The roundoff tolerances of numerator and denominator should match:
  ok &= num.getRoundoffTolerance() == den.getRoundoffTolerance();

  return ok;


  // ToDo:
  //
  // - Maybe we are too strict here. Maybe we should allow for empty numerators. I'm also not sure
  //   if an exact comparison to 1 is appropriate the const coeff of the denominator. Maybe we 
  //   should use an inexact comparison, i.e. use some tolerance. We'll see...
  //
  // - Maybe put the tests that access getPower(0), getCoeff(0) into a conditional to avoid access
  //   violations when the den actually is empty.
}

template<class T, class TTol>
void rsSparseDigitalTransferFunction<T, TTol>::invert()
{
  rsAssert(_isCanonical());

  // A filter with predelay cannot be inverted (at least not if it operates in realtime). The 
  // best thing we can do in this case is to invert the filter up to the predelay. Removing the
  // predelay ensures that the 0-th term in the numerator has power of 0, i.e. it's a  b0 * z^-0  
  // term:
  removePreDelay();
  rsAssert(num.getPower(0) == 0); // Numerator is already asserted to be non-empty in isCanoncial
  rsAssert(num.getCoeff(0) != 0); // ..so we can access the 0-th element without risk here

  // Swap numerator and denominator while maintaining the a0 = 0 normalization condition by 
  // appropriately scaling the numerator before and after the swap:
  T s = T(1) / num.getCoeff(0);  // Desired scaler s is 1/b0
  _scale(s);                     // Scale old numerator to achieve b0 = 1 before swap
  std::swap(num, den);           // Swap numerator and denominator. b0 is now 1 because a0 was.
  _scale(s);                     // Scale new numerator to achieve desired overall gain

  // I'm pretty sure it doesnt' allocate. The swap of the underyling std::vectors should use move 
  // semantics. Verify and document this. We need to be able to call this function on a realtime 
  // thread in some damped comb allpass filters, so it is important that this function is 
  // non-allocating.
}

template<class T, class TTol>
void rsSparseDigitalTransferFunction<T, TTol>::reflectZeros()
{
  int deg = num.getDegree();
  for(int i = 0; i < num.getNumTerms(); i++)
    num._setPower(i, deg - num.getPower(i));
  num._reverse();                               // Order array by ascending powers again

  // How about a reflectPoles() function? But that would turn stable filters into unstable ones,
  // so it's usefulness is questionable. For the time being, we can do without. Maybe we should 
  // also apply complex conjugation of the coeffs in case of complex coeffs? If so, maybe
  // do it in a function num.conjugateCoeffs() which just calls rsConj() on each coeff (which is
  // an empty function for real types)
}