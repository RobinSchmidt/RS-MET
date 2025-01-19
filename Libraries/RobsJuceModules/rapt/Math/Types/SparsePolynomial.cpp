
template<class T>
void rsSparsePolynomial<T>::setupFromDenseCoeffs(const T* newCoeffs, int newNumTerms, T tol)
{
  terms.clear();
  terms.reserve(newNumTerms);
  for(int i = 0; i < newNumTerms; i++)
    if(rsAbs(newCoeffs[i]) > tol)
      terms.emplace_back(rsMonomial<T>(newCoeffs[i], i));

  //canonicalize(); // Not sure, if we should do this automatically...maybe not
  // ...wait - the result is actually ensured to be canonical already anyway.

  // It's really important to use  >  rather than  >=  in the conditional. Consider tol = 0. If we
  // would use  >=  then  >= 0  would return true when the coeff is zero, so zero coeffs would get 
  // accepted which is not what we want.
}

template<class T>
void rsSparsePolynomial<T>::addTerm(T coeff, int power, T tol)
{
  // We assume that this polynomial is in canonical representation:
  rsAssert(isCanonical());

  int i = 0;
  while(i < getNumTerms())
  {
    if(getPower(i) == power)
    {
      shiftCoeff(i, coeff);
      if(rsAbs(getCoeff(i)) <= tol)
        rsRemove(terms, (size_t) i);
      return;
    }
    else if(getPower(i) < power)
    {
      i++;
    }
    else
    {
      break;
    }
  }
  rsInsert(terms, rsMonomial<T>(coeff, power), (size_t) i);

  // After the operation, it should still be in canonical representation
  rsAssert(isCanonical());
}

template<class T>
void rsSparsePolynomial<T>::addScaled(
  const rsSparsePolynomial<T>& q, const rsMonomial<T>& s, T tol)
{
  rsAssert(rsAreAddressesDistinct(*this, q), 
           "rsSparsePolynomial::addScaled() can't be used in place.");

  for(int i = 0; i < q.getNumTerms(); i++)
    addTerm(s.getCoeff() * q.getCoeff(i), s.getPower() + q.getPower(i), tol);

  // ToDo:
  //
  // - The algorithm above that calls addTerm() in a loop may potentially trigger a lot of data 
  //   movement because each call potentially moves data. Maybe try to implement a different 
  //   algorithm that just appends the (scaled) content of q to our terms array and then calls 
  //   canonicalize(). Benchmark both variants and then choose the faster (but keep the slower 
  //   around for reference and unit tests).
}

template<class T>
void rsSparsePolynomial<T>::canonicalize(T tol)
{
  // In the empty case, we have nothing to do and we really *need* to return early in order to not 
  // get an access violation in the code below (in the  int p = getPower(0);  line):
  if(isEmpty())
    return;

  // Sort the terms by power/exponent:
  using Mon = rsMonomial<T>;
  std::sort(terms.begin(), terms.end(), 
            [](const Mon& lhs, const Mon& rhs){ return lhs.getPower() < rhs.getPower(); });

  // Consolidate multiple terms with equal power/exponent into single term: 
  int numTerms = getNumTerms();
  int p = getPower(0);              // Current power
  int r = 1;                        // Read index
  int w = 0;                        // Write index
  while(r < numTerms) 
  {
    if(getPower(r) == p)
      shiftCoeff(w, getCoeff(r));
    else 
    {
      w++;
      setTerm(w, getCoeff(r), getPower(r));
      p = getPower(r);
    }
    r++;
  }
  setNumTerms(w+1);                 // Possibly shorten the terms array
  // This algorithm works only when the terms are sorted by exponent so it doesn't really make 
  // sense to factor it out into a function in its own right. Doing so could invite calling it on 
  // unsorted term arrays in which case we would have a bug.

  // Remove terms with coefficient zero:
  rsRemoveIf(terms, [&tol](const Mon& term){ return rsAbs(term.getCoeff()) <= tol; });

  // Check postcondition:
  rsAssert(isCanonical(), "Canonicalization failed");
  // If this triggers, there's a bug in the canonicalization code above and/or in the 
  // implementation of isCanonical().


  // ToDo:
  //
  // - Maybe try using  rsHeapSort()  instead of  std::sort(). Do benchmarks with different sorting
  //   algorithms and choose the best. The sorting algorithm should have good performance 
  //   especially in the size range of a handful up to dozens or maybe hundreds. That's what we 
  //   typically deal with in digital filters which is the main intended use for this class. I'm 
  //   not sure, if we should care about the algo being a stable sort or not. But if it's not 
  //   stable, the roundoff behavior may be different and - what is more important - unpredictable.
  //   With a stable sort, the rounding that occurs in the consolidation of coeffs for terms with 
  //   equal powers, would be the same every time. That might be thing that we might want to have. 
  //   Or maybe it doesn't matter. We'll see.....
}

template<class T>
bool rsSparsePolynomial<T>::isCloseTo(const rsSparsePolynomial<T>& q, T tol) const
{
  if(getNumTerms() != q.getNumTerms())
    return false;

  for(int i = 0; i < getNumTerms(); i++)
  {
    if(getPower(i) != q.getPower(i))
      return false;
    if(rsAbs(getCoeff(i) - q.getCoeff(i)) > tol)
      return false;
  }

  return true;
}

template<class T>
int rsSparsePolynomial<T>::getMinPower() const
{
  if(isEmpty())
    return 0;
  int minPower = std::numeric_limits<int>::max();
  for(auto& term : terms)
    minPower = rsMin(minPower, term.getPower());
  return minPower;
}

template<class T>
int rsSparsePolynomial<T>::getMaxPower() const
{
  if(isEmpty())
    return 0;
  int maxPower = std::numeric_limits<int>::min();
  for(auto& term : terms)
    maxPower = rsMax(maxPower, term.getPower());
  return maxPower;

  // The implementation is written in such a way that it should still work reasonably when the
  // client code sets up terms with negative powers. The empty polynomial will still have a max
  // power (aka degree) of zero. I'm not yet sure, if we should allow for negative powers, though.
  // For causal filters, it's not needed. But maybe there could be other applications where it 
  // makes more sense. We'll see...
  //
  // Maybe init with maxPower = terms[0].getPower(). We can do this because at that point, we know
  // that the terms array is not empty. Then we can start the loop at 1, i.e. don't use a 
  // range-based loop. The range based loop should still work though - but it does one superfluous
  // iteration.
}

template<class T>
int rsSparsePolynomial<T>::getMaxPowerIndex() const
{
  rsAssert(isCanonical());
  // The output of this function is not well defined when there are multiple terms with the highest
  // power, so this function should really only be used on canonical representations. But wait:
  // In a canonical representation, the index of the maximum power is already known so we don't 
  // need to do a search in this case. It's always at terms.size()-1. Maybe we should do a test 
  // like: rsAssert(arePowersUnique()) - but such a check would be expensive (O(N^2)) on an 
  // unsorted terms array. It would even need temporary memory. On the other hand, it's only 
  // compiled into debug versions anyway.

  if(isEmpty())
    return -1;

  int maxIndex = 0;
  int maxPower = getPower(0);
  for(int i = 1; i < getNumTerms(); i++)
  {
    if(getPower(i) > maxPower)
    {
      maxPower = getPower(i);
      maxIndex = i;
    }
  }

  return maxIndex;
}

template<class T>
T rsSparsePolynomial<T>::getLeadingCoeff() const 
{ 
  return getLeadingTerm().getCoeff();
}

template<class T>
rsMonomial<T> rsSparsePolynomial<T>::getLeadingTerm() const 
{ 
  int i = getMaxPowerIndex();
  if(i != -1)
    return getTerm(i);
  else
    return rsMonomial<T>(T(0), 0);  // This branch has no test coverage yet
}

template<class T>
bool rsSparsePolynomial<T>::isCanonical(T tol) const
{
  // An empty polynomial is the canonical representation of the zero polynomial:
  if(isEmpty())
    return true;

  // Check that 0-th coeff is nonzero:
  if(rsAbs(getCoeff(0)) <= tol)
    return false;

  // Check that all other coeffs are also nonzero and that the powers are strictly increasing:
  int prevPow = getPower(0);                     // Previous power
  for(int i = 1; i < getNumTerms(); i++)
  {
    // Coeffs should be nonzero:
    if(rsAbs(getCoeff(i)) <= tol)
      return false;

    // Powers should be strictly increasing:
    int curPow = getPower(i);                    // Current power..
    if(curPow <= prevPow)
      return false;
    prevPow = curPow;                            // ..becomes previous power for next iteration.
  }

  return true;
}

template<class T>
T rsSparsePolynomial<T>::evaluateAt(T x) const 
{ 
  T y = 0;
  for(auto& term : terms)
    y += term.evaluateAt(x);
  return y;
}

template<class T>
void rsSparsePolynomial<T>::add(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();      // Number of terms in left operand p
  int Nq = q.getNumTerms();      // Number of terms in right operand q
  int Nr = Np + Nq;              // Number of terms in result r (before canonicalization)

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::subtract(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, -q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::weightedSum(
  const rsSparsePolynomial<T>& p, T wp,
  const rsSparsePolynomial<T>& q, T wq,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->setTerm(i, wp * p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->setTerm(Np + i, wq * q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::multiply(
  const rsSparsePolynomial<T>& p,
  const rsSparsePolynomial<T>& q,
  rsSparsePolynomial<T>* r, T tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np * Nq;
  r->setNumTerms(Nr);

  // Running through the loops backwards allows us to use it in place, i.e. the polynomial r can
  // point to the location of p and/or q:
  for(int i = Np-1; i >= 0; i--)
    for(int j = Nq-1; j >= 0; j--)
      r->setTerm(i*Nq+j, p.getCoeff(i) * q.getCoeff(j), p.getPower(i) + q.getPower(j));

  r->canonicalize(tol);
}

template<class T>
void rsSparsePolynomial<T>::divide(
  const rsSparsePolynomial<T>& num,
  const rsSparsePolynomial<T>& den,
  rsSparsePolynomial<T>* quot,
  rsSparsePolynomial<T>* rem,
  T tol)
{
  // Sanity checks:
  rsAssert(rsAreAddressesDistinct(num,   *quot));
  rsAssert(rsAreAddressesDistinct(num,   *rem ));
  rsAssert(rsAreAddressesDistinct(den,   *quot));
  rsAssert(rsAreAddressesDistinct(den,   *rem ));
  rsAssert(rsAreAddressesDistinct(*quot, *rem ));
  rsAssert(!den.isZero(tol));
  // What about num == den (address-wise)? I think, we should also check that this is not the case.
  // But in such a case, we can just assign quot to 1 and rem to 0 and return early. Right? Also, 
  // maybe num == rem could be ok - except for the verification of the loop invariant.

  // Initialization:
  quot->clear();             // q = 0
  rem->copyDataFrom(num);    // r = n, Invariant holds: n = d*q + r = d*0 + r = r

  // Main loop:
  while(!rem->isZero(tol) && rem->getDegree() >= den.getDegree())
  {
    rsMonomial<T> t = rem->getLeadingTerm() / den.getLeadingTerm();    // t = lead(r) / lead(d)
    quot->addTerm(t, tol);                                             // q = q + t
    rem->addScaled(den, -t, tol);                                      // r = r - t * d

    // Check the loop invariant n = d*q + r:
    rsAssert(num.isCloseTo(den * *quot + *rem, tol), "Loop invariant violated");
  }

  // The algorithm has been adapted from: 
  //
  //   https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
  //
  // In the following pseudocode, all variables (n,d,q,r,t) are polynomials (t is actually a 
  // monomial, though). Wikipedia says:
  //
  // Inputs:    n: numerator, d: denominator
  // Outputs:   q: quotient,  r: remainder
  // Require:   d != 0
  // Invariant: n = d * q + r                  # This holds at each step
  //
  // q = 0                                     # Init quotient to zero
  // r = n                                     # Init remainder to numerator
  // while( r != 0 and deg(r) >= deg(d) )
  // {
  //    t = lead(r) / lead(d)                  # t is a monomial
  //    q = q + t
  //    r = r - t * d
  // }
  // return (q, r)
  //
  //
  // ToDo:
  //
  // - Maybe at some point, when the function is battle tested well enough, we can get rid of the
  //   code that checks the loop invariant. But maybe leave it in. It helped me a lot to find a bug
  //   that I had initially in the computation of the greatest common divisor, i.e. a bug 
  //   elsewhere. It had to do with attempting to do in place processing. It would now be caught by
  //   rsAssert(rsAreAddressesDistinct(den, *rem);  which I didn't have back then.
  //
  // - Figure out and document, if it can be used in place in certain cases. If this is not the 
  //   case, explicitly document that too and maybe explain why it's not possible. In the loop, we
  //   do not seem to read from num, so maybe it's ok if num aliases to quot or rem? Maybe 
  //   num == rem is ok because rem gest initialized to num anyway. But then we will violate the 
  //   loop invariant when we allow this kind of aliasing. But that may be ok - we verify it only 
  //   for sanity checking purposes anyway. It also looks like we only ever access 
  //   den.getLeadingTerm() and never change den. That means, we could extract the leading term 
  //   (and the degree) once outside the loop and should the be free to do whatever we want with 
  //   den (or its alias) inside the loop without affecting the result. Maybe it means that rem is 
  //   allowed to alias to num and den is allowed to alias to quot after we make that change? Check
  //   this! Aaah...noo...wrong! We do read from den in rem->addScaled(dem, ...). OK - so den must
  //   be distinct. It could still make sense to drag den.getDegree() and den.getLeadingTerm out of
  //   the loop as the operations are O(N) in non-canonical representations.
}

template<class T>
void rsSparsePolynomial<T>::greatestCommonDivisorInPlace(
  rsSparsePolynomial<T>* a, 
  rsSparsePolynomial<T>* b,
  rsSparsePolynomial<T>* tmp1,
  rsSparsePolynomial<T>* tmp2,
  T tol, bool monic)
{
  while(!b->isZero(tol))
  {
    tmp1->copyDataFrom(*b);
    rsSparsePolynomial<T>::divide(*a, *tmp1, tmp2, b, tol);
    a->copyDataFrom(*tmp1);
  }
  if(monic)
    a->makeMonic();

  // Algorithm implementation has been adapted from rsRationalFunction<T>::polyGCD. 
}


/**================================================================================================


ToDo:

- Figure out what happens if client code uses negative powers. Currently, there's nothing that
  prevents this and maybe it could even make sense to allow it. But then the notion of degree
  gets murky. Maybe then there is indeed a difference between the degree and the max power in
  the case of an empty polynomial? Maybe, for the time being, we should trap attempts to set up
  terms with negative powers. This can later be relaxed, if needed.

- Maybe keep the class invariant that the polynomial is in canonical representation. 
  Implementing algorithms for both cases is a mess. Maybe prepend a __ to those member functions
  that could destroy the canonical representation to signal to the caller that they are now 
  doing something low level and potentially dangerous. Like __shiftPower(int index, int amount). 
  The regular shiftPower function can still be present. It would just call __shiftPower() and 
  then canonicalize(). Or maybe just scan through the terms to find a term with the same power
  and if one is found, consolidate the two terms into one.

- Implement a class rsSparseRationalFunction. Model it after rsRationalFunction.

- Use class rsSparseRationalFunction in rsSparseFilter (maybe as a member H). We can then 
  implement getTransferFunctionAt() as H.evaluateTyped(z)...or maybe just H(z). That would be 
  neat.

 
Notes:

- It might be tempting to write a constructor and/or setup function that takes a dense 
  polynomial, i.e. an object of type rsPolynomial<T>. But I think, that's not a good idea 
  because it would introduce unnecessary coupling.

*/