


template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::addTerm(T coeff, int power)
{
  // We assume that this polynomial is in canonical representation:
  rsAssert(isCanonical());

  int i = 0;
  while(i < getNumTerms())
  {
    if(getPower(i) == power)
    {
      _shiftCoeff(i, coeff);
      if( rsIsNegligible(getCoeff(i), tol) )
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

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::addScaled(
  const rsSparsePolynomial<T, TTol>& q, const rsMonomial<T>& s)
{
  rsAssert(rsAreAddressesDistinct(*this, q), 
           "rsSparsePolynomial::addScaled() can't be used in place.");

  for(int i = 0; i < q.getNumTerms(); i++)
    addTerm(s.getCoeff() * q.getCoeff(i), s.getPower() + q.getPower(i));

  // ToDo:
  //
  // - The algorithm above that calls addTerm() in a loop may potentially trigger a lot of data 
  //   movement because each call potentially moves data. Maybe try to implement a different 
  //   algorithm that just appends the (scaled) content of q to our terms array and then calls 
  //   canonicalize(). Benchmark both variants and then choose the faster (but keep the slower 
  //   around for reference and unit tests). Maybe implement an _addScaled or _appendScaled()
  //   function that client code can call (perhaps in combination with canonicalize())
  //
  // - Make it work in place, i.e. when this == &q. We may need to write a special case handler
  //   for that.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::canonicalize(TTol tol)
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
      _shiftCoeff(w, getCoeff(r));
    else 
    {
      w++;
      _setTerm(w, getCoeff(r), getPower(r));
      p = getPower(r);
    }
    r++;
  }
  _setNumTerms(w+1);                // Possibly shorten the terms array
  // This algorithm works only when the terms are sorted by exponent so it doesn't really make 
  // sense to factor it out into a function in its own right. Doing so could invite calling it on 
  // unsorted term arrays in which case we would have a bug. Or maybe if we split it out, we should
  // assert that the terms are sorted.

  // Remove terms with coefficient zero:
  //rsRemoveIf(terms, [&tol](const Mon& term){ return rsAbs(term.getCoeff()) <= tol; });  // old
  rsRemoveIf(terms, [&tol](const Mon& term){ return rsIsNegligible(term.getCoeff(), tol); }); // new


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
  //   equal powers, would be the same every time. That might be a thing that we might want to 
  //   have. Or maybe it doesn't matter. We'll see.....
  //
  // - Maybe we should split the canonicalize function into sortTermsByPower(), 
  //   combineTermsWithSamePower(), removeTermsWithZeroCoeff(). The combineTermsWithSamePower 
  //   should assert that the terms array is sorted. Maybe we should have a function 
  //   areTermsSorted or arePowersAscending/arePowersStrictlyAscending
}

template<class T, class TTol>
bool rsSparsePolynomial<T, TTol>::isCloseTo(const rsSparsePolynomial<T, TTol>& q, TTol tol) const
{
  if(getNumTerms() != q.getNumTerms())
    return false;

  for(int i = 0; i < getNumTerms(); i++)
  {
    if(getPower(i) != q.getPower(i))
      return false;

    if( !rsIsNegligible(getCoeff(i) - q.getCoeff(i), tol) )
      return false;
    // ToDo: Maybe use rsIsCloseTo(getCoeff(i), q.getCoeff(i), tol). But we may need a new 
    // implementation for that - one that takes a TTol template parameter

  }

  return true;

  // Maybe assert that *this and q are canonical. Maybe we should add such assertion everywhere 
  // where we assume a canonical representation. That's quite an overhead because the test is 
  // (moderately) costly - but only in debug versions.
}

template<class T, class TTol>
int rsSparsePolynomial<T, TTol>::_getMinPower() const
{
  if(isEmpty())
    return 0;
  int minPower = std::numeric_limits<int>::max();
  for(auto& term : terms)
    minPower = rsMin(minPower, term.getPower());
  return minPower;
}

template<class T, class TTol>
int rsSparsePolynomial<T, TTol>::_getMaxPower() const
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

template<class T, class TTol>
int rsSparsePolynomial<T, TTol>::_getMaxPowerIndex() const
{
  //rsAssert(isCanonical());
  // The output of this function is not well defined when there are multiple terms with the highest
  // power, so this function should really only be used on canonical representations. But wait:
  // In a canonical representation, the index of the maximum power is already known so we don't 
  // need to do a search in this case. It's always at terms.size()-1. Maybe we should do a test 
  // like: rsAssert(arePowersUnique()) - but such a check would be expensive (O(N^2)) on an 
  // unsorted terms array. It would even need temporary memory. On the other hand, it's only 
  // compiled into debug versions anyway.
  //
  // I think, we should remove this assertion. We can make the function well defined even in case 
  // of mutliple terms with same exponent. It should just return the index of the first or last of 
  // these terms then. I think, it currently returns the first. If we would use 
  // "if(getPower(i) >= maxPower)" rather than "if(getPower(i) > maxPower)", it would return the 
  // last, I think. In an underscore-prefixed method, it is not appropriate to assume a canonical
  // representation - that's what the underscore means!


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

template<class T, class TTol>
T rsSparsePolynomial<T, TTol>::_getLeadingCoeff() const 
{ 
  return _getLeadingTerm().getCoeff();
}

template<class T, class TTol>
rsMonomial<T> rsSparsePolynomial<T, TTol>::_getLeadingTerm() const 
{ 
  int i = _getMaxPowerIndex();
  if(i != -1)
    return getTerm(i);
  else
    return rsMonomial<T>(T(0), 0);  // This branch has no test coverage yet
}

template<class T, class TTol>
bool rsSparsePolynomial<T, TTol>::isCanonical(TTol tol) const
{
  // TODO: use rsIsNegligible instead of direct comparisons with tol

  // An empty polynomial is the canonical representation of the zero polynomial:
  if(isEmpty())
    return true;

  // Check that 0-th coeff is nonzero:
  //if(rsAbs(getCoeff(0)) <= tol)         // old
  if( rsIsNegligible(getCoeff(0), tol) )  // new
    return false;

  // Check that all other coeffs are also nonzero and that the powers are strictly increasing:
  int prevPow = getPower(0);                     // Previous power
  for(int i = 1; i < getNumTerms(); i++)
  {
    // Coeffs should be nonzero:
    //if(rsAbs(getCoeff(i)) <= tol)         // old
    if( rsIsNegligible(getCoeff(i), tol) )  // new
      return false;

    // Powers should be strictly increasing:
    int curPow = getPower(i);                    // Current power..
    if(curPow <= prevPow)
      return false;
    prevPow = curPow;                            // ..becomes previous power for next iteration.
  }

  return true;
}

template<class T, class TTol>
T rsSparsePolynomial<T, TTol>::evaluateAt(T x) const 
{ 
  T y = 0;
  for(auto& term : terms)
    y += term.evaluateAt(x);
  return y;
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::add(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r, TTol tol)
{
  int Np = p.getNumTerms();      // Number of terms in left operand p
  int Nq = q.getNumTerms();      // Number of terms in right operand q
  int Nr = Np + Nq;              // Number of terms in result r (before canonicalization)

  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::subtract(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r, TTol tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, -q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::weightedSum(
  const rsSparsePolynomial<T, TTol>& p, T wp,
  const rsSparsePolynomial<T, TTol>& q, T wq,
  rsSparsePolynomial<T, TTol>* r, TTol tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, wp * p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, wq * q.getCoeff(i), q.getPower(i));

  r->canonicalize(tol);
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::multiply(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r, TTol tol)
{
  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np * Nq;
  r->_setNumTerms(Nr);

  // Running through the loops backwards allows us to use it in place, i.e. the polynomial r can
  // point to the location of p and/or q:
  for(int i = Np-1; i >= 0; i--)
    for(int j = Nq-1; j >= 0; j--)
      r->_setTerm(i*Nq+j, p.getCoeff(i) * q.getCoeff(j), p.getPower(i) + q.getPower(j));

  // We may have to re-canonicalize to combine terms with equal exponent:
  r->canonicalize(tol);
  // Maybe a full canonicalization is not needed. Maybe the first step (the sorting) is superfluous
  // if we can assume that p and q are canonical (or even just sorted)?
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::multiplyByDenseCoeffs(const T* coeffs, int numTerms, TTol tol)
{
  int Np = getNumTerms();
  int Nq = numTerms;
  int Nr = Np * Nq;
  this->_setNumTerms(Nr);

  for(int i = Np-1; i >= 0; i--)
    for(int j = Nq-1; j >= 0; j--)
      this->_setTerm(i*Nq+j, getCoeff(i) * coeffs[j], getPower(i) + j);

  this->canonicalize(tol);
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::divide(
  const rsSparsePolynomial<T, TTol>& num,
  const rsSparsePolynomial<T, TTol>& den,
  rsSparsePolynomial<T, TTol>* quot,
  rsSparsePolynomial<T, TTol>* rem,
  TTol tol)
{
  // Sanity checks:
  rsAssert(rsAreAddressesDistinct(num,   *quot));
  rsAssert(rsAreAddressesDistinct(num,   *rem ));
  rsAssert(rsAreAddressesDistinct(den,   *quot));
  rsAssert(rsAreAddressesDistinct(den,   *rem ));
  rsAssert(rsAreAddressesDistinct(*quot, *rem ));
  rsAssert(!den._isZero(tol));   // ToDo: use canonical isZero ..or maybe not
  rsAssert(num.isCanonical());
  rsAssert(den.isCanonical());
  // What about num == den (address-wise)? I think, we should also check that this is not the case.
  // But in such a case, we can just assign quot to 1 and rem to 0 and return early. Right? Also, 
  // maybe num == rem could be ok - except for the verification of the loop invariant.

  // Initialization:
  quot->clear();             // q = 0
  rem->copyDataFrom(num);    // r = n, Invariant holds: n = d*q + r = d*0 + r = r

  // Main loop:
  while(!rem->_isZero(tol) && rem->_getDegree() >= den._getDegree())  // ToDo: use canonical isZero()/getDegree()...or should we not?
  {
    rsMonomial<T> t = rem->_getLeadingTerm() / den._getLeadingTerm();  // t = lead(r) / lead(d)
    quot->addTerm(t);                                                  // q = q + t
    rem->addScaled(den, -t);                                           // r = r - t * d

    // Check the loop invariant n = d*q + r:
    SparsePoly test = den * *quot + *rem;  // For inspection in debugger
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
  //    t = lead(r) / lead(d)                  # t is a monomial, not just a coeff (!)
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
  //   rsAssert(rsAreAddressesDistinct(den, *rem);  which I didn't have back then. But the learning
  //   is that the assertion may actually catch bugs in higher level code.
  //
  // - Figure out and document, if it can be used in place in certain cases. If this is not the 
  //   case, explicitly document that too and maybe explain why it's not possible. In the loop, we
  //   do not seem to read from num, so maybe it's ok if num aliases to quot or rem? Maybe 
  //   num == rem is ok because rem gets initialized to num anyway. But then we will violate the 
  //   loop invariant when we allow this kind of aliasing. But that may be ok - we verify it only 
  //   for sanity checking purposes anyway. It also looks like we only ever access 
  //   den.getLeadingTerm() and never change den. That means, we could extract the leading term 
  //   (and the degree) once outside the loop and should the be free to do whatever we want with 
  //   den (or its alias) inside the loop without affecting the result. Maybe it means that rem is 
  //   allowed to alias to num and den is allowed to alias to quot after we make that change? Check
  //   this! Aaah...noo...wrong! We do read from den in rem->addScaled(den, ...). OK - so den must
  //   be distinct. It could still make sense to drag den.getDegree() and den.getLeadingTerm out of
  //   the loop as the operations are O(N) in non-canonical representations.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::greatestCommonDivisorInPlace(
  rsSparsePolynomial<T, TTol>* a, 
  rsSparsePolynomial<T, TTol>* b,
  rsSparsePolynomial<T, TTol>* tmp1,
  rsSparsePolynomial<T, TTol>* tmp2,
  TTol tol, bool monic)
{
  rsAssert(a->isCanonical());
  rsAssert(b->isCanonical());
  while(!b->_isZero(tol))        // ToDo: use canonical isZero
  {
    tmp1->copyDataFrom(*b);
    rsSparsePolynomial<T, TTol>::divide(*a, *tmp1, tmp2, b, tol);
    a->copyDataFrom(*tmp1);
  }
  if(monic)
    a->makeMonic();

  // Algorithm implementation has been adapted from rsRationalFunction<T>::polyGCD. I'm not sure, 
  // if we strictly require a,b to be canonical, but let's err to the conservative side. If we use
  // the non-canonical b->isZero() function, we need at least b to be canonical, I think.
}


/**================================================================================================


ToDo:

- Verify the usage pattern of the tolerance tol. I think, many member functions that receive a tol
  parameter should now not receive the tolerance as parameter anymore. If they are non-static,
  they should use the tol member. If they are static but receive at least one sparse polynomial
  as parameter, they should retriever the tolerance from there. If they receive more than one 
  polynomial as parameter, they should use the max of all tolerances. Also, the non-static 
  functions that receive an additional sparse polynomial p should use  max(this->tol, p.tol), etc.

- We also do not yet use a relative tolerance anywhere. Maybe to facilitate this, we should provide 
  a member getScaledTolerance() or getAbsoluteTolerance that returns tol * getMaxAbsCoeff() where 
  getMaxAbsCoeff() should find the maximum absolute value of all of the coeffs. Maybe that function 
  should return a value of type TTol - not of type T. Then, whenever we need to actually use the 
  tolerance, we should retrieve it by calling getScaledTolerance(). Maybe we should use a new
  function name like rsMaxNorm<TNorm>(TArg x) that returns the maximum norm of the given x, e.g.
  for complex type, it would return max(abs(re), abs(im)) ...but maybe instead of invoking abs, it 
  should actually invoke a single argument variant rsMaxNorm<TNorm>(TArg x) that we may specialize
  for float, double, complex, etc explicitly

- We may also want to implement a getter for the unscaled tolerance (maybe getTolerance() or 
  getRelativeTolerance()) and a setter. And maybe constructors that can (optionally) take the 
  tolerance to be used. Ah - and copyFataFrom should also copy the tolerance. Check, if we need to
  do this also in some copy constructors and/or assignment operators.

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

- Implement root finding/factorization. Maybe we first need evaluation of derivatives. Maybe we can
  add an optimized function for evaluating the polynomial itself along with its 1st and 2nd 
  derivative. I think, the Laguerre root finding algorithm needs these. Look up the implementation
  in rsPolynomial where we implement it for the case of dense polynomials. Maybe try to adapt the
  implementation for the sparse case.

- Implement a function that evaluates p(x) and a given number M of its derivatives. The API should
  be like  evalWithDerivatives(T x, int M, T* f)  where f is an array of M (or maybe M+1) values
  for p(x), p'(x), p''(x), ... up to the M-th derivative. 

 
Notes:

- It might be tempting to write a constructor and/or setup function that takes a dense 
  polynomial, i.e. an object of type rsPolynomial<T>. But I think, that's not a good idea 
  because it would introduce unnecessary coupling.

- When using the potentially decanonicalizing setup methods (prefixed by an underscore), there are
  3 options:

    (1) You know exactly what you are doing and that this is in fact ok, i.e. doesn't actually
        decanonicalize.

    (2) You re-canonicalize after you have finished with your operations by calling e.g.
        canonicalize().

    (3) You don't really care if the representation is canonical or not. For many purposes, a
        non-canonical representation should work just fine, although being suboptimal.



*/