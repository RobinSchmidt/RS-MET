
template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::setupFromDenseCoeffs(const T* newCoeffs, int newNumTerms)
{
  terms.clear();
  terms.reserve(newNumTerms);
  for(int i = 0; i < newNumTerms; i++)
    if( !rsIsNegligible(newCoeffs[i], tol) ) 
      terms.emplace_back(rsMonomial<T>(newCoeffs[i], i));

  rsAssert(_isCanonical()); // Dense coeffs are naturally ordered and we skip the zeros.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::addTerm(T coeff, int power)
{
  // We assume that this polynomial is in canonical representation:
  rsAssert(_isCanonical());

  // Find the point where we have to insert the term or update the coeff. In the latter case, do 
  // the update and return early
  int i = 0;
  while(i < getNumTerms())
  {
    if(getPower(i) == power)
    {
      _shiftCoeff(i, coeff);
      if( rsIsNegligible(getCoeff(i), tol) )
        rsRemove(terms, (size_t) i);
      rsAssert(_isCanonical());               // Make sure we didn't mess up canonicalness
      return;
      // Check, if this has test coverage!
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
  // Maybe move this into the else branch instead of breaking there. Then return directly from this 
  // branch. I think, this makes the code more readable. ...but no! That would be wrong because the
  // insertion may also happen at the end. Make sure to have unit test coverage for all 3 possible 
  // cases: insert, remove, append. Make sure that one of the unit tests would fail, if we would 
  // put the rsInsert() into the else branch.

  // After the operation, it should still be in canonical representation:
  rsAssert(_isCanonical());

  // Notes:
  //
  // - The assertion at the bottom is not always reached because we have this early return 
  //   statement in the if statement. 
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::addScaled(
  const rsSparsePolynomial<T, TTol>& q, const rsMonomial<T>& s)
{
  rsAssert(rsAreAddressesDistinct(*this, q), 
           "rsSparsePolynomial::addScaled() can't be used in place.");

  tol = rsMax(tol, q.tol);
  for(int i = 0; i < q.getNumTerms(); i++)
    addTerm(s.getCoeff() * q.getCoeff(i), s.getPower() + q.getPower(i));

  // ToDo:
  //
  // - The algorithm above that calls addTerm() in a loop may potentially trigger a lot of data 
  //   movement because each call potentially moves data. Maybe try to implement a different 
  //   algorithm that just appends the (scaled) content of q to our terms array and then calls 
  //   canonicalize(). Benchmark both variants and then choose the faster (but keep the slower 
  //   around for reference and unit tests). Maybe implement an _addScaled or _appendScaled()
  //   function that client code can call (perhaps in combination with canonicalize()). This may 
  //   result in less data movement - but it may blow up the required memory temporarily. So: 
  //   no - let's not do that in general. It may even lead to allocations when we really don't 
  //   want them.
  //
  // - Make it work in place, i.e. when this == &q. We may need to write a special case handler
  //   for that.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::addScaled(
  const rsSparsePolynomial<T, TTol>&  p, T scaler)
{
  tol = rsMax(tol, p.tol);
  for(int i = 0; i < p.getNumTerms(); i++)
    addTerm(scaler * p.getCoeff(i), p.getPower(i));
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::_canonicalize()
{
  // In the empty case, we have nothing to do and we really *need* to return early in order to not 
  // get an access violation in the code below (in the  int p = getPower(0);  line):
  if(terms.empty())
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
  rsRemoveIf(terms, [this](const Mon& term){ return rsIsNegligible(term.getCoeff(), tol); });
  // Factor out inot function so we can call it from setRoundoffTolerance(), too

  // Check postcondition:
  rsAssert(_isCanonical(), "Canonicalization failed");
  // If this triggers, there's a bug in the canonicalization code above and/or in the 
  // implementation of _isCanonical().


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
  rsAssert(_isCanonical() && q._isCanonical());

  if(getNumTerms() != q.getNumTerms())
    return false;

  for(int i = 0; i < getNumTerms(); i++)
  {
    if(getPower(i) != q.getPower(i))
      return false;
    if( !rsIsNegligible(getCoeff(i) - q.getCoeff(i), tol) )
      return false;
  }

  return true;

  // Maybe assert that *this and q are canonical. Maybe we should add such assertion everywhere 
  // where we assume a canonical representation. That's quite an overhead because the test is 
  // (moderately) costly - but only in debug versions.
}

template<class T, class TTol>
bool rsSparsePolynomial<T, TTol>::_isCanonical() const
{
  // An empty polynomial is the canonical representation of the zero polynomial:
  if(terms.empty())
    return true;

  // Check that 0-th coeff is nonzero:
  if( rsIsNegligible(getCoeff(0), tol) )
    return false;

  // Check that all other coeffs are also nonzero and that the powers are strictly increasing:
  int prevPow = getPower(0);                     // Previous power
  for(int i = 1; i < getNumTerms(); i++)
  {
    // Coeffs should be nonzero:
    if( rsIsNegligible(getCoeff(i), tol) )
      return false;

    // Powers should be strictly increasing:
    int curPow = getPower(i);                    // Current power..
    if(curPow <= prevPow)
      return false;
    prevPow = curPow;                            // ..becomes previous power for next iteration.
  }

  // Check that all the powers are nonnegative. Maybe this restriction can be lifted later but
  // for the time being, let's be conservative:
  for(int i = 0; i < getNumTerms(); i++)
    if(getPower(i) < 0)
      return false;

  return true;
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::add(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r)
{
  rsAssert(p._isCanonical());
  rsAssert(q._isCanonical());

  int Np = p.getNumTerms();      // Number of terms in left operand p
  int Nq = q.getNumTerms();      // Number of terms in right operand q
  int Nr = Np + Nq;              // Number of terms in result r (before canonicalization)

  r->tol = rsMax(p.tol, q.tol);
  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, q.getCoeff(i), q.getPower(i));

  r->_canonicalize();
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::subtract(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r)
{
  rsAssert(p._isCanonical());
  rsAssert(q._isCanonical());

  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->tol = rsMax(p.tol, q.tol);
  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, -q.getCoeff(i), q.getPower(i));

  r->_canonicalize();
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::weightedSum(
  const rsSparsePolynomial<T, TTol>& p, T wp,
  const rsSparsePolynomial<T, TTol>& q, T wq,
  rsSparsePolynomial<T, TTol>* r)
{
  rsAssert(p._isCanonical());
  rsAssert(q._isCanonical());

  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np + Nq;

  r->tol = rsMax(p.tol, q.tol);
  r->_setNumTerms(Nr);
  for(int i = 0; i < Np; i++)
    r->_setTerm(i, wp * p.getCoeff(i), p.getPower(i));
  for(int i = 0; i < Nq; i++)
    r->_setTerm(Np + i, wq * q.getCoeff(i), q.getPower(i));

  r->_canonicalize();
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::multiply(
  const rsSparsePolynomial<T, TTol>& p,
  const rsSparsePolynomial<T, TTol>& q,
  rsSparsePolynomial<T, TTol>* r)
{
  rsAssert(p._isCanonical());
  rsAssert(q._isCanonical());

  int Np = p.getNumTerms();
  int Nq = q.getNumTerms();
  int Nr = Np * Nq;

  r->tol = rsMax(p.tol, q.tol);
  r->_setNumTerms(Nr);

  // Running through the loops backwards allows us to use it in place, i.e. the polynomial r can
  // point to the location of p and/or q:
  for(int i = Np-1; i >= 0; i--)
    for(int j = Nq-1; j >= 0; j--)
      r->_setTerm(i*Nq+j, p.getCoeff(i) * q.getCoeff(j), p.getPower(i) + q.getPower(j));

  // We may have to re-canonicalize to combine terms with equal exponent:
  r->_canonicalize();
  // Maybe a full canonicalization is not needed. Maybe the first step (the sorting) is superfluous
  // if we can assume that p and q are canonical (or even just sorted)? Maybe factor out the 
  // partial steps of the canonicalization and think about, if we can get away with less steps 
  // here. But the create thorough unit tests for that.
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::multiplyByDenseCoeffs(const T* coeffs, int numTerms)
{
  rsAssert(_isCanonical());

  int Np = getNumTerms();
  int Nq = numTerms;
  int Nr = Np * Nq;

  _setNumTerms(Nr);

  for(int i = Np-1; i >= 0; i--)
    for(int j = Nq-1; j >= 0; j--)
      _setTerm(i*Nq+j, getCoeff(i) * coeffs[j], getPower(i) + j);

  _canonicalize();
  // Do we need this? If so, document why. I think in the loop above, it will typically happen that
  // we produce multiple terms with the same power. "getPower(i) + j" will take on the same value
  // multiple times, so _setTerm(..) will set multiple different terms to the same power. Verify 
  // this!
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::divide(
  const rsSparsePolynomial<T, TTol>& num,
  const rsSparsePolynomial<T, TTol>& den,
  rsSparsePolynomial<T, TTol>* quot,
  rsSparsePolynomial<T, TTol>* rem)
{
  // Sanity checks:
  rsAssert(num._isCanonical());
  rsAssert(den._isCanonical());
  rsAssert(!den.isZero());
  rsAssert(rsAreAddressesDistinct(num,   *quot));
  rsAssert(rsAreAddressesDistinct(num,   *rem ));
  rsAssert(rsAreAddressesDistinct(den,   *quot));
  rsAssert(rsAreAddressesDistinct(den,   *rem ));
  rsAssert(rsAreAddressesDistinct(*quot, *rem ));
  // What about num == den (address-wise)? I think, we should also check that this is not the case.
  // But in such a case, we can just assign quot to 1 and rem to 0 and return early. Right? Also, 
  // maybe num == rem could be ok - except for the verification of the loop invariant.

  // Initialization:
  TTol newTol = rsMax(num.tol, den.tol);  // Tolerance of the results
  quot->clear();                          // q = 0. Quotient is empty/zero.
  quot->tol = newTol;                     // Set up tolerance of quotient.
  *rem = num;                             // r = n. Invariant holds: n = d*q + r = d*0 + r = r
  rem->tol = newTol;                      // Important to do this after *rem = num
  
  // Main loop:
  while(!rem->isZero() && rem->getDegree() >= den.getDegree())
  {
    rsMonomial<T> t = rem->getLeadingTerm() / den.getLeadingTerm();  // t = lead(r) / lead(d)
    quot->addTerm(t);                                                // q = q + t
    rem->addScaled(den, -t);                                         // r = r - t * d

    // Check sanity and loop invariant n = d*q + r:
    rsAssert(quot->_isCanonical());
    rsAssert( rem->_isCanonical());
    rsAssert(num.isCloseTo(den * *quot + *rem, num.tol), "Loop invariant violated");
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
  // - Maybe call den.getDegree() and den.getLeadingTerm() outside the loop and use variables like
  //   denDeg, denLead in the loop. These do not change during the loop. The calls are cheap, but 
  //   still.
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
  //   the loop.
}

template<class T, class TTol>
rsSparsePolynomial<T, TTol> rsSparsePolynomial<T, TTol>::greatestCommonDivisor(
  const rsSparsePolynomial<T, TTol>& p, const rsSparsePolynomial<T, TTol>& q, bool monic)
{
  SparsePoly a = p, b = q, tmp1, tmp2;
  greatestCommonDivisorInPlace(&a, &b, &tmp1, &tmp2, monic);
  return a;
}

template<class T, class TTol>
void rsSparsePolynomial<T, TTol>::greatestCommonDivisorInPlace(
  rsSparsePolynomial<T, TTol>* a,    rsSparsePolynomial<T, TTol>* b,
  rsSparsePolynomial<T, TTol>* tmp1, rsSparsePolynomial<T, TTol>* tmp2, bool monic)
{
  rsAssert(a->_isCanonical());
  rsAssert(b->_isCanonical());

  a->tol = rsMax(a->tol, b->tol);
  while(!b->isZero())
  {
    *tmp1 = *b;
    rsSparsePolynomial<T, TTol>::divide(*a, *tmp1, tmp2, b);
    *a = *tmp1;
  }
  if(monic)
    a->makeMonic();

  // Notes:
  //
  // - Algorithm implementation has been adapted from rsRationalFunction<T>::polyGCD. I'm not sure,
  //   if we strictly require a,b to be canonical, but let's err to the conservative side. 
  //
  // - I think, the significance of the leading coeff of the result of the gcd algo may be: Assume 
  //   p and q have been produced via  p = g*a, q = g*b  where polynomials a,b have no common 
  //   divisors such that g is the gcd of p and q. If g happens to be non-monic, then calling 
  //   gcd(p, q, false) will restore g correctly including its leading coeff. ...I think. Verify! 
}


/**================================================================================================

ToDo:

- Implement unary plus. It's trivial but sometimes, we may want to use it at call sites for 
  clarity. But maybe it should return a (const?) reference rather than a value? Is that even 
  possible? In any case, we should make sure, that the unary + operator doesn't create a copy. A 
  unit test should verify that - perhaps by looking at the addresses of objects.

- Sprinkle in rsAssert(_isCanonical()); calls in all functions that assume a canonical 
  reprensentation in the spirit of defensive programming and contract based programming. Client
  code that uses the low level API and thereby messes up the canonical representation will fail
  early when we do this. ...ok done - not everywhere, though - but in a lot of places.

- We also do not yet use a relative tolerance anywhere. Maybe to facilitate this, we should provide 
  a member getScaledTolerance() or getAbsoluteTolerance that returns tol * getMaxAbsCoeff() where 
  getMaxAbsCoeff() should find the maximum absolute value of all of the coeffs. Maybe that function 
  should return a value of type TTol - not of type T. Then, whenever we need to actually use the 
  tolerance, we should retrieve it by calling getScaledTolerance(). Maybe we should use a new
  function name like rsMaxNorm<TNorm>(TArg x) that returns the maximum norm of the given x, e.g.
  for complex type, it would return max(abs(re), abs(im)) ...but maybe instead of invoking abs, it 
  should actually invoke a single argument variant rsMaxNorm<TNorm>(TArg x) that we may specialize
  for float, double, complex, etc explicitly

- Figure out what happens if client code uses negative powers. Currently, there's nothing that
  prevents this and maybe it could even make sense to allow it. But then the notion of degree
  gets murky. Maybe then there is indeed a difference between the degree and the max power in
  the case of an empty polynomial? Maybe, for the time being, we should trap attempts to set up
  terms with negative powers. This can later be relaxed, if needed. ...ok done: _isCanonical() now
  also verifies that the powers are all nonnegative.

- Implement root finding/factorization. Maybe we first need evaluation of derivatives. Maybe we can
  add an optimized function for evaluating the polynomial itself along with its 1st and 2nd 
  derivative. I think, the Laguerre root finding algorithm needs these. Look up the implementation
  in rsPolynomial where we implement it for the case of dense polynomials. Maybe try to adapt the
  implementation for the sparse case.

- Implement a function that evaluates p(x) and a given number M of its derivatives. The API should
  be like  evalWithDerivatives(T x, int M, T* f)  where f is an array of M (or maybe M+1) values
  for p(x), p'(x), p''(x), ... up to the M-th derivative. 

- What about copy- and move constructors and copy- and move assignment operators? Do we need to
  define them or can we rely on the auto-generated ones? It's important that swapping two
  sparse polynomials can be done allocation free. This is needed for inverting sparse filters by
  swapping numerator and denominator of their transfer functions (plus some extra stuff to 
  maintain the a0 = 1 normalization). This is an an operation that we need to do in a realtime 
  safe manner. Verify and document this! We are currently relying on the auto-generated copy- and
  move constructors and assignment operators and I think, this is totally fine.

- Implement composition of sparse polynomials (see free function rsComposeNaive() in Prototypes.h 
  file)

- Implement a lowestCommonMultiple() function, aka lcm - the cousin of gcd.

- Implement back and forth conversions between rsPolynomial and rsSparsePolynomial. In principle,
  both classes can be used for both purposes. It's just that the performance of the two 
  implementations is optimized for different cases. The conversion functions should be free 
  functions, I think. I don't want to couple the two classes too tightly to one another. Maybe
  We should have a class rsTypeConverter that has static functions for various type conversions.
  The function to convert matrices from one element type to another could then also go there.

- Make sure that for all member functions without underscore, at least one of the 3 things is 
  true:

    (1) We know that they don't mess up the canonical representation. In this case they should
        call rsAssert(_isCanonical()) at the end to document that. Well, maybe only in those cases
        where this is not trivially obvious.

    (2) They call _canonicalize() at the end. This is needed, if they potentially do destroy a
        canonical representation.

    (3) They call only other member functions without underscore, i.e. other members that are 
        already known to be safe.

  Then we can be sure that they always maintain a canonical representation.
   
- Implement += operator for right operand being another polynomial, a monomial, a constant. Do 
  the same for -=, *=, /=

- Maybe make the tolerance parameter for the constructors optional. I'm not sure about that, 
  though. It may invite forgetting to set it when it's really needed. But on the other hand, some
  types T don't need any tolerance at all. Maybe keep it mandatory for a while and make it optional
  later.


Notes:

- It might be tempting to write a constructor and/or setup function that takes a dense 
  polynomial, i.e. an object of type rsPolynomial<T>. But I think, that's not a good idea 
  because it would introduce unnecessary coupling.

*/