
template<class T>
void rsPiecewisePolynomial<T>::addPiece(const rsPolynomial<T>& p, T pL, T pU, T weight)
{
  T tol(0);  // tolerance for (floating point) equality comparisons -> make member
  auto match = [&](T x, T y) -> bool { return rsAbs(x-y) <= tol; };

  rsAssert(pL < pU);
  if(pieces.empty()) {             // initialize with the first piece
    pieces.push_back(p);
    domains.push_back(pL);
    domains.push_back(pU);
    return;  }

  // I think, these cases can now be subsumed by the code below...but maybe it's neverless useful 
  // to handle them specially for efficiency - these cases here (especially the first) come up in 
  // the supposedly common case of building up a piecewise polynomial from scratch by appending one
  // segment after another:
  if(match(pL, rsLast(domains))) { // append piece at the right end
    pieces.push_back(p);
    domains.push_back(pU);
    return;  }
  if(match(pU, domains[0])) {      // prepend piece at the left end
    rsPrepend(pieces,  p);
    rsPrepend(domains, pL);
    return; }

  // Figure out start- and end indices for segment:
  int numPieces = getNumPieces();
  int iL = getIndex(pL);
  int iU;
  if(pU >= rsLast(domains))
    iU = numPieces;
  else
    iU = getIndex(pU);

  // (don't) handle gaps:
  if(iL == numPieces || iU == -1) {  // p starts after this or ends before this
    rsError("Gaps are not allowed"); // we can't handle them with the current implementation
    return; }

  // Function to split the piece at index i into two pieces at x0:
  auto split = [&](int i, T x0) 
  { 
    rsAssert(x0 > domains[i] && x0 < domains[i+1], "x0 outside domain of piece i");
    rsInsert(domains, x0, i+1);
    rsInsert(pieces, pieces[i], i);
  };
  // maybe move out and make it a member function, maybe also have a merge function.

  // Function to accumulate polynomial q into polynomial p:
  auto accumulate = [&](rsPolynomial<T>& p, const rsPolynomial<T>& q)
  { p.addWithWeight(q, weight); };
  // todo: use a std::function that is passed as parameter, so we can use the same function also 
  // for subtraction, multiplication, etc. - but for multiplication or division, we have to do 
  // something else in the case where we prepend or append the pieces here...we'll see...
  // maybe just implement addition, subtraction and scalar multiplication (as operators) - then
  // we have at least a vector-space. elementwise multiplication and division can come later. 
  // mulitplication may actually shrink the domain because the result is nonzero only where both
  // factors are nonzero, so we need some additional logic to update the domains - but maybe we 
  // don't need to mainpulate the domains - we can just let some sections be the zero polynomial
  // ...or we could cut off zero sections at start and end as post-processing
  // ...maybe this function should have an additional "mode" parameter, 0: add, 1: multiply, 
  // 2: divide...but no - divide doe not really make sense - from a mathematical perspective, we
  // would expect a piecewise rational function to come out and not whatever results from 
  // polynomial division (the results are equal only if this is divisible by p)


  // Handle case when left boundary of p is to the left of our left boundary:
  if(iL == -1) {
    rsPrepend(domains, pL);
    rsPrepend(pieces,  p);   // this works only for addition, not multiplication
    iL = 1;
    pL = domains[iL];        // so subsequent code can be used as if we had a match
    iU++;                    // because we have a piece more now
    numPieces++; }

  // Handle case when left boundary of p is in the middle of one of our existing pieces:
  if(!match(pL, domains[iL])) {
    split(iL, pL); iL++; iU++; numPieces++; }

  // Accumulate the new polynomial into the exisitng pieces that are now in full overlap:
  for(int i = iL; i < iU; i++)
    accumulate(pieces[i], p);

  // If end of new piece is aligned with end of an existing piece, we are done:
  if(match(pU, domains[iU]))
    return;

  // ...and if it's not aligned, we have two cases to consider:
  if(iU == numPieces) {          // End of new piece is beyond or right boundary 
    domains.push_back(pU);
    pieces.push_back(p);   }     // this works only for addition, not multiplication
  else  {                        // End of new piece is in the middle of some existing piece
    split(iU, pU);
    accumulate(pieces[iU], p); }
}

template<class T>
void rsPiecewisePolynomial<T>::scale(T factor)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].scale(factor);
}

template<class T>
void rsPiecewisePolynomial<T>::stretch(T factor)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].stretch(factor);
  for(size_t i = 0; i < domains.size(); i++)
    domains[i] *= factor;
}

template<class T>
void rsPiecewisePolynomial<T>::integrate(T c)
{
  for(size_t i = 0; i < pieces.size(); i++)
    pieces[i].integrate();
  T x  = domains[0];
  T yL = pieces[0](x);
  pieces[0].shiftY(c-yL);  // adjust start value
  makeContinuous();
}

template<class T>
void rsPiecewisePolynomial<T>::makeContinuous()
{
  for(size_t i = 1; i < pieces.size(); i++)
  {
    T x  = domains[i];
    T yL = pieces[i-1](x);
    T yR = pieces[i](x);
    pieces[i].shiftY(yL-yR);
  }
}

template<class T>
int rsPiecewisePolynomial<T>::getIndex(T x) const
{
  if(pieces.empty() || x < domains[0])
    return -1;
  if( x >= rsLast(domains))
    return getNumPieces();
  int i = rsArrayTools::findSplitIndex(&domains[0], (int) domains.size(), x);
  if(x < domains[i])
    i--;
  rsAssert(i >= -1 && i < getNumPieces());
  return i;
}

template<class T>
T rsPiecewisePolynomial<T>::evaluate(T x) const
{
  int i = getIndex(x);
  if(i < 0 || i >= getNumPieces()) 
    return T(0);
  else        
    return pieces[i](x);  // use evaluate function (needs to be written)
}

template<class T>
void rsPiecewisePolynomial<T>::convolvePieces(
  const rsPolynomial<T>& p, T pL, T pU, const rsPolynomial<T>& q, T qL, T qU,
  rsPolynomial<T>& rL, T& rLL, T& rLU, rsPolynomial<T>& rM, rsPolynomial<T>& rR, T& rRL, T& rRU)
{
  // Check, which one of p,q has higher degree and invoke the function recursively with swapped 
  // inputs if they are in disadvantageous order. We don't want q to be the higher degree 
  // polynomial because it's q that gets blown up to a bivariate polynomial. Mathematically, the 
  // results should be the same, but using the lower degree polynomial as q should be more 
  // efficient and maybe also more precise. The unit tests should also pass when we comment this 
  // optimization out:
  if(q.getDegree() > p.getDegree()) {
    convolvePieces(q, qL, qU, p, pL, pU, rL, rLL, rLU, rM, rR, rRL, rRU);
    return; }

  // Compute the domains for the 3 output segments:
  T wp = pU - pL;   // width of domain of p
  T wq = qU - qL;   // width of domain of q
  rLL  = pL + qL;
  rLU  = pL + qU;   // == rLL + wq
  rRL  = pU + qL;
  rRU  = pU + qU;   // == rRL + wq
  if(wq > wp)
    rsSwap(rLU, rRL);

  // Create the bivariate polynomial PQ(x,y) = p(y)*q(x-y) which is our integrand in the 
  // convolution integral:
  using BiPoly = rsBivariatePolynomial<T>;
  BiPoly Q  = BiPoly::composeWithLinear(q, T(1), T(-1)); //  Q(x,y) = q(x-y)
  BiPoly PQ = Q.multiplyY(p);                            // PQ(x,y) = p(y)*q(x-y)

  // Integrate out the dummy variable y (typically tau in literature about convolution, and our x 
  // is their t), leaving a polynomial only in x (aka t). We get 3 segments:
  using Poly = rsPolynomial<T>;
  Poly a({-qU, 1});                           // lower integration limit (in some cases)
  Poly b({-qL, 1});                           // upper integration limit (in some cases)
  rL                  = PQ.integralY(pL, b);  // left segment
  if(     wp > wq) rM = PQ.integralY(a,  b);  // middle segment when wp longer than wq
  else if(wq > wp) rM = PQ.integralY(pL, pU); // middle segment when wq longer than wp
  else             rM = Poly(0);              // p,q have same length -> no middle segment
  rR                  = PQ.integralY(a,  pU); // right segment

  // Adjust degrees by truncating the trailing zero coefficients. I think, they arise because the
  // matrix of coefficients of the bivariate polynomial Q is triangular. The integral could 
  // potentially produce higher order nonzero coeffs but doesn't due to the special structure of Q.
  int deg = 0;
  if(wp > wq) deg = p.getDegree();
  if(wp < wq) deg = q.getDegree();
  rM.setAllocatedDegree(deg);
  deg = p.getDegree() + q.getDegree() + 1;
  rL.setAllocatedDegree(deg);
  rR.setAllocatedDegree(deg);

  // old:
  // Cut off trailing zero coefficients in the produced segments:
  //rL.truncateTrailingZeros();
  //rM.truncateTrailingZeros();
  //rR.truncateTrailingZeros();
  // The middle section may still have close-to-zero trailing coeffs (i think, it happens only 
  // when wp > wq but this needs more tests)
  // i think, the degrees are degP+degQ+1 for the L/R sections and degP or deQ for the middle M
  // section and which one of the two it is, is determined by which one has the longer domain

  // Notes:
  // -The expressions for the integration limits were found by trial and error and need more tests,
  //  especially, when q has longer support than p. Perhaps, we should switch the roles of p and q
  //  in such a case, but maybe it's more advisable to select the roles of p and q by their 
  //  degrees. We should probably create the bivariate polynomial from whichever has lower degree.
  //  ...we'll see
  // ToDo:
  // -make a function that uses a workspace - if we later need to convolve many pieces in a 
  //  piecewise polynomial, we'll call this in a double-loop: each piece from one input is 
  //  convolved with each piece from the other (and then all the results get added up), so we want
  //  the operation to be efficient
}

template<class T>
rsPiecewisePolynomial<T> rsPiecewisePolynomial<T>::convolve(const rsPiecewisePolynomial<T>& q)
{
  rsPiecewisePolynomial<T> r;
  using Poly = rsPolynomial<T>;
  Poly rL, rM, rR;
  T rLL, rLU, rRL, rRU;
  for(int i = 0; i < getNumPieces(); i++) {
    for(int j = 0; j < q.getNumPieces(); j++) {
      const Poly& pi = getPieceConstRef(i);
      T pL = domains[i];
      T pU = domains[i+1];
      const Poly& qj = q.getPieceConstRef(j);
      T qL = q.domains[j];
      T qU = q.domains[j+1];
      convolvePieces(pi, pL, pU, qj, qL, qU, rL, rLL, rLU, rM, rR, rRL, rRU);
      r.addPiece(rL, rLL, rLU);
      if(rLU < rRL)
        r.addPiece(rM, rLU, rRL);
      r.addPiece(rR, rRL, rRU);   }}
  return r;
}

template<class T>
rsPiecewisePolynomial<T> rsPiecewisePolynomial<T>::irwinHall(int order, T a, T b)
{
  rsAssert(order >= 0);
  rsPiecewisePolynomial<T> p0; 
  p0.addPiece(rsPolynomial<T>({ T(1)/(b-a) }), a, b);  // our seed function
  rsPiecewisePolynomial<T> p = p0;
  for(int i = 1; i <= order; i++)
    p = p.convolve(p0);
  return p;
}
