#ifndef RAPT_PIECEWISEPOLYNOMIAL_H
#define RAPT_PIECEWISEPOLYNOMIAL_H

/** A class for representing and performing computations with functions that are defined as 
piecewise polynomials. */

template<class T>
class rsPiecewisePolynomial
{

public:



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of pieces. */
  int getNumPieces() const { return (int) pieces.size(); }

  /** Returns a constant reference to the piece at index i. */
  const rsPolynomial<T>& getPieceConstRef(int i) const 
  { rsAssert(i >= 0 && i < getNumPieces(), "Index out of range"); return pieces[i]; }

  /** Returns the index of the piece, where x belongs or -1, if x is out of range to the left or
  getNumPieces() if x is out of range to the right. */
  int getIndex(T x) const;

  /** Returns the left boundary of the domain, where the function is defined. */
  T getDomainMinimum() const
  {
    if(domains.empty())
      return T(0);
    return domains[0];
  }

  /** Returns the right boundary of the domain, where the function is defined. */
  T getDomainMaximum() const
  {
    if(domains.empty())
      return T(0);
    return rsLast(domains);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Evaluation

  /** Evaluates the function at the given x. If x is outside the range where we have defined 
  pieces, it returns zero. */
  T evaluate(T x) const;

  /** Evaluation as oprator. */
  T operator()(T x) const { return evaluate(x); }


  //-----------------------------------------------------------------------------------------------
  // \name Manipulations

  /** Adds another piece to the object, possibly with a scalar weight applied to the new piece. */
  void addPiece(const rsPolynomial<T>& p, T pL, T pU, T weight = T(1));

  /** Clears the polynomial, setting it back into a freshly constructed state. */
  void clear() { domains.clear(); pieces.clear(); }

  /** Scales the whole function in the y-direction by the given factor. */
  void scale(T factor);

  /** Stretches the whole function in the x-direction by the given factor. */
  void stretch(T factor);

  /** Integrates the function. The integration constant determines the function value at the left 
  boundary. */
  void integrate(T c = T(0));

  rsPiecewisePolynomial<T> integral(T c = T(0)) const
  {
    rsPiecewisePolynomial<T> p = *this;
    p.integrate(c);
    return p;
  }


  // todo: derivative

  /** Shifts the pieces up or down in the y-direction such that they match at the segment 
  boundaries. */
  void makeContinuous();


  //-----------------------------------------------------------------------------------------------
  // \name Combination

  rsPiecewisePolynomial<T>& operator+=(const rsPiecewisePolynomial<T>& q) 
  { 
    for(int i = 0; i < q.getNumPieces(); i++)
      addPiece(q.pieces[i], q.domains[i], q.domains[i+1]);
    return *this;
  }

  rsPiecewisePolynomial<T> operator+(const rsPiecewisePolynomial<T>& q) const 
  { rsPiecewisePolynomial<T> r = *this; r += q; return r; }

  rsPiecewisePolynomial<T>& operator-=(const rsPiecewisePolynomial<T>& q) 
  { 
    for(int i = 0; i < q.getNumPieces(); i++)
      addPiece(q.pieces[i], q.domains[i], q.domains[i+1], T(-1));
    return *this;
  }

  rsPiecewisePolynomial<T> operator-(const rsPiecewisePolynomial<T>& q) const 
  { rsPiecewisePolynomial<T> r = *this; r -= q; return r; }

  // todo: implement multiplication, maybe the addPiece function can be extended to also handle 
  // multiplication, such that we do not have to replicate the splitting logic


  /** Convolves two polynomial pieces p(x) and q(x) that are defined on the domains pL..pU and 
  qL..qU respectively and assumed to be zero outside these domains (L and U stand for lower and 
  upper boundaries of the domains). The result are 3 polynomial pieces rL,rM,rR that are adjacent 
  to each other (L,M,R stand for left, middle, right). These polynomials are output parameters and
  assigend by the function. These output pieces are defined on the domains rLL..rLU, rLU..rRL, 
  rRL..rRR respectively which are also output parameters. The left/right boundaries of the middle 
  segment are the same as the right/left boundaries of the adjacent left and right pieces, so there 
  are no output parameters for the middle section's domain boundaries because they would be 
  redundant. We get 3 pieces as output because for the left piece, the two input segments that are 
  convolved, are not yet fully overlapping as q gets shifted over p and for the right piece, they 
  are not fully overlapping anymore. Only in the middle segment, there's full overlap, i.e. q is 
  fully inside p or the other way around. */
  static void convolvePieces(
    const rsPolynomial<T>& p, T pL, T pU, const rsPolynomial<T>& q, T qL, T qU,
    rsPolynomial<T>& rL, T& rLL, T& rLU, rsPolynomial<T>& rM, rsPolynomial<T>& rR, T& rRL, T& rRU);

  /** Convolves this piecewise polynomial with another piecewise polynomial p and returns the 
  result which is again a piecewise polynomial. */
  rsPiecewisePolynomial<T> convolve(const rsPiecewisePolynomial<T>& p);


  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Creates an Irwin-Hall distribution of given order N between a and b. This is the uniform 
  distribution convolved with itself N times. It arises as amplitude distribution, when you add up 
  the outputs of N independent noise generators with uniform distributions between a..b. N=0 gives 
  the uniform distribution, N=1 gives a triangular (piecewise linear) distribution, N=2 a piecewise
  parabolic distribution and so on. The higher order distributions approach a Gaussian distribution
  due to the central limit theorem. */
  static rsPiecewisePolynomial<T> irwinHall(int N, T a = T(0), T b = T(1));


protected:

  std::vector<T> domains;  // maybe rename to boundaries, ends, limits, borders
  std::vector<rsPolynomial<T>> pieces;
  // -the pieces are adjacent (no verlap, no gaps)
  // -piece[i] goes from domains[i] to domains[i+1]

  // todo: maybe have an evaluation mode parameter that determines what happens when the user wants 
  // evaluate outside the domain. possible values are: zero (as it is now), just use the left/right
  // polynomials for extrapolation, clamp output at whatever left/right polynomials give at the
  // boundaries - so modes could be named: zero, extrapolate, clamp ..maybe, the extrapolation 
  // could be restricted to use only the lower order coeffs
  // 

};

#endif