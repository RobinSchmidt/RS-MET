#pragma once 

// file is obsolete - code has been moved into class rsRationalFunction

// Functions that operate on std::vectors to perform polynomial coefficient array manipulations,
// translated from my python implementation. They should be moved into rsRationalFunction as static
// meber functions. They are sort of low-level, although they use std::vector...maybe mid-level, but
// they may need to resize the vectors. maybe factor out true low-level functions (operating on raw 
// arrays) and if the need to resize, they don't actually resize anything but just inform the caller
// about the new size by a return value.

/* Evaluates the polynomial p a the given x using Horner's algorithm */
template<class T>
T polyEval(std::vector<T>& p, T x);

/* Truncates trailing zeros of the list p */
template<class T>
void polyTrunc(std::vector<T>& p, T tol = 0.0);

/* Makes the polynomial a monic, i.e. divides all coefficients by the  leading coefficient to make
the leading coefficient 1. Will result in division by zero error, if p is the zero polynomial. It 
works in place and will return the leading coefficient (which may or may not be of interest to the 
caller) */
template<class T>
T makeMonic(std::vector<T>& p);

/* Forms a weighted sum of the two coefficient lists p and q with weights wp and wq 
respectively. If the resulting list will have trailing zeros, these will be truncated. */
template<class T>
std::vector<T> polyAdd(
  const std::vector<T>& p, const std::vector<T>& q, 
  T tol = 0.0, T wp = 1, T wq = 1);

/* Subtracts the coefficient list q from the coefficient list p. If the result has trailing zeros, 
these will be truncated. */
template<class T>
std::vector<T> polySub(const std::vector<T>& p, const std::vector<T>& q,
  T tol = 0.0);

/* Multiplies two lists of polynomial coefficients by convolution. */
template<class T>
std::vector<T> polyMul(const std::vector<T>& p, const std::vector<T>& q,
  T tol = 0.0);

/* Divides polynomial p (product) by polynomial d (divisor) and returns 
the quotient in q and remainder in r */
template<class T>
void polyDivMod(std::vector<T> p, std::vector<T> d, 
  std::vector<T>& q, std::vector<T>& r, T tol = 0.0);

/* Quotient of polynomial division - this corresponds to the integer part of the division of 
natural numbers. */
template<class T>
std::vector<T> polyDiv(std::vector<T> p, std::vector<T> d, T tol);

/* Remainder of polynomial division */
template<class T>
std::vector<T> polyMod(std::vector<T> p, std::vector<T> d, T tol);

/* Checks, if vector v contains only zeros. */
template<class T>
bool isAllZeros(const std::vector<T>& v, T tol); // move to rapt

/* Computes the greatest common divisor of polynomials p and q which is defined as the polynomial 
of highest degree that divides both p and q. Such a polynomial is unique only up to multiplication
by a constant, so it is often additionally required to be a monic polynomial to make it unique. 
This normalization can be controlled by by the monic parameter. */
template<class T>
std::vector<T> polyGCD(
  const std::vector<T>& p, const std::vector<T>& q, T tol, bool monic = true);

/* Given the coefficient lists of two polynomials a(x) and b(x), this function computes the 
coefficient list of the polynomial c(x) that results from nesting a(x) and b(x) where a(x) is the
inner and b(x) the outer polynomial such that: c(x) = b(a(x)) */
template<class T>
std::vector<T> polyNest(const std::vector<T>& a, const std::vector<T>& b);

/* Reduces rational function p/q to the lowest possible denominator. */
template<class T>
void ratReduce(
  const std::vector<T>& pIn, const std::vector<T>& qIn,
  std::vector<T>& pOut, std::vector<T>& qOut, T tol);

/* Multiplies two rational functions represented as lists of coefficients for	numerator and 
denominator. Computes u/v = (p/q) * (r/s). By default, it will reduce the result to the lowest 
possible denominator but you can turn that off via the reduced parameter. */
template<class T>
void ratMul(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol = 0.0, bool reduced = true);

/* Divides two rational functions */
template<class T>
void ratDiv(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol = 0.0, bool reduced = true);

/* Adds two rational functions represented as lists of coefficients for numerator and denominator 
with optional weighting. It computes numerator and denominator of nr/dr = w1*n1/d1 + w2*n2/d2. */
template<class T>
void ratAdd(
  const std::vector<T>& n1, const std::vector<T>& d1,
  const std::vector<T>& n2, const std::vector<T>& d2,
  std::vector<T>& nr, std::vector<T>& dr, 
  T tol = T(0), T w1 = T(1), T w2 = T(1));

/* Nesting of an inner rational function ni/di with an outer polynomial po. */
template<class T>
void ratPolyNest(
  const std::vector<T>& ni, const std::vector<T>& di,
  const std::vector<T>& po,
  std::vector<T>& nr, std::vector<T>& dr, T tol = 0.0);

/** Nesting of two rational functions. The inner functions numerator and denominator are given by
nI, dI, likewise for the outer and nO, dO. The result is returned in nR, dR. */
template<class T>
void ratNest(
  const std::vector<T>& nI, const std::vector<T>& dI,
  const std::vector<T>& nO, const std::vector<T>& dO,
  std::vector<T>& nR, std::vector<T>& dR, T tol = 0.0);