#pragma once 

// Functions that operate on std::vectors to perform polynomial coefficient array manipulations,
// translated from my python implementation

/* Evaluates the polynomial p a the given x using Horner's algorithm */
double polyEval(std::vector<double>& p, double x);

/* Truncates trailing zeros of the list p */
void polyTrunc(std::vector<double>& p, double tol = 0.0);

/* Makes the polynomial a monic, i.e. divides all coefficients by the  leading coefficient to make
the leading coefficient 1. Will result in division by zero error, if p is the zero polynomial. It 
works in place and will return the leading coefficient (which may or may not be of interest to the 
caller) */
double makeMonic(std::vector<double>& p);

/* Forms a weighted sum of the two coefficient lists p and q with weights wp and wq 
respectively. If the resulting list will have trailing zeros, these will be truncated. */
std::vector<double> polyAdd(
  const std::vector<double>& p, const std::vector<double>& q, 
  double tol = 0.0, double wp = 1, double wq = 1);

/* Subtracts the coefficient list q from the coefficient list p. If the result has trailing zeros, 
these will be truncated. */
std::vector<double> polySub(const std::vector<double>& p, const std::vector<double>& q,
  double tol = 0.0);

/* Multiplies two lists of polynomial coefficients by convolution. */
std::vector<double> polyMul(const std::vector<double>& p, const std::vector<double>& q,
  double tol = 0.0);

/* Divides polynomial p (product) by polynomial d (divisor) and returns 
the quotient in q and remainder in r */
void polyDivMod(std::vector<double> p, std::vector<double> d, 
  std::vector<double>& q, std::vector<double>& r, double tol = 0.0);

/* Quotient of polynomial division - this corresponds to the integer part of the division of 
natural numbers. */
std::vector<double> polyDiv(std::vector<double> p, std::vector<double> d, double tol);

/* Remainder of polynomial division */
std::vector<double> polyMod(std::vector<double> p, std::vector<double> d, double tol);

/* Checks, if vector v contains only zeros. */
bool isAllZeros(const std::vector<double>& v, double tol); // move to rapt

/* Computes the greatest common divisor of polynomials p and q which is defined as the polynomial 
of highest degree that divides both p and q. Such a polynomial is unique only up to multiplication
by a constant, so it is often additionally required to be a monic polynomial to make it unique. 
This normalization can be controlled by by the monic parameter. */
std::vector<double> polyGCD(
  const std::vector<double>& p, const std::vector<double>& q, double tol, bool monic = true);

/* Given the coefficient lists of two polynomials a(x) and b(x), this function computes the 
coefficient list of the polynomial c(x) that results from nesting a(x) and b(x) where a(x) is the
inner and b(x) the outer polynomial such that: c(x) = b(a(x)) */
std::vector<double> polyNest(const std::vector<double>& a, const std::vector<double>& b);

/* Reduces rational function p/q to the lowest possible denominator. */
void ratReduce(
  const std::vector<double>& pIn, const std::vector<double>& qIn,
  std::vector<double>& pOut, std::vector<double>& qOut, double tol);

/* Multiplies two rational functions represented as lists of coefficients for	numerator and 
denominator. Computes u/v = (p/q) * (r/s). By default, it will reduce the result to the lowest 
possible denominator but you can turn that off via the reduced parameter. */
void ratMul(
  const std::vector<double>& p, const std::vector<double>& q,
  const std::vector<double>& r, const std::vector<double>& s,
  std::vector<double>& u, std::vector<double>& v, double tol = 0.0, bool reduced = true);

/* Divides two rational functions */
void ratDiv(
  const std::vector<double>& p, const std::vector<double>& q,
  const std::vector<double>& r, const std::vector<double>& s,
  std::vector<double>& u, std::vector<double>& v, double tol = 0.0, bool reduced = true);

/* Adds two rational functions represented as lists of coefficients for numerator and denominator 
with optional weighting. It computes numerator and denominator of nr/dr = w1*n1/d1 + w2*n2/d2. */
void ratAdd(
  const std::vector<double>& n1, const std::vector<double>& d1,
  const std::vector<double>& n2, const std::vector<double>& d2,
  std::vector<double>& nr, std::vector<double>& dr, 
  double tol = 0.0, double w1 = 1, double w2 = 1);

/* Nesting of an inner rational function ni/di with an outer polynomial po. */
void ratPolyNest(
  const std::vector<double>& ni, const std::vector<double>& di,
  const std::vector<double>& po,
  std::vector<double>& nr, std::vector<double>& dr, double tol);

/** Nesting of two rational functions. The inner functions numerator and denominator are given by
nI, dI, likewise for the outer and nO, dO. The result is returned in nR, dR. */
void ratNest(
  const std::vector<double>& nI, const std::vector<double>& dI,
  const std::vector<double>& nO, const std::vector<double>& dO,
  std::vector<double>& nR, std::vector<double>& dR, double tol);