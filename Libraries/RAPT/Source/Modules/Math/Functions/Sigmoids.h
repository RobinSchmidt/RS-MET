#ifndef RAPT_SIGMOIDS_H_INCLUDED
#define RAPT_SIGMOIDS_H_INCLUDED

/** This class is a collection of saturation functions that are supposed to be applied to positive
input values only. To get a value for a negative input value x, use -f(-x), i.e. symmetrize it 
manually. */

template<class T>
class PositiveSigmoids
{

public:

  /** ... */
  static T linear(T x);

  /** Implements y = x / (1+x). */
  static T rational(T x);

  /** Implements y = (x + x^2 + x^3) / (x + x^2 + x^3 + 1). This is the special case of the 
  parametric sigmoid function implemented in class ParametricSigmoid with the parameter value
  y1 = 0.75, which is the critical value at which the 2nd derivative of the core function becomes 
  0 (see comments of the class for more details). */
  static T cubicRational(T x);

  /** This function implements a cubic polynomial y = x + a3*x^3. It satisfies the conditions:
  f(0)=0, f'(0)=1, f''(0)=0 by design, regardless of the value of a3. The coefficient a3 is chosen
  according to the condition: f(s)=1 where s is the maximum at which f'(s)=0 holds. The location of
  this maximum comes out as s=1.5 and a3=-4/27. Input values larger than 1.5 are mapped to unity, 
  i.e. saturated at unity. */
  static T cubic(T x);

  /** This function implements a quartic polynomial y = x + a3*x^3 + a4*x^4. In addition to the 
  conditions that are satisfied by the cubic, this function also satisfies f''(s)=0 which
  makes the function smoother at the junction between nonlinear and constant zone. The maximum of
  the polynomial (i.e. the saturation level) occurs at x=2 and we have a3=-1/4, a4=1/16. */
  static T quartic(T x);

  /** This function implements a hexic polynomial y = x + a4*x^4 + a5*x^5 + a6*x^6. In addition to 
  the conditions that are satisfied by the quartic, it also satisfies f'''(0)=f'''(s)=0 where the 
  saturation level is given by s=2. */
  static T hexic(T x);

  /** This function is an identity function for x < 0.5, equal to 1 for x > 1 and a third order 
  order polynomial in between, such that at x=0.5, the function value and 1st derivative is 
  matched. */
  //static double softClip(double x);
    // this function is not good - it goes above the identity function - remove


  /** Soft clipping using a hexic polynomial with threshold t. Below t, the output is the identity
  function, above t, a properly scaled and shifted sextic polynomial will be used. */
  static T softClipHexic(T x, T t);

  /** Special case of softClipHexic(double x, double t) where t = 0.5. */
  static T softClipHexic(T x);

};

//-------------------------------------------------------------------------------------------------

/** This class is a collection of normalized sigmoid-shaped functions. They are normalized in the 
sense that they go through the output range from -1 to +1 as the input x goes from -inf to +inf and 
they also have unit slope at the origin. */

template<class T>
class NormalizedSigmoids
{

public:

  /** Clips the input to the range -1..+1. */
  static T clip(T x);

  /** Normalized arc-tangent function. */
  static T atan(T x);

  /** Normalized arc-tangent function. */
  static T tanh(T x);

  /** A family of sigmoid curves with a parameter p, realizing the function: 
  y = x / (1 + |x|^p)^(1/p) = x * (1 + |x|^p)^(-1/p)
  The parameter p determines, how fast the saturation level will be reached. It should be p >= 1.
  When p approaches infinity, the function approaches the hard clipper. ...beware of numerical 
  problems with large p though. */
  static T powRatio(T x, T p);

  // symmetrized versions of the corresponding positive-range functions:
  static T rational(T x);
  static T cubicRational(T x);
  static T cubic(T x);
  static T quartic(T x);
  static T hexic(T x);
  static T softClipHexic(T x);

};

//=================================================================================================

/**

This class implements a parametrized sigmoid function. The parameter y1 adjusts the value at x=1, 
so y1 = f(x=1). The value y1 has to in the range 0.5..1.0 (ends inclusive). If you set y1 = 1, you 
will get a hardclipper. If you set y1 = 0.5, you will get the function y = x / (1+|x|). There are 
some deeper parameters, but if you leave them at their default settings and just adjust y1, the 
function will be 2nd order continuous for y1 >= 0.75 and 1st order continuous for y1 < 0.75.

The function is based on the core-function:

        x + a*(b*x^2 + (1-b)*x^3)
f(x) = -------------------------------
        x + a*(b*x^2 + (1-b)*x^3) + 1

where a and b are adjustable parameters. The a-parameter is directly related to the function value 
at unity y1 defined as y = f(x=1). Because it's more intuitive, the a-parameter will typically set 
in terms of that value using setValueAt1() which can take value in the range 0.5..1 
(ends inclusive). As y1 gets closer to 1, the function f(x) becomes kinda strange looking (getting 
turning points etc.), that's why when y1 is set above some critical breakpoint value yb, we 
actually use a piecewise function that is the identity-function below a threshold and an 
appropriately scaled-and-shifted version of the the original function f above the threshhold. That 
breakpoint value is also adjustable, the default is yb = 0.75.  The secondary parameter b adjusts 
how the function is distributed between the range x < 1 and x > 1. High values of b give more 
weight to the left side, low values give more weight to the right side. We should have b <= 1 for 
the function to be bounded for all x. Currently, the b-parameter is not exposed as user-parameter 
but instead computed automatically from the a-parameter via the formula b = 1 / rsMax(1.0, a). This
formula ensures that the 2nd derivative of f(x) vanishes at x=0 (for a >= 1, corresponding to 
y1 >= 0.75).

ToDo:
-maybe have a function to compute the derivative and/or value-and-derivative at the same time. This
 might be useful for Newton iteration, when the function is used in a ZDF-filter feedback path. Use
 quotient rule for computation (might be more efficient than computing the derivative from scratch)

*/

template<class T>
class ParametricSigmoid
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  ParametricSigmoid();


  /** \name Setup */

  /** Sets the desired output-value when the input value is unity. The value should be in the
  range 0.5..1 (inclusive). */
  void setValueAt1(T newValue);

  /** Sets the threshold above which the nonlinear range begins. Can be used alternatively to
  setValueAt1. */
  void setThreshold(T newThreshold);

  /** Beyond some breakpoint for the value passed into setValueAt1, we will switch to a piecewise
  defined function that is an identity function in the lower range and an appropriately scaled
  and shifted version of the original function in the upper range. This breakpoint is set with this
  function. The most reasonable values are in the range 0.7..0.8. */
  void setPiecewiseBreakpoint(T newBreakpoint);


  /** \name Function Evaluation */

  /** Implements the core function 
  y = f(x) = (x + a*(b*x^2 + (1-b)*x^3)) / (x + a*(b*x^2 + (1-b)*x^3) + 1) for some given parameters
  a and b. This is mainly for internal use. */
  static T coreFunction(T x, T a, T b);

  /** Returns an output value for given x. */
  inline T getValue(T x);

  // todo: getValueAndDerivative


protected:

  /** \name Misc */

  /** Computes the a-parameter form the desired value y1 = f(x=1) */
  T getA(T y1);

  /** Computes the b-parameter from a given a-parameter. */
  T getB(T a);

  /** Computes the coefficient "a" and the scale/shift coefficients sx,sy,ty from y1,yb. */
  void computeCoeffs();

  /** Implements the core-function (for postive x, without possible scale-and shift) with parameters
  given by our member variables. */
  inline T coreFunction(T x);


  /** \name Data */

  T y1;         // value y at x=1
  T yb;         // breakpoint for y1 beyond which we switch to a piecewise function
  T a, b;       // parameters
  T c2, c3;     // coeffs for x^2 and x^3
  T sx, sy, ty; // scale/shift coefficients for the piecewise case

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
inline T ParametricSigmoid<T>::coreFunction(T x)
{
  x *= 1 + c2*x + c3*x*x;    // x + c2*x^2 + c3*x^3 = x + a*(b*x^2 + (1-b)*x^3)
  return x / (x+1);          // (x + a*(b*x^2 + (1-b)*x^3)) / (x + a*(b*x^2 + (1-b)*x^3) + 1)
}

template<class T>
inline T ParametricSigmoid<T>::getValue(T x)
{
  T t = rsAbs(x);
  if(t < ty)
    return x;                                                // identity below threshold
  else
    return rsSign(x) * (ty + sy * coreFunction(sx*(t-ty)));  // scaled/shifted core function
}

#endif
