#ifndef RAPT_ROOTFINDER_H_INCLUDED
#define RAPT_ROOTFINDER_H_INCLUDED

/** This class implements various (one-dimensional) root finding algorithms, i.e. it finds 
solutions to the equation f(x) = 0 for an arbitrary given function f(x). For more generality, the
right hand side does not actually need to be equal to zero but can be any target value, so it 
actually finds solutions to f(x) = y for given f(x) and given y (which defaults to zero). This 
little addition to standard textbook root-finding algorithms doesn't change much algorithmically 
(we just need to subtract the target value in each function evaluation), yet adds a lot to the 
flexibility and ease of use. Some algorithms require the user to pass an initial interval that is 
assumed to bracket the root, others require an initial estimate of the root. Some higher level 
functions don't require anything like that - but these functions will need to make a guess for the 
bracket internally which may lead to suboptimal performance - because, you know, good guesses are 
hard to come by when you don't have any information to work with. So, these high-level functions
are mostly meant for quick-and-dirty prototype implementations. If you know a bit more about your 
particular function at hand, it's a good idea to use that knowledge together with the lower level 
root-finding functions.

References:

  (1) Numerical Recipies in C (2nd Edition), Chapter 9

*/

template<class T>
class rsRootFinder
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Bracketing Methods

  /** A high-level function that needs only the function, the desired y-value (defaulting to zero) 
  and an optional initial guess for x where the root/y-value might be found (also defaulting to 
  zero). It's recommended to be used only when you don't have any guess for the initial root 
  bracket available. It uses a rather dumb bracket-guessing algorithm internally, namely the one
  implemented in the findBracket() function. Also, it currently uses the simple bisection method 
  for the actual root finding whose convergence rate is not really great. The function is meant for 
  quick and dirty prototype stuff to get something running. Later, you may want to refine the code 
  by coming up with a better initial guess/bracket and picking the most suitable algorithm for the 
  problem at hand yourself. For this, you'll have to use the lower level functions of this class 
  yourself.  */
  static T findRoot(const std::function<T(T)>& f, T y = T(0), T x0 = T(0));

  /** Bisection takes an initial interval xLeft, xRight (assumed to bracket the root) and evaluates 
  the function at the midpoint of the interval. Depending on the function value, the midpoint 
  becomes either the new left or the new right border of the bracketing interval. So, in each 
  iteration, the size of the interval is halved. The order of convergence is linear (in each 
  iteration, we get one more correct binary digit in the root estimate) and convergence is 
  guaranteed. That means: the method is slow but safe. */
  static T bisection(const std::function<T(T)>& func, T xLeft, T xRight, T y = 0);
  // -Maybe instead of std::function use a second template parameter F
  // -Maybe declare the template parameters in front of the functions, not the class (like in 
  //  rsArrayTools)
  // -Let the function take a tolerance parameter (maybe defaulting 
  //  std::numeric_limits<T>::epsilon)
  // -Let the function take a maxNumIterations paremeter.
  // -When adding these additional parameters, make sure that their API and semantics matches that
  //  in rsMinimizer1D (in Optimization.h/cpp). 

  /** Similar to bisection but doesn't use the midpoint of the current bracketing interval, but the
  point where a line between (xLeft,yLeft), (xRight,yRight) crosses the x-axis. Convergence is
  typically (for smooth functions, that are well approximated by a line near the root) faster than
  for bisection and convergence is guranteed. But for some functions, the convergence may also be 
  slower than bisection. The method is also known as "regula falsi". */
  static T falsePosition(const std::function<T(T)>& func, T xLeft, T xRight, T y = 0);

  // ToDo: 
  // -Let caller pass a tolerance (maybe separate for x and y), state pathological conditions
  //  when no convergence can be epxected (like when the function has poles like 1/x)
  // -secant, ridders, brent, ... see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp


  //static T modifiedFalsePosition(std::function<T(T)>& func, T xLeft, T xRight, T y = 0);

  /** Many root finding algorithms require the user to provide an initial guess for the interval in
  which the root is to be found. This is called "root bracketing". The root is assumed to lie 
  somewhere in the bracket between xL and xR. If you don't have any such guess for the bracket, you 
  can use this function to produce one. It's using a simple algorithm that starts at an initial 
  point x0 (defaulting to 0) and then progressively expands an interval [xL, xR] around that point 
  until f(xL) <= y <= f(xR). It's not very sophisticated and guaranteed to work only for functions 
  that satisfy (f(-inf) = -inf and f(+inf) = +inf) or (f(-inf) = +inf and f(+inf) = -inf). I 
  recommend to try to avoid using it if you have any better way of coming up with an initial guess 
  for the interval. But sometimes, this dirty guesswork is just needed. Numerical analysis can be a 
  messy business. */
  static void findBracket(const std::function<T(T)>& f, T* xL, T* xR, T y = T(0), T x0 = T(0));
  // -Maybe return num-iterations (or error-code).
  // -Maybe have a function isBracket that can be used in unit tests as well as in assertions
  //  inside higher level functions (like findRoot) that use findBracket



  // I think, this code is obsolete now - we should now use findBracket() instead
  /** For a given function f = f(x) and target vaule y, finds a reasonable left bracket xL such 
  that f(xL) <= y. You can pass an initial guess for xL as well as an initial distance d by which
  xL will be decremented, in case, the condition f(xL) doesn't hold. */
  //static inline T findLeftBracket(const std::function<T(T)>& f, T y, T xL = T(0), T d = T(1))
  //{
  //  while(f(xL) > y) { xL -= d; d *= 2; }
  //  return xL;
  //}
  //// hmm - but this function works only for increasing functions
  //// maybe the growth multiplier for d should be an optional user parameter defaulting to 2. 
  //// When doing so, we should assert that it is >= 1.

  //static inline T findRightBracket(const std::function<T(T)>& f, T y, T xR = T(0), T d = T(1))
  //{
  //  while(f(xR) < y) { xR += d; d *= 2; }
  //  return xR;
  //}



  //-----------------------------------------------------------------------------------------------
  // \name Derivative Methods

  /** Computes the delta for one update step in the Newton iteration method. The update step uses 
  the formula  xNew = xOld + dx  with: 
  
          -f(x)
    dx = -------
          f'(x)
  
  This function computes the dx and can be used within the Newton iteration like  
  x += newtonStep(f, fp)  after you have computed function value f and derivative fp at the 
  current estimate for x. In practice, you'll probably want to assign the dx to a variable, though
  - so you can check the convergence  criterion. 
  See:  https://en.wikipedia.org/wiki/Newton%27s_method  */
  static T stepNewton(const T& f, const T& fp) { return -f/fp; }

  /** Implements root finding via Newton iteration. The function func should take the x value as 
  1st parameter and produce the value and derivative in 2nd and 3rd parameter respectively. These
  are output parameters and passed by pointer. The reason to use a single function to compute value
  and derivative is that it often happens that it's more efficient to evaluate a function and its
  derivative at the same time rather than starting completely from scratch for evaluating the 
  derivative. */
  static T newton(const std::function<void(T, T*, T*)>& func, T xGuess, T y = 0);


  /** Computes the delta for one update step in the Halley iteration method.  The 
  update step uses the formula xNew = xOld + dx with:
  
             2 * f * f1
    dx = ---------------------
          f * f2  -  2 * f1^2
  
  where f = f(x), f1 = f'(x), f2 = f''(x). See https://en.wikipedia.org/wiki/Halley%27s_method
  The method converges faster than Newton's - at least in theory. But it has the additional cost
  per step to compute the 2nd derivative. If it's practically advantageous to use Halley over 
  Newton may depend on the problem at hand. */
  static T stepHalley(const T& f, const T& f1, const T& f2) { return (2*f*f1) / (f*f2 - 2*f1*f1); }

  /** Implements root finding via Halley iteration. See newton() - the API is analoguous with the
  obvious difference that func must now also produce the 2nd derivative. */
  static T halley(const std::function<void(T, T*, T*, T*)>& func, T xGuess, T y = 0);


  /** Computes the delta for one update step in the 3rd order Householder iteration method. The 
  update step uses the formula xNew = xOld + dx with:
  
                3 * f^2 * f2  -  6 * f * f1^2
    dx = -------------------------------------------
          6 * f1^3  -  6 * f * f1 * f2  +  f^2 * f3

  where f = f(x), f1 = f'(x), f2 = f''(x), f3 = f'''(x). See:
  https://en.wikipedia.org/wiki/Householder%27s_method#Example 
  The Newton and Halley methods are the 1st and 2nd order Householder methods respectively. They 
  have their own names but belong conceptually to the Householder family. In theory, the higher the
  order, the faster the convergence. In practice - well...dunno - I guess, there are confounding 
  factors like rounding errors. ...TBC...figure out! */
  static T stepHouseholder3(const T& f, const T& f1, const T& f2, const T& f3) 
  { return (3*f*f*f2 - 6*f*f1*f1) / (6*f1*f1*f1 - 6*f*f1*f2 + f*f*f3); }
  // NEEDS VERIFICATION. In the unit tests that we have so far, the practical convergence is not 
  // much faster than that of Halley's method - but the numbers are very small anyway (like 3 
  // iterations), so that could be noise/coincidence. More tests are needed.

  /** Implements root finding via the 3rd order Householder method. The function "func" must take 
  the input a first parameter and produce the 0th, 1st, 2nd and 3rd derivative in the following 
  output parameters which are passed by pointer. */
  static T householder3(const std::function<void(T, T*, T*, T*, T*)>& func, T xGuess, T y = 0);


};


#endif
