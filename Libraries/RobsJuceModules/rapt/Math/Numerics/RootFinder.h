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

References
(1) Numerical Recipies in C (2nd Edition), Chapter 9

*/

template<class T>
class rsRootFinder
{
public:

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

  // todo: let caller pass a tolerance (maybe separate for x and y), state pathological conditions
  // when no convergence can be epxected (like when the function has poles like 1/x)


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


  //static T modifiedFalsePosition(std::function<T(T)>& func, T xLeft, T xRight, T y = 0);

  static T newtonStep(const T& f, const T& fp) { return -f/fp; }
  // can be used like x += newtonStep(f, fp)
  // https://en.wikipedia.org/wiki/Newton%27s_method
  // https://de.wikipedia.org/wiki/Newtonverfahren

  /** Under Construction */
  static T newton(const std::function<void(T, T*, T*)>& func, T xGuess, T y = 0);

  // func should take the x value as 1st parameter and produce value and derivative in 2nd and 3rd
  // parameter respectively



  static T halleyStep(const T& f, const T& f1, const T& f2) { return (2*f*f1) / (f*f2 - 2*f1*f1); }
  // can be used like x += halleyStep(f, f1, f2) where f,f1,f2 are: value, 1st derivative, 
  // 2nd derivative respectively
  // https://en.wikipedia.org/wiki/Halley%27s_method
  // https://de.wikipedia.org/wiki/Halley-Verfahren


  static T halley(const std::function<void(T, T*, T*, T*)>& func, T xGuess, T y = 0);



  static T householder3Step(const T& f, const T& f1, const T& f2, const T& f3) 
  { 
    return (3*f*f*f2 - 6*f*f1*f1) / (6*f1*f1*f1 - 6*f*f1*f2 + f*f*f3);
  }


  static T householder3(const std::function<void(T, T*, T*, T*, T*)>& func, T xGuess, T y = 0);


  // Has formula with 3 derivatives:
  // https://en.wikipedia.org/wiki/Householder%27s_method#Example
  // return (3*f*f*f2 - 6*f*f1*f1) / (6*f1*f1*f1 - 6*f*f1*f2 + f*f*f3);

  // Higher order variants of Newton iteration:
  // http://numbers.computation.free.fr/Constants/Algorithms/newton.html
  // https://tminka.github.io/papers/minka-newton.pdf
  // https://www.researchgate.net/publication/268555974_Beyond_Newton's_Method_Generalized_Higher-Order_Approximation_Methods
  //  -> discusses Chebychev's method as a vector generalization of Halley's method
  // https://en.wikipedia.org/wiki/Householder%27s_method
  // https://archive.org/details/numericaltreatme0000hous

  // secant, newton, ridders, brent, ...

  // see RSLib::UnivariateScalarFunction in file FunctionObjects.h/cpp
  

};


#endif
