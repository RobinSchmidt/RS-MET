
template<class T>
std::vector<T> rsOdeCoeffs<T>::adamsBashforth(int s)
{
  using Poly = rsPolynomial<T>;
  std::vector<T> b(s);
  T sign = T(1);                        // for the (-1)^j
  for(int j = 0; j < s; j++) {
    Poly p({ T(1) });
    for(int i = 0; i < s; i++) 
      if(i != j)
        p = p * Poly({T(i), T(1)});
    T d = p.definiteIntegral(T(0), T(1));
    b[s-j-1] = (sign*d) / (rsFactorial(j) * rsFactorial(s-j-1));
    sign *= T(-1);  }
  return b;
}

template<class T>
std::vector<T> rsOdeCoeffs<T>::adamsMoulton(int s)
{
  using Poly = rsPolynomial<T>;
  std::vector<T> b(s+1);
  T sign = T(1);
  for(int j = 0; j <= s; j++) {
    Poly p({ T(1) });
    for(int i = 0; i <= s; i++) 
      if(i != j)
        p = p * Poly({T(i-1), T(1)}); 
    T d = p.definiteIntegral(T(0), T(1));
    b[s-j] = (sign*d) / (rsFactorial(j) * rsFactorial(s-j));
    sign *= T(-1); }
  return b;
}
// Maybe the algos can be turned into an O(N) algos by not creating the polynomial p from 
// scratch leaving out the i=j factor each time but instead construction a "master" polynomial and 
// dividing out the i=j factor in each iteration. Oh, and the factorials could be computed more 
// efficiently on the fly, too. And we should avoid all these pesky re-allocations related to the
// polynomials. We should just have a function adamsMoulton(int order, T* coeffs, T* work) and 
// don't allocate any heap memory. The workspace is required for the polynomials.
// Maybe it's more convenient to return the arrays in reversed order - we'll see. That would just 
// require to replace b[s-j-1] = ... and b[s-j] = ... by b[j] in A-Bashforth and A-Moulton 
// respectively. Depends on which format is most convenient to use in the solvers themselves.



/*

ToDo: 
-implement similar methods to generate coeffs for other types of solvers such as Runge-Kutta:
 https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods 
 https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf


-Maybe have a baseclass that implements explicit solvers, a subclass that implements implicit 
 solvers (maybe using an explicit fomrula for the initial guess which may also serve for error
 estimation using the explicit method as predicttor and the implicit method as corrector.
-Maybe factopr out another class that contains only the explicit single-step methods
-Make subclass with stepsize control...or maybe that should be a driver class that has a pointer
 to the stepper, so it may work with the baseclass or the subclass

-Baseclass: rsOdeSolverBase, contains a virtual method doStep() and/or doStep(Tx dx)
-rsOdeSolverMultistep : rsOdeSolverBase: contains pointer to rsOdeSolver for the initial section 
 (which should be a single-step solver)
-rsOdeSolver should be a convenience class (using the Facade pattern, sort of) which contains all
 the available ODE-solving functionality in a single entity. it should provide functions:
 doStep, getOutput, reset, setMethod, setStepSize(Tx s, bool fixed), setFunction(Func f)

-For the implicit solvers, we should not just use Newton iteration but instead give the user a way
 to supply the root-finder function. And that should also work when y is a vector

-Maybe we should also have a way to solve implicit ODEs not given as y' = F(x,y) but as 
 F(x,y,y') = 0? Don't confuse this sort of implicitness with the implicitness of a solver method!
 Here, it applies to the ODE itself, not the solver. Maybe try as examples things like 
 (x+y)*y' = 0, x*y*y' = 0, (x+y)/y' = 0, (x-y)*y' = 0 etc. - it's important that y' can't be 
 isolated on one side. ..Well, for testing we could actually use one, where it could be isolated 
 but we don't do so and nevertheless express it in implicit form, like y' = x*y becomes 
 x*y - y' = 0. Would that require fundamentally different methods? Or maybe the same methods can be
 used and it's just that the implementation of F will itself have to use an iterative solver? So 
 maybe the API can be prepared for that right from the start by giving F the right type of 
 function. Well, actually, using std::function<void(Tx x, const Ty& y, Ty& yd)> (as it is now) 
 should work perfectly fine. We don't really care here, how yd is computed - whether there's an 
 explicit formula or a root-finder at work (or even something else - maybe a table or some weird 
 algo could also be options). Maybe we should define that when it's a root finder, the passed yd
 should be used an initial guess and we should probably pass the value from the evaluation of the
 previous datapoint. ...or maybe we could give it a better guess by passing a polynomial through 
 our previous yd datapoints and use it to extrapolate? ...but that should be optional because it's
 an extra cost. Maybe we could give the user the option to choose the order of the extrapolation
 polynomial where 0 means to just hold the value from the previous datapoint (the default 
 behavior). Of course, the client coud also record past values of yd himself and do the 
 extrapolation for the guess himself, but we have these yd values available here anyway, so that
 would be (a small amount of) redundant data storage..hmmm...



https://en.wikipedia.org/wiki/Ordinary_differential_equation





*/