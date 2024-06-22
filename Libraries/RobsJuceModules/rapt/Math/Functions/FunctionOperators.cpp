

/*=================================================================================================


Notes:

These functionals can be helpful to create new functions with certain desried properties from known
functions, i.e. they can be used as building blocks in the construction of functions.

Basic transformations:

  g(x) = f(x-a)                        Shift to the right by a 
  g(x) = f(a*x)                        Squeeze by factor a in x-directon
  g(x) = f(x) + a                      Shift upward by a
  g(x) = a*f(x)                        Stretch by factor a in y-direction
  g(x) = f(-x)                         Reflect along y-axis
  g(x) = -f(x)                         Reflect along x-axis
  g(x) = -f(-x)                        Reflect about origin
  g(x) = (f(x) + f(-x)) / 2            Extract symmetric (aka even) part
  g(x) = (f(x) - f(-x)) / 2            Extract antisymmetric (aka odd) part


Functionals involving a sigmoid function s(x) with range 0..1 like the logistic function:

  s(x) = 1 / (1 + exp(-x)) = exp(x) / (1 + exp(x)) = 1 - s(-x)


  h(x) = (1-s(x))*f(x) + s(x)*g(x)     Crossfade between f(x) and g(x). At x=0, mix is 50/50.


...TBC...Include things like

g(x) = f(x) * x*sigmoid(x)

crossfades between two functions

Ideas:

- The crossfade could be generalized to 3 functions by shifting the sigmoids for f and g by -1 and
  +1 in the x-direction and using 1 - sum-of-sigmoids for a weighting function in the middle

- f(x) / (1 + f(x)), 1 / f(x), inv(f(x)) - inverse function - maybe locally

- What about g(x) = f(x-a) * f(x-b) * f(x-c) * .... If f is the identity, this should create a 
  polynomial with roots at a,b,c,...

- We can construct rational functions from their desired poles and zeros. For example:
  f(x) = x / ((x-1)*(x+1)) has one zero at x = 0, and two poles at x = -1 and +1.

- Multiplication by a bell function can smoothly crossfade between zero, f(x) and back to zero if
  the bell falls off quickly enough, i.e. more quickly than f increases, if it does increase.


*/