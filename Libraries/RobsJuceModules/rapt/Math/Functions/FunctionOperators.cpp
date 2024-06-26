

/*=================================================================================================


Notes:

These functionals can be helpful to create new functions with certain desried properties from known
functions, i.e. they can be used as building blocks in the construction of functions.

Basic transformations:

  g(x) = f(x-a)                        Shift to the right by a 
  g(x) = f(a*x)                        Squeeze by factor a in x-direction
  g(x) = f(x) + a                      Shift upward by a
  g(x) = a*f(x)                        Stretch by factor a in y-direction
  g(x) = f(-x)                         Reflect along y-axis
  g(x) = -f(x)                         Reflect along x-axis
  g(x) = -f(-x)                        Reflect about origin
  g(x) = (f(x) + f(-x)) / 2            Extract symmetric (aka even) part
  g(x) = (f(x) - f(-x)) / 2            Extract antisymmetric (aka odd) part
  g(x) = 1/f(x)                        Reciprocal function
  g(x) = f(1/x)                        Turns the function "inside out", x=0 maps to x=inf etc.
  g(x) = 1/f(1/x)                      Reciprocal of inside-out
  g(x) = f^-1 (x)                      Inverse function
  g(x) = f^n (x)                       Iterated function (n times)
  g(x) = f^(n) (x)                     n-th derivative (notation a bit ambiguous)
  g(x) = f(1-x)                        Mirror the unit interval 0..1 on input side
  g(x) = 1-f(x)                        Mirror the unit interval 0..1 on output side


Functionals involving a fixed sigmoid function s(x) with range 0..1 like the logistic function:

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

- g(x) = (f(x-a) + f(x+a)) / 2

- To transfrom a function by a matrix, see the experiment in functionOperators() in the test 
  project. The difficult part is that we want to express the resulting graph again as a function
  rather than an implicit or parametric equation. In the latter case, case we can just explicitly 
  produce (x,y)-pairs and then transform them explicitly by our matrix. But when we want to 
  produce a transformed function g(x), we are not allowed to transform the x-variable. It is given 
  to us and fixed and we must find the point on the *transformed* graph with corresponding 
  y-coordinate. This might not even be possible and if it is, the result may not be unique. 
  ...TBC..., See:
  https://math.stackexchange.com/questions/17246/is-there-a-way-to-rotate-the-graph-of-a-function
  https://www.youtube.com/watch?v=h9OWnuarYuc

- We could transform a function by a matrix using something like the class rsInterpolatingFunction.
  We can transform all the (x,y)-pairs in such an object - but the trafo may destroy the property
  that the x-values are strictly monotonically increasing which is essential for the look-up table
  like mode of operation. We might clean this up by removing the problematic portions of the 
  resulting non-function to turn it into a proper function again. Such a feature may actually be 
  useful in FuncShaper. Maybe give the user some ways of constructing a transform matrix - maybe in 
  terms of a pre-rotation -> axis-parallel scalings for x,y -> post-rotation. But that feature 
  would then also have to apply the "clean-up" algorithm to the result. The algo could work as 
  follows:
  -Produce the raw (x,y) pairs
  -Transform all the pairs via the matrix
  -If x[0] > x[N-1], reverse x- and y-arrays
  -Go through x-array and if x[n] <= x[n-1], remove the pair (x[n],y[n]). In practice, this should
   be optimized by not removing the pairs one by one but by first determining portions to be 
   removed and then removing them in one go (to avoid excessive copying around of values)
  I think, in the context of rsInterpolating function, such a "removal" would correspond to replace 
  the portion with an interpolant between the endpoints of the "removed" section - which seems to 
  be reasonable behavior. Or maybe it would lead to a jump discontinuity at the end of the section?
  Yeah - guess, that's what's going to happen. But maybe a more sophisticated algorithm could give
  rise to the (perhaps more desirable) "replace with interpolant" behavior.


*/