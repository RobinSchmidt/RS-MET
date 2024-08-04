//-------------------------------------------------------------------------------------------------
// positive range saturation functions:

template<class T>
T rsPositiveSigmoids<T>::linear(T x)
{
  if(x > 2.0)
    return 1.0;
  else
    return 0.5*x;
}

template<class T>
T rsPositiveSigmoids<T>::rational(T x)
{
  return x / (1+x);
}

template<class T>
T rsPositiveSigmoids<T>::cubicRational(T x)
{
  x *= 1 + x + x*x;    //  x + x^2 + x^3
  return x / (x+1);    // (x + x^2 + x^3) / (x + x^2 + x^3 + 1)
}

template<class T>
T rsPositiveSigmoids<T>::cubic(T x)
{
  // The coefficient for the cubic term. If we wanted f'(0)=k, we'd get a3 = -k*(k-3)^2/27 which
  // reduces to -4/27 for k=1:
  static const T a3 = T(-4) / T(27);

  if(x > 1.5)
    return 1.0;
  else
    return x + a3*x*x*x;
}

template<class T>
T rsPositiveSigmoids<T>::quartic(T x)
{
  if(x > 2.0)
    return 1.0;
  else
    return T(0.0625*x*((x-4)*x*x + 16)); // y = x - x^3/4 + x^4/16
}

template<class T>
T rsPositiveSigmoids<T>::hexic(T x)
{
  if(x > 2.0)
    return 1.0;
  else
  {
    T x2 = x*x;    // x^2
    T x4 = x2*x2;  // x^4
    return  T(x - 0.3125*x4 + 0.1875*x4*x - 0.03125*x4*x2);
  }
}

template<class T>
T rsPositiveSigmoids<T>::softClipHexic(T x, T t)
{
  if(x <= t)
    return x;
  else
    return t + (1-t) * hexic((x-t)/(1-t));
}

template<class T>
T rsPositiveSigmoids<T>::softClipHexic(T x)
{
  if(x <= 0.5)
    return x;
  else
    return T(0.5 + 0.5 * hexic(2*(x-(T)0.5)));
}

template<class T>
T rsPositiveSigmoids<T>::invRational(T x)
{
  T p2 = T(1) / (T(2)*x);              //  p/2 in pq-formula
  return -p2 + rsSqrt(p2*p2 + T(1));   // -p/2 + sqrt( (p/2)^2 - q )
}

// saturation polynomials:
//
// f(x) = k*x + a_3*x^3
// f(0)=0,f'(0)=k,f''(0)=0 satisfied by design
// f'(s)=0,f(s)=1 additionally required
// 2 equations for 2 unknowns s, a_3:
// 0 = k + 3*a_3*s^2,
// 1 = k*s + a_3*s^3
// solution: s = 3/(2*k), a_3 = -(4*k^3)/27
//
// f(x) = k*x + a_3*x^3 + a_4*x^4
// f(0)=0,f'(0)=k,f''(0)=0 satisfied by design
// f'(s)=0,f(s)=1,f''(s)=0 additionally required
// 3 equations for 3 unknowns s, a_3, a_4:
// 0 = k + 3*a_3*s^2 + 4*a_4*s^3,
// 1 = k*s + a_3*s^3 + a_4*s^4,
// 0 = 6*a_3*s + 12*a_4*s^2
// solution: s = 2/k, a_3 = -k^3/4, a_4 = k^4/16
//
// f(x) = k*x + a_4*x^4 + a_5*x^5 + a_6*x^6
// f(0)=0,f'(0)=k,f''(0)=0,f'''(0)=0 satisfied by design
// f'(s)=0,f(s)=1,f''(s)=0,f'''(s)=0 additionally required
// 4 equations for 4 unknowns s, a_4, a_5, a_6:
// 0 = k + 4*a_4*s^3 + 5*a_5*s^4 + 6*a_6*s^5,
// 1 = k*s + a_4*s^4 + a_5*s^5 + a_6*s^6,
// 0 = 12*a_4*s^2 + 20*a_5*s^3 + 30*a_6*s^4,
// 0 = 24*a_4*s + 60*a_5*s^2 + 120*a_6*s^3
// ...wolfram alpha doesn't understand your query, even if we set k=1 - maybe the s^6 terms are too
// much, i.e. no analytic solution exists

// what's also good is something like y = clip(x + a*x^N, -1, +1) where N is an odd number and a is
// chosen such that the maximum is at y = 1, the higher the exponent, the more linear the range 
// around the origin for example (factors have been eyballed)
// y = x - 0.035*x^11, x - 0.043*x^9, x - 0.057*x^7, x - 0.082*x^5, x - 0.148*x^3
// maybe N should be (2^k)+1 for fast evaluation by successive squaring, x^5 or x^9 seems to look 
// best
// maybe, we can find a formula for a, given N? 
// we have f(x) = x + a*x^N, f'(x) = 1 + N*a*x^(N-1) and we want to find a such that f(x) = 1 when
// f'(x) = 0. To find the peak p, we set f'(p) = 0 = 1 + N*a*p^(N-1), solving this for p gives
// p = (-1/(N*a))^(1/(N-1)) -> put this in f(p) = 1 ....
// (-1/(N*a))^(1/(N-1)) + a * (-1/(N*a))^(N/(N-1)) = 1 ...solve for a
// nah - it's simpler: solve f(p) = 1 for a and use that in f'(p) = 0, it gives
// p = N/(N-1), a = (1-p)/p^N
// oh - and we need to hardclip BEFORE the polynomial at x=p, not after at 1!!!
// ok - implemented in rsNormalizedSigmoids<T>::clippedOddPower ...needs testing


// maybe give other options for saturation-functions: sin, atan, x/(1+x), etc.

//-------------------------------------------------------------------------------------------------
// normalized, symmetric saturation functions:

template<class T>
T rsNormalizedSigmoids<T>::clip(T x)
{
  if(x < -1.0)
    return -1.0;
  if(x > +1.0)
    return +1.0;
  return x;
}
// todo: do it branch free via 0.5 * (abs(x+1) - abs(x-1))
// https://www.desmos.com/calculator/4hkxkw7xnv

template<class T>
T rsNormalizedSigmoids<T>::atan(T x)
{
  return (T) (::atan(0.5*PI*x) / (0.5*PI));  // optimize: precompute PI/2 and 1/(PI/2)
}

template<class T>
T rsNormalizedSigmoids<T>::tanh(T x)
{
  return ::tanh(x);
  //return rsTanh(x); // use the exp-based version later (more efficient), whe the real functions 
                    // have been added
}

template<class T>
T rsNormalizedSigmoids<T>::clippedOddPower(T x, int N)
{
  T p = T(N) / T(N-1);
  T a = (T(1)-p) / rsPow(p, N);
  x = rsClip(x, -p, p);
  return x - rsPow(x, N);
  // this is a very inefficient prototype - todo: choose some power, like 5 or 9 (if N=(2^k)-1, the
  // power can be evaluated by repeated squaring) and precompute the p,a coeffs
}

template<class T>
T rsNormalizedSigmoids<T>::powRatio(T x, T p)
{
  T tmp = pow(fabs(x), p);
  if(tmp == RS_INF(T))
    return rsSign(x);
  return x * pow(1 + tmp, -1/p);
}

// symmetrized positive-range functions (some boilerplate code):
// maybe try rsSign(x) * positiveSigmoid(fabs(x)) - make performance test

template<class T>
T rsNormalizedSigmoids<T>::rational(T x)
{
  return x / (1+rsAbs(x));
}

template<class T>
T rsNormalizedSigmoids<T>::cubicRational(T x)
{
  return rsSign(x) * rsPositiveSigmoids<T>::cubicRational(rsAbs(x));
}

template<class T>
T rsNormalizedSigmoids<T>::cubic(T x)
{
  if(x >= 0.0) return  rsPositiveSigmoids<T>::cubic( x);
  else         return -rsPositiveSigmoids<T>::cubic(-x);
}

template<class T>
T rsNormalizedSigmoids<T>::quartic(T x)
{
  if(x >= 0.0) return  rsPositiveSigmoids<T>::quartic( x);
  else         return -rsPositiveSigmoids<T>::quartic(-x);
}

template<class T>
T rsNormalizedSigmoids<T>::hexic(T x)
{
  if(x >= 0.0) return  rsPositiveSigmoids<T>::hexic( x);
  else         return -rsPositiveSigmoids<T>::hexic(-x);
}

template<class T>
T rsNormalizedSigmoids<T>::softClipHexic(T x)
{
  if(x >= 0.0) return  rsPositiveSigmoids<T>::softClipHexic( x);
  else         return -rsPositiveSigmoids<T>::softClipHexic(-x);
}

//-------------------------------------------------------------------------------------------------
// class rsParametricSigmoid:

template<class T>
rsParametricSigmoid<T>::rsParametricSigmoid()
{
  y1 = 0.75;         // value y at x=1
  yb = 0.75;         // breakpoint for y1, above which we switch to piecewise function
  computeCoeffs();
}

template<class T>
void rsParametricSigmoid<T>::setValueAt1(T newValue)
{
  y1 = newValue;
  computeCoeffs();
}

template<class T>
void rsParametricSigmoid<T>::setThreshold(T newThreshold)
{
  setValueAt1(newThreshold * (T)0.25 + (T)0.75);
  // nope - formula is wrong
}

template<class T>
void rsParametricSigmoid<T>::setPiecewiseBreakpoint(T newBreakpoint)
{
  yb = newBreakpoint;
  computeCoeffs();
}

template<class T>
T rsParametricSigmoid<T>::coreFunction(T x, T a, T b)
{
  T t = x*x;                   // t = x^2
  t = x + a*(b*t + (1-b)*t*x); // t = x + a*(b*x^2 + (1-b)*x^3) 
  return t / (t+1);            //    (x + a*(b*x^2 + (1-b)*x^3)) / (x + a*(b*x^2 + (1-b)*x^3) + 1)
}

template<class T>
T rsParametricSigmoid<T>::getA(T y1)
{
  return (1-2*y1)/(y1-1);
}

template<class T>
T rsParametricSigmoid<T>::getB(T a)
{
  //return 1 / (a+1);
  return 1 / rsMax((T)1, a);

  // The formula: b = 1 / (a+1) is motivated as follows: the 2nd derivative (curvature) of the 
  // core function f at the origin x=0 is given by c = 2*a*b - 2. It would seem desirable to set the 
  // curvature to 0 at x=0 which would result in b = 1/a. However, we need to restrict b to
  // b <= 1 (otherwise the function becomes unbounded), that's why 1 is added in the denominator (a
  // may go down to 0). Another option would be to use b = min(1, 1/a) or to expose the desired 
  // curvature at the origin as user parameter and then use b = min(1, (2+c)/(2*a). Maybe this can be
  // explored further at some stage...

  // Update: i think, b = 1 / rsMax(1.0, a) is better. Then, the piecewise function will have 
  // matched 1st and 2nd derivative at the junction. But we need yb >= 0.75 then to ensure the 
  // function to be below the identity function everywhere (i have no formula for this, just
  // found this value by trial and error)

  // So, that also means, when we set y1 >= yb = 0.75, the function becomes 2nd order continuous
  // for y1 < 0.75, it's only 1st order continuous

  // Intuitive interpretation: we minimize the 2nd derivative at the origin subject to the 
  // constraint that the function doesn't bulge beyond the identity function. With yb = 0.75, as 
  // soon as the 2nd derivative reaches 0, we enter the regime of a piecewise defined function.
  // ...at least, i think so - anyway, the results look good.

  // I'm not 100% sure, but i think it behaves as follows:
  // Assume yb = 0.75. When y1 = 0.5, the 2nd derivative at the origin jumps from 2 to -2. As we 
  // increase y1, this jump gets smaller and smaller until it disappears (the 2nd derivative 
  // reaches 0) which happens at y1 = 0.75. When increasing y1 further, we insert an appropriate
  // segment of the identity function from -ty to +ty, keeping the 2nd derivative zero, and 
  // scaling/shifting the core function appropriately. If yb would be > 0.75, we would get a 
  // nonzero 2nd derivative when increasing y1 > 0.75 (but i think, without jump) - take this
  // with a grain of salt - i didn't really verify this.
}

template<class T>
void rsParametricSigmoid<T>::computeCoeffs()
{
  if(y1 > yb)
  {
    a  = getA(yb);
    b  = getB(a);
    ty = (y1 - coreFunction(1, a, b)) / (1 - yb);
    sy = 1-ty;
    if(sy >= RS_EPS(T))
      sx = 1/sy;
    else
      sx = 1.0;
  }
  else
  {
    a  = getA(y1);
    b  = getB(a);
    ty = 0;
    sy = 1;
    sx = 1;
  }
  c2 = a*b;
  c3 = a*(1-b);
}

//-------------------------------------------------------------------------------------------------
// class ScaledAndShiftedSigmoid:

template<class T>
void rsScaledAndShiftedSigmoid<T>::setCenter(T newCenter)
{
  center = newCenter;
  updateCoeffs();
}

template<class T>
void rsScaledAndShiftedSigmoid<T>::setWidth(T newWidth)
{
  width = newWidth;
  updateCoeffs();
}

template<class T>
void rsScaledAndShiftedSigmoid<T>::setSigmoid(T (*newSigmoid)(T))
{
  sigmoid = newSigmoid;
}

template<class T>
void rsScaledAndShiftedSigmoid<T>::updateCoeffs()
{
  //rsRangeConversionCoefficients(center-T(0.5)*width, center+T(0.5)*width, T(-1), T(+1), 
  //  &scaleX, &shiftX);
  // This is numerically not good to use rsRangeConversionCoefficients with 
  // center-T(0.5)*width, center+T(0.5)*width because inside the function, we will subtract
  // two very close numbers if width is small. this results in precision loss. We should drag
  // the function body in here and simplify the formulas for better numeric precision..
  // ..done?

  //*scale = (outMax-outMin) / (inMax-inMin); // outMax=1, outMin=-1
  //*shift = outMin - (*scale * inMin);

  scaleX =  2 / width;
  shiftX = -1 - (scaleX * (center-T(0.5)*width));
  scaleY =  1 / scaleX; 
  shiftY = -shiftX * scaleY;
}

/*
ToDo: implement the smooth-transition function from here:

https://en.wikipedia.org/wiki/Non-analytic_smooth_function

defined as
  g(x) = f(x) / (f(x) - f(1-x)) where f(x) = exp(-1/x) for x > 0 and 0 for x <= 0
f(x) is the standard example for a smooth, non-analytic function. A smooth, multi-dimensional bump 
function can also be constructed form f as:
  h(x) = f( r^2 - |x|^2 )
where |x| is the Euclidean norm of the vector x and r is the radius of the bump


Here's another interesting cheap sigmoid:
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=521377


May also be relevant:
https://en.wikipedia.org/wiki/Smoothstep


https://www.youtube.com/watch?v=60VoL-F-jIQ Smoothstep: The most useful function
https://www.youtube.com/watch?v=YJ4iyff7zbk&list=PLGmrMu-IwbgsY3onv9rrzHvm7OpG43Uvk Function Smoothing Explained!


Implement also:
https://en.wikipedia.org/wiki/Gudermannian_function  f(x) = (2/pi) * atan(sinh((pi/2)*x))
See also:
https://www.youtube.com/watch?v=oIChUOV_0w4  at 5:55


Maybe make a similar file for functions of the "1-cycle", "wobble" kind. Examples:

x * exp(-x^2), x - x*tanh(x^2), ...

We can construct it by multiplying a bell shape with x, taking a derivative of a bell, subtracting
x*sigmoid(x^2) from x, ... See:

https://www.desmos.com/calculator/qnueorgb0z

Maybe define the length as two times the length between minimum and maximum. These shapes could be 
useful for some sort of time-scale transform.

The inverse function of f(x) = x + x^3:
https://www.wolframalpha.com/input?i=inverse+function+of+x+%2B+x%5E3
could be useful as a function that has an s-shape but is unbounded, but grows slowly - namely, 
asymptotically like the cube-root. Maybe such a behavior can sometimes be useful for stabilizing
recursive systems without setting hard saturation limits. asinh would grow even slower 
(logarithmically). We could also consider inverses of x + x^n where n is any odd integer


Implemented as invRational:

The simplest function that I could think of to map the interval [-1,+1] to [-inf,+inf] is
y = -x / ((x-1)*(x+1)) = -x / (x^2 - 1), see https://www.desmos.com/calculator/1h1evzhyay
I just created a rational function with a zero at 0 and two poles at +-1. How about using the 
inverse of that as sigmoid? We just need to solve it for x in terms of y:
y = -x / (x^2 - 1)  ->  y * (x^2 - 1) = -x  ->  y*x^2 + x - y = 0  ->  x^2 + x/y - 1 = 0
...use pq-formula with p = 1/y, q = -1  ->  x = -(1/(2y)) +- sqrt( (1/(4y^2)) + 1 )
we need the +branch for positive y and the -branch for negative y
see https://www.desmos.com/calculator/yxlnijspzg  
comparison with other sigmoids: https://www.desmos.com/calculator/iiml3yqzf6
...it approaches +-1 more slowly. How about raising x to an odd power? In the inverse, we could
use the same formula and at the end, extract the n-th root. How about using poles and/or zeros
with multiplicities, like  y = -x^k / ((x-1)^m * (x+1)^n) ? This way, we coul pehaps make
asymmetric sigmoids. 
Here is it with the branch implemented via sign: sign(x) * sqrt(1 + 1/(2x)^2) - 1/(2x)
https://www.desmos.com/calculator/ziw32b5ovy
A variation: (sign(x) * sqrt(a + 1/(2x)^2) - 1/(2x)) / a
https://www.desmos.com/calculator/7mhot2u03b
with a = 1, we recover the old one.

How about considering the function  1 / (1 - sigmoid(x))  and using its growth behavior to
categorize sigmoids? The idea is that the decay of  1 - sigmoid(x)  mesures how quickly the
sigmoid approaches one. Using its reciprocal therefore expressed this approach in a way in 
which a faster growth means a quicker approach. Maybe we can then design sigmoids with desired
approach-rates? I think, the asymptotic approach rate by which a sigmoid approaches one could be
a meaningful measurement on a sigmoid. Is it linear, quadratic, exponential, etc.?

Looks like the approach-rate of the invRat function is 2*(x + 1/4), see:
https://www.desmos.com/calculator/iugvzyoeqs
...so it's indeed linear. Btw: the approach-rate function itself looks like a smooth version of
the ReLU function - so maybe it could be useful as activation function in neural networks. See:
https://www.desmos.com/calculator/mb982w6mkw

Variations:
y = -x    / ((x  -1)*(x  +1))
y = -x    / ((x^3-1)*(x^3+1))
y = -x^3  / ((x^3-1)*(x^3+1))
y = -4x^3 / ((x^3-1)*(x^3+1))
y = -x^3  / ((x-1)^3*(x+1)^3)



The approach-rate for tanh looks more like an exponential - but not quite. The approach-rate of
x / sqrt(1+x^2) looks kinda parabolic - a bit like (1.3 x + 1)^2 + 1. Try to find more accurate
asymptotic expressions in a systematic way - these were found by trial and error. Maybe to find
an asymptotic expression for some function f(x), we could to a Taylor expansion of f(1/x) around 0?
See: 
https://en.wikipedia.org/wiki/Asymptotology
https://en.wikipedia.org/wiki/Asymptotic_analysis
https://en.wikipedia.org/wiki/Asymptotic_expansion



Maybe implement these in C++
https://easings.net/
Clicking on the graphics will lead to a page where javascript code is given for the functions
Ah - here is the code, I guess:
https://github.com/ai/easings.net
OK - I have copied all the code into Notes/FadeInOut.txt



*/