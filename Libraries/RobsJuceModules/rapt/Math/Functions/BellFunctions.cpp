//-------------------------------------------------------------------------------------------------
// positive range bell functions:

template<class T>
T rsPositiveBellFunctions<T>::linear(T x)
{
  if(x > 1)
    return 0;
  else
    return T(1) - x;
}

template<class T>
T rsPositiveBellFunctions<T>::cubic(T x)
{
  if(x > 1)
    return 0;
  else
    return 1 + (2*x - 3) * x*x;
}

template<class T>
T rsPositiveBellFunctions<T>::quintic(T x)
{
  if(x > 1)
    return 0;
  else {
    T x2 = x*x;
    return 1 + (-10 + 15*x - 6*x2) * x*x2; // 1 - 10*x^3 + 15*x^4 - 6*x^5
  }
}

template<class T>
T rsPositiveBellFunctions<T>::heptic(T x)
{
  if(x > 1)
    return 0;
  else {
    T x2 = x*x;
    return 1 + (-35 + 84*x - 70*x2 + 20*x2*x) * x2*x2; // 1 - 35*x^4 + 84*x^5 - 70*x^6 + 20*x^7
  }
}

template<class T>
T rsPositiveBellFunctions<T>::bump(T x)
{
  if(x > 1)
    return 0;
  else
    return T(EULER) * exp(T(-1)/(T(1)-x*x));
  // can we find its antiderivative for use as sigmoid?
  // http://www.wolframalpha.com/input/?i=integral+of+exp(-1%2F(1-x%5E2))
  // hmm...only in terms of a series - but maybe that's good enough
  // the derivative is also interesting - maybe for use as a smooth sawtooth wave:
  // http://www.wolframalpha.com/input/?i=derivative+of+exp(-1%2F(1-x%5E2))

  // Derivation of a recursion of for the derivatives of e^(-1/x) 
  // https://www.youtube.com/watch?v=tvCWi1cCqYM
}

template<class T>
T rsPositiveBellFunctions<T>::bump(T x, T p)
{
  if(x > 1)
    return 0;
  else
    return T(EULER) * exp(T(-1)/(T(1) - pow(rsAbs(x), p))); // exp(-1/(1-|x|^p))
}


//-------------------------------------------------------------------------------------------------
// class rsParametricBellFunction:

template<class T>
rsParametricBellFunction<T>::rsParametricBellFunction()
{
  bell   = &rsPositiveBellFunctions<T>::cubic;
  center = 0; 
  setWidth(2);
  setFlatTopWidth(0);
}

template<class T>
void rsParametricBellFunction<T>::setCenter(T newCenter)
{
  center = newCenter; 
}

template<class T>
void rsParametricBellFunction<T>::setWidth(T newWidth)
{
  a = 2 / newWidth;
}

template<class T>
void rsParametricBellFunction<T>::setFlatTopWidth(T newWidth)
{
  flat = newWidth;
  b = 1 / (1 - flat);
  //if(flat != 1)
  //  b = 1 / (1 - flat);
  //else
  //  b = 1;
}

template<class T>
void rsParametricBellFunction<T>::setPrototypeBell(T (*newFunction)(T))
{
  bell = newFunction;
}


/*


Ideas for other (roughly) bell-shaped functions:

f = (sech(x))^a       sech is 1/cosh.
g = exp(-x)           for reference in plot, should approach f for x -> inf when a=1
h = (f-g)/(f+g)       difference between f and g, renormalized by sum
  -> approaches 1/3 for a=1, 1 for a < 1, -1 for a > 1
  https://www.desmos.com/calculator/di1eccbwuk
  https://www.desmos.com/calculator/kwjphmhlcd (better)
  https://www.wolframalpha.com/input?i=derivative+of+tanh%28x%29
  -> f for a=1 is the derivative of tanh and gives a nice, smooth bell shape that falls off like 
  e^(-x) when a=1, i.e. has fatter tails than a Gaussian. The h function might actually be useful 
  as well (maybe it should be shifted up by 1). With the parameter a, we can fade between sigmoid 
  and bell shapes with asymmetric/skewed bells in between. f/g is also interesting. Maybe plot f 
  against a Gaussian bell to see the difference in tail behavior.

Maybe use f(x) = 1 / (cosh(a*x)) and adjust a such that f(x) meets exp(-x^2) exactly at x=1. Or 
normalize f to unit area (maybe that gives the same a-value? Figure out!). That way, we would have
a nice smooth function that can be interpreted as probability density function but with fatter 
tails than the Gaussian. See:  https://www.desmos.com/calculator/er2azb0rkz
The meet at x=1 happens for a = acosh(e) = 1.657...

How can we generalize this sort of functions to exp(-x^n) shapes but to both sides? A simple 
exp(-|x|^n) is possible and is called the Minkowski distribution but it's not necessarily smooth at 
the origin. Maybe (exp(-(+x)^n) + exp(-(-x)^n)) / 2  could work, i.e. just take an arbitrary 
natural exponent for the sinde of the positive x-axis and then take the symmetric part of the 
resulting function. That's how cosh is created from exp. ..oh - but maybe then take the reciprocal 
or something -> figure out! See: https://www.desmos.com/calculator/kthmhd4r5x
Interestingly, the even and odd functions


*/