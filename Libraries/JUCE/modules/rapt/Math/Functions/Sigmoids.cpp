//-------------------------------------------------------------------------------------------------
// positive range saturation functions:

template<class T>
T PositiveSigmoids<T>::linear(T x)
{
  if(x > 2.0)
    return 1.0;
  else
    return 0.5*x;
}

template<class T>
T PositiveSigmoids<T>::rational(T x)
{
  return x / (1+x);
}

template<class T>
T PositiveSigmoids<T>::cubicRational(T x)
{
  x *= 1 + x + x*x;    //  x + x^2 + x^3
  return x / (x+1);    // (x + x^2 + x^3) / (x + x^2 + x^3 + 1)
}

template<class T>
T PositiveSigmoids<T>::cubic(T x)
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
T PositiveSigmoids<T>::quartic(T x)
{
  if(x > 2.0)
    return 1.0;
  else
    return T(0.0625*x*((x-4)*x*x + 16)); // y = x - x^3/4 + x^4/16
}

template<class T>
T PositiveSigmoids<T>::hexic(T x)
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
T PositiveSigmoids<T>::softClipHexic(T x, T t)
{
  if(x <= t)
    return x;
  else
    return t + (1-t) * hexic((x-t)/(1-t));
}

template<class T>
T PositiveSigmoids<T>::softClipHexic(T x)
{
  if(x <= 0.5)
    return x;
  else
    return T(0.5 + 0.5 * hexic(2*(x-(T)0.5)));
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

// maybe give other options for saturation-functions: sin, atan, x/(1+x), etc.

//-------------------------------------------------------------------------------------------------
// normalized, symmetric saturation functions:

template<class T>
T NormalizedSigmoids<T>::clip(T x)
{
  if(x < -1.0)
    return -1.0;
  if(x > +1.0)
    return +1.0;
  return x;
}

template<class T>
T NormalizedSigmoids<T>::atan(T x)
{
  return (T) (::atan(0.5*PI*x) / (0.5*PI));  // optimize: precompute PI/2 and 1/(PI/2)
}

template<class T>
T NormalizedSigmoids<T>::tanh(T x)
{
  return ::tanh(x);
  //return rsTanh(x); // use the exp-based version later (more efficient), whe the real functions 
                    // have been added
}

template<class T>
T NormalizedSigmoids<T>::powRatio(T x, T p)
{
  T tmp = pow(fabs(x), p);
  if(tmp == RS_INF(T))
    return rsSign(x);
  return x * pow(1 + tmp, -1/p);
}

// symmetrized positive-range functions (some boilerplate code):
// maybe try rsSign(x) * positiveSigmoid(fabs(x)) - make performance test

template<class T>
T NormalizedSigmoids<T>::rational(T x)
{
  return x / (1+rsAbs(x));
}

template<class T>
T NormalizedSigmoids<T>::cubicRational(T x)
{
  return rsSign(x) * PositiveSigmoids<T>::cubicRational(rsAbs(x));
}

template<class T>
T NormalizedSigmoids<T>::cubic(T x)
{
  if(x >= 0.0)
    return PositiveSigmoids<T>::cubic(x);
  else
    return -PositiveSigmoids<T>::cubic(-x);
}

template<class T>
T NormalizedSigmoids<T>::quartic(T x)
{
  if(x >= 0.0)
    return PositiveSigmoids<T>::quartic(x);
  else
    return -PositiveSigmoids<T>::quartic(-x);
}

template<class T>
T NormalizedSigmoids<T>::hexic(T x)
{
  if(x >= 0.0)
    return PositiveSigmoids<T>::hexic(x);
  else
    return -PositiveSigmoids<T>::hexic(-x);
}

template<class T>
T NormalizedSigmoids<T>::softClipHexic(T x)
{
  if(x >= 0.0)
    return PositiveSigmoids<T>::softClipHexic(x);
  else
    return -PositiveSigmoids<T>::softClipHexic(-x);
}

//-------------------------------------------------------------------------------------------------
// class rsParametricSigmoid:

template<class T>
ParametricSigmoid<T>::ParametricSigmoid()
{
  y1 = 0.75;         // value y at x=1
  yb = 0.75;         // breakpoint for y1, above which we switch to piecewise function
  computeCoeffs();
}

template<class T>
void ParametricSigmoid<T>::setValueAt1(T newValue)
{
  y1 = newValue;
  computeCoeffs();
}

template<class T>
void ParametricSigmoid<T>::setThreshold(T newThreshold)
{
  setValueAt1(newThreshold * (T)0.25 + (T)0.75);
  // nope - formula is wrong
}

template<class T>
void ParametricSigmoid<T>::setPiecewiseBreakpoint(T newBreakpoint)
{
  yb = newBreakpoint;
  computeCoeffs();
}

template<class T>
T ParametricSigmoid<T>::coreFunction(T x, T a, T b)
{
  T t = x*x;                   // t = x^2
  t = x + a*(b*t + (1-b)*t*x); // t = x + a*(b*x^2 + (1-b)*x^3) 
  return t / (t+1);            //    (x + a*(b*x^2 + (1-b)*x^3)) / (x + a*(b*x^2 + (1-b)*x^3) + 1)
}

template<class T>
T ParametricSigmoid<T>::getA(T y1)
{
  return (1-2*y1)/(y1-1);
}

template<class T>
T ParametricSigmoid<T>::getB(T a)
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
void ParametricSigmoid<T>::computeCoeffs()
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
void ScaledAndShiftedSigmoid<T>::setCenter(T newCenter)
{
  center = newCenter;
  updateCoeffs();
}

template<class T>
void ScaledAndShiftedSigmoid<T>::setWidth(T newWidth)
{
  width = newWidth;
  updateCoeffs();
}

template<class T>
void ScaledAndShiftedSigmoid<T>::setPrototypeSigmoid(T (*newSigmoid)(T))
{
  sigmoid = newSigmoid;
}

template<class T>
void ScaledAndShiftedSigmoid<T>::updateCoeffs()
{
  //rsRangeConversionCoefficients(center-T(0.5)*width, center+T(0.5)*width, T(-1), T(+1), 
  //  &scaleX, &shiftX);
  // This is numerically not good to use rsRangeConversionCoefficients with 
  // center-T(0.5)*width, center+T(0.5)*width because inside the function, we will subtract
  // two very close numbers if width is small. this results in precision loss. We should drag
  // the function body in here and simplify the formulas for better numeric precision..

  //*scale = (outMax-outMin) / (inMax-inMin); // outMax=1, outMin=-1
  //*shift = outMin - (*scale * inMin);

  scaleX =  2 / width;
  shiftX = -1 - (scaleX * (center-T(0.5)*width));
  scaleY =  1 / scaleX; 
  shiftY = -shiftX * scaleY;
}