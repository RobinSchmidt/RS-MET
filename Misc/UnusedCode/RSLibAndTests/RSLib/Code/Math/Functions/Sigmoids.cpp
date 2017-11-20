//-------------------------------------------------------------------------------------------------
// positive range saturation functions:

double rsPositiveSigmoids::linear(double x)
{
  if(x > 2.0)
    return 1.0;
  else
    return 0.5*x;
}

double rsPositiveSigmoids::rational(double x)
{
  return x / (1+x);
}

double rsPositiveSigmoids::cubicRational(double x)
{
  x *= 1 + x + x*x;    //  x + x^2 + x^3
  return x / (x+1);    // (x + x^2 + x^3) + (x + x^2 + x^3 + 1)
}

double rsPositiveSigmoids::cubic(double x)
{
  // The coefficient for the cubic term. If we wanted f'(0)=k, we'd get a3 = -k*(k-3)^2/27 which
  // reduces to -4/27 for k=1:
  static const double a3 = -4.0/27.0;

  if(x > 1.5)
    return 1.0;
  else
    return x + a3*x*x*x;
}

double rsPositiveSigmoids::quartic(double x)
{
  if(x > 2.0)
    return 1.0;
  else
    return 0.0625*x*((x-4)*x*x + 16); // y = x - x^3/4 + x^4/16
}

double rsPositiveSigmoids::hexic(double x)
{
  if(x > 2.0)
    return 1.0;
  else
  {
    double x2 = x*x;    // x^2
    double x4 = x2*x2;  // x^4
    return x - 0.3125*x4 + 0.1875*x4*x - 0.03125*x4*x2;
  }
}

double rsPositiveSigmoids::softClipHexic(double x, double t)
{
  if(x <= t)
    return x;
  else
    return t + (1-t) * hexic((x-t)/(1-t));
}

double rsPositiveSigmoids::softClipHexic(double x)
{
  if(x <= 0.5)
    return x;
  else
    return 0.5 + 0.5 * hexic(2*(x-0.5));
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

double rsNormalizedSigmoids::clip(double x)
{
  if(x < -1.0)
    return -1.0;
  if(x > +1.0)
    return +1.0;
  return x;
}

double rsNormalizedSigmoids::atan(double x)
{
  return ::atan(0.5*PI*x) / (0.5*PI);  // optimize: precompute PI/2 and 1/(PI/2)
}

double rsNormalizedSigmoids::tanh(double x)
{
  return rsTanh(x);
}

double rsNormalizedSigmoids::powRatio(double x, double p)
{
  double tmp = pow(fabs(x), p);
  if(tmp == RS_INF(double))
    return rsSign(x);
  return x * pow(1 + tmp, -1/p);
}

// symmetrized positive-range functions (some boilerplate code):
// maybe try rsSign(x) * positiveSigmoid(fabs(x)) - make performance test

double rsNormalizedSigmoids::rational(double x)
{
  return x / (1+rsAbs(x));
}

double rsNormalizedSigmoids::cubicRational(double x)
{
  return rsSign(x) * rsPositiveSigmoids::cubicRational(rsAbs(x));
}

double rsNormalizedSigmoids::cubic(double x)
{
  if(x >= 0.0)
    return rsPositiveSigmoids::cubic(x);
  else
    return -rsPositiveSigmoids::cubic(-x);
}

double rsNormalizedSigmoids::quartic(double x)
{
  if(x >= 0.0)
    return rsPositiveSigmoids::quartic(x);
  else
    return -rsPositiveSigmoids::quartic(-x);
}

double rsNormalizedSigmoids::hexic(double x)
{
  if(x >= 0.0)
    return rsPositiveSigmoids::hexic(x);
  else
    return -rsPositiveSigmoids::hexic(-x);
}

double rsNormalizedSigmoids::softClipHexic(double x)
{
  if(x >= 0.0)
    return rsPositiveSigmoids::softClipHexic(x);
  else
    return -rsPositiveSigmoids::softClipHexic(-x);
}

//-------------------------------------------------------------------------------------------------
// class rsParametricSigmoid:

rsParametricSigmoid::rsParametricSigmoid()
{
  y1 = 0.75;         // value y at x=1
  yb = 0.75;         // breakpoint for y1, above which we switch to piecewise function
  computeCoeffs();
}

void rsParametricSigmoid::setValueAt1(double newValue)
{
  y1 = newValue;
  computeCoeffs();
}

void rsParametricSigmoid::setThreshold(double newThreshold)
{
  setValueAt1(newThreshold * 0.25 + 0.75);
  // nope - formula is wrong
}

void rsParametricSigmoid::setPiecewiseBreakpoint(double newBreakpoint)
{
  yb = newBreakpoint;
  computeCoeffs();
}

double rsParametricSigmoid::coreFunction(double x, double a, double b)
{
  double t = x*x;              // t = x^2
  t = x + a*(b*t + (1-b)*t*x); // t = x + a*(b*x^2 + (1-b)*x^3) 
  return t / (t+1);            //    (x + a*(b*x^2 + (1-b)*x^3)) / (x + a*(b*x^2 + (1-b)*x^3) + 1)
}

double rsParametricSigmoid::getA(double y1)
{
  return (1-2*y1)/(y1-1);
}

double rsParametricSigmoid::getB(double a)
{
  //return 1 / (a+1);
  return 1 / rsMax(1.0, a);

  // The formula: b = 1 / (a+1) is motivated as follows: the 2nd derivative (curvature) of the 
  // core function f at the origin x=0 is given by c = 2*a*b - 2. It would seem desirable to set the 
  // curvature to 0 at x=0 which would result in b = 1/a. However, we need to restrict b to
  // b <= 1 (otherwise the function becomes unbouded), that's why 1 is added in the denominator (a
  // may go down to 0) . Another option would be to use b = min(1, 1/a) or to expose the desired 
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

void rsParametricSigmoid::computeCoeffs()
{
  if(y1 > yb)
  {
    a  = getA(yb);
    b  = getB(a);
    ty = (y1 - coreFunction(1, a, b)) / (1 - yb);
    sy = 1-ty;
    if(sy >= RS_EPS(double))
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