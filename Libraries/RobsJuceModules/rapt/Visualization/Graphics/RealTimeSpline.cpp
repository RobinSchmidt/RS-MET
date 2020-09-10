template<class TCor, class TWgt>
rsRealTimeSpline<TCor, TWgt>::rsRealTimeSpline()
{
  reset();
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::setDotBuffers(TCor* bufX, TCor* bufY, TWgt* bufWeights, 
  int bufLengths)
{
  dotsX = bufX;
  dotsY = bufY;
  dotsW = bufWeights;
  dotBufferLength = bufLengths;
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::updateDotBuffers(TCor inX, TCor inY)
{
  updatePointBuffers(inX, inY);

  if(density == 0.f) {
    dotsX[0] = inX;
    dotsY[0] = inY;
    dotsW[0] = brightness;
    return 1;
  }

  switch(drawMode)
  {
  case LINEAR:        return dotsLinear();
  case CUBIC_HERMITE: return dotsCubicHermite();
  case QUADRATIC:     return dotsQuadratic();
  default:            return 0;
  }
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::updatePointBuffers(TCor newX, TCor newY)
{
  rsArrayTools::pushFrontPopBack4(newX, x);
  rsArrayTools::pushFrontPopBack4(newY, y);
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::shiftPointBuffers(TCor shiftX, TCor shiftY)
{
  for(int i = 0; i < 3; i++) {
    x[i] += shiftX;
    y[i] += shiftY;
  }
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::reset(TCor x_, TCor y_)
{
  rsArrayTools::fillWithValue(x, 4, x_);
  rsArrayTools::fillWithValue(y, 4, y_);
  wOld = TWgt(0);
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::numDotsForSegment(TCor segmentLength)
{
  TCor minDotDistance = 1; // maybe make user-adjustable
  int numDots = rsMax(1, (int)floor(density*segmentLength/minDotDistance));
  if(maxNumDots > -1)
    numDots = rsMin(numDots, maxNumDots);
  return rsMin(numDots, dotBufferLength);
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::getStartAndEndWeights(int numDots, TWgt* wStart, TWgt* wEnd)
{
  TWgt scaler = TWgt(1);
  if(scaleWeightsByNumDots)
    scaler /= (TWgt)numDots;

  if(useGradient) {
    TWgt wt = brightness * scaler;   // target weight that would be used if we don't do gradients
    *wEnd = (2.f*wt) - wOld;         // desired endpoint color
    *wEnd = rsMax(*wEnd, wt);        // wEnd could come out negative, use wt as lower bound
    *wStart = wOld;
  }
  else {
    *wStart = brightness * scaler;
    *wEnd   = brightness * scaler;
  }
  wOld = *wEnd;
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::dotsLinear()
{
  // distance computation:
  TCor dx = x[1]-x[2];
  TCor dy = y[1]-y[2];
  TCor pixelDistance = sqrt(dx*dx + dy*dy);
  int numDots = numDotsForSegment(pixelDistance);

  // start/end weight computation:
  TWgt w1, w2, dw;
  getStartAndEndWeights(numDots, &w1, &w2);
  dw = w2-w1;

  // dot computation:
  TCor scaler = (TCor)(1.0 / numDots);
  TCor k;
  for(int i = 0; i < numDots; i++) {
    k = scaler * (i+1);     // == (i+1) / numDots
    dotsX[i] = x[2] + k*dx;
    dotsY[i] = y[2] + k*dy;
    dotsW[i] = w1 + TWgt(k)*dw;
  }

  return numDots;
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::dotsCubicHermite()
{
  // compute polynomial coeffs:
  static const TCor h = TCor(0.5);
  cubicSplineArcCoeffs2D(
    x[2], h*(x[1]-x[3]), y[2], h*(y[1]-y[3]),  // older inner point comes first
    x[1], h*(x[0]-x[2]), y[1], h*(y[0]-y[2]),  // newer inner point comes second
    a, b);

  // factor out (duplicated in dotsQuadratic):
  int numDots = prepareParameterArray();
  TWgt w1, w2;
  getStartAndEndWeights(numDots, &w1, &w2);  // start/end weight computation
  dotsCubic(w1, w2, numDots);                // dot computation
  return numDots;
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::dotsQuadratic()
{
  static const TCor h = TCor(0.5);
  quadraticOrCubicSplineArcCoeffs2D(
    x[2], h*(x[1]-x[3]), y[2], h*(y[1]-y[3]), 
    x[1], h*(x[0]-x[2]), y[1], h*(y[0]-y[2]), 
    a, b);

  // factor out:
  int numDots = prepareParameterArray();
  TWgt w1, w2;
  getStartAndEndWeights(numDots, &w1, &w2);
  dotsCubic(w1, w2, numDots);
  return numDots;
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::prepareParameterArray()
{
  // distance computation:
  TCor splineLength;
  int numDots;
  if(normalizeDensity) {
    // obtain arc-length s(t) as function of parameter t at N sample points:
    static const int N = 17; // make this user-adjustable (setDensityCompensationPrecision)
    r.resize(N);  // move into setDensityCompensationPrecision
    s.resize(N);
    rsArrayTools::fillWithRangeLinear(&r[0], N, TCor(0), TCor(1));
    cubicSplineArcLength2D(a, b, &r[0], &s[0], N);
    splineLength = s[N-1]; // last value in s is total length: s(t=1)
    numDots = numDotsForSegment(splineLength);

#ifdef RS_DEBUG_PLOTTING
    GNUPlotter::plot(N, &r[0], &s[0]);
#endif

    // compute t-array by evaluating t(s), i.e. the function that is inverse to s(t):
    u.resize(numDots); // equally spaced values of s at which we compute t(s)
    t.resize(numDots); // t(s), i.e t(u[i])
    TCor scaler = splineLength / TCor(numDots);
    for(int i = 0; i < numDots; i++)
      u[i] = scaler * i;
    resampleNonUniformLinear(&s[0], &r[0], N, &u[0], &t[0], numDots);

#ifdef RS_DEBUG_PLOTTING
    GNUPlotter::plot(numDots, &u[0], &t[0]);
#endif
  }
  else {
    // crude estimate of arc-length as distance between endpoints:
    TCor dx = x[1]-x[2]; 
    TCor dy = y[1]-y[2]; 
    splineLength = sqrt(dx*dx + dy*dy); 
    numDots = numDotsForSegment(splineLength);

    // compute t-array:
    t.resize(numDots);
    TCor scaler = (TCor)(1.0 / numDots);
    for(int i = 0; i < numDots; i++)
      t[i] = scaler * i;
  }
  return numDots;
}

template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::dotsCubic(TWgt w1, TWgt w2, int numDots)
{
  TWgt dw = w2-w1;                      // weight difference
  TCor scaler = (TCor)(1.0 / numDots);  // not 1/(numDots-1) because last dot of this call is drawn 
                                        // as first dot in next call? ..avoids drawing it twice?
  for(int i = 0; i < numDots; i++) {
    dotsX[i] = rsPolynomial<TCor>::evaluate(t[i], a, 3);
    dotsY[i] = rsPolynomial<TCor>::evaluate(t[i], b, 3);
    dotsW[i] = w1 + TWgt(t[i])*dw;

    //rsAssert(rsIsFiniteNumber(dotsX[i]));
    //rsAssert(rsIsFiniteNumber(dotsY[i]));
    //rsAssert(rsIsFiniteNumber(dotsW[i]));
  }
}



// Ideas:
// -maybe let the client retrieve the polynomial coefficients, so it can produce the dot-buffers 
//  itself if it wants to
// -instead of the polynomial coefficients, the client shpuld also be able to retrieve a 
//  parametrization in terms of control points as in (cubic) Bezier curves - we need functions like
//  hermiteToBezier/bezierToHermite (converting control points from/to avlue and deribative) and/or 
//  polyToBezier/bezierToPoly (convert cotrol points from/to polynomial coeffs). these conversion 
//  functions may be put as static functions into a class rsCubicSpline...perhaps templatized on 
//  the vector type (may be 2D or 3D)
// -maybe this class could also be refactored to work with vectors - the a,b coeff arrays could be 
//  replaced with a single, vector-valued array of coeffs - also the x,y buffers could be replaced
//  by a single vector valued buffer
// -try to draw the spline with short line-segments
//  -maybe in this case, the higher density in regions of high curvature is actually is actually 
//   beneficial, so no density compensation is required?
// -try quartic interpolation spline 
//  -maybe less overshoot?
// -try to fix the 2nd derivative at the joints to 0
//  -needs a quinitc spline
// -or try to estimate meaningful values for the 2nd derivatives from the data
// -try a simpler (to compute) density compensation:
//  -evaluate instantaneous speed sqrt( (dx/dt)^2 + (dy/dt)^2 ) at each spline point and skip dots or
//   insert extra dots depending on the value
//  -maybe use |dx/dt| + |dy/dt| instead of the sqrt
//  -maybe instead of skip/insert dots, change the brightness of the dot to be drawn
//  -maybe instead of keeping track of the actual distance (requiring a square-root for each 
//   dot), we could use something based on the 2nd derivative (which is easy to compute), maybe
//   brightness *= 1 / (1+cx*cx+cy*cy), cx = d2x/dt2, cy = d2y/dt2 - 2nd derivatives of x,y 
//   with respect to t

// try a quadratic spline (...done in dotsQuadratic):
// x(t) = a0 + a1*t + a2*t^2, x'(t) = a1 + 2*a2*t
// y(t) = b0 + b1*t + b2*t^2, y'(t) = b1 + 2*b2*t
// and require only dy/dx = x'/y' to be equal at the joints. This is less restrictive than 
// requiring dx/dt and dy/dt to be simultaneously equal. Does that also ensure parametric 
// continuity or only geometric continuity (i think, the former)? The constraint equations are:
// x(t=0) = x0 = a0
// y(t=0) = y0 = b0
// x(t=1) = x1 = a0 + a1 + a2
// y(t=1) = y1 = b0 + b1 + b2
// s(t=0) = s0 = dy/dx|t=0 -> s0 = (b1+2*b2*0)/(a1+2*a2*0) = b1/a1 -> s0*a1 = b1
// s(t=1) = s1 = dy/dx|t=1 -> s1 = (b1+2*b2*1)/(a1+2*a2*1) -> s1*(a1+2*a2) = b1+2*b2
// where dy/dx = (dy/dt)/(dx/dt). So, we have 6 equations for 6 unknowns, so this should work - but 
// what if dx/dt = 0? We could use r0 = 1/s0 = dx/dy in this case, but what if both dx/dt and dy/dt
// are both zero? ...maybe treat such special cases separately and set dy/dx = 1
// Maybe this quadratic interpolation could also be improved by normalizing the integral under
// x(t) and y(t) to be equal to the integral of a linear interpolant. This would give two 
// additional constraint equations calling for cubics. But these would be different cubics than 
// the Hermite cubics that we have now....perhaps better? less wiggle/overshoot? -> try it!
// Also, we could use as additional constraints that the 2nd derivatives d2y/dx2|t=0, d2y/dx2|t=1 
// should be equal at the joints - 2 additional constraints -> two additional parameters (a3,b3)
// -> 8 coupled equations
// solving the system with sage:
/*
# variables and function definitions:
var("a0 a1 a2 b0 b1 b2 t x0 x1 y0 y1 s0 s1 r0 r1")
x(t)  = a0 + a1*t + a2*t^2
y(t)  = b0 + b1*t + b2*t^2
dx(t) = diff(x(t), t)         # dx/dt
dy(t) = diff(y(t), t)         # dy/dt

# constraint equations and solution:
e1 = x(0) == x0
e2 = y(0) == y0
e3 = x(1) == x1
e4 = y(1) == y1
e5 = dy(0)/dx(0) == s0
e6 = dy(1)/dx(1) == s1
solve([e1,e2,e3,e4,e5,e6], a0,a1,a2,b0,b1,b2)

gives:
a0 == x0 
a1 == 2*(s1*x0 - s1*x1 - y0 + y1)/(s0 - s1)
a2 == -((s0 + s1)*x0 - (s0 + s1)*x1 - 2*y0 + 2*y1)/(s0 - s1)
b0 == y0
b1 == 2*(s0*s1*x0 - s0*s1*x1 - s0*y0 + s0*y1)/(s0 - s1)
b2 == -(2*s0*s1*x0 - 2*s0*s1*x1 - (s0 + s1)*y0 + (s0 + s1)*y1)/(s0 - s1)

the case s0 == s1 needs to be treated separately to avoid div-by-zero - maybe use a line in this 
case, precompute: 1/(s0-s1), (s0+s1), 

using s0, r1:
a0 == x0, 
a1 == -2*(r1*y0 - r1*y1 - x0 + x1)/(r1*s0 - 1), 
a2 == -((r1*s0 + 1)*x0 - (r1*s0 + 1)*x1 - 2*r1*y0 + 2*r1*y1)/(r1*s0 - 1), 
b0 == y0, 
b1 == -2*(r1*s0*y0 - r1*s0*y1 - s0*x0 + s0*x1)/(r1*s0 - 1), 
b2 == -(2*s0*x0 - 2*s0*x1 - (r1*s0 + 1)*y0 + (r1*s0 + 1)*y1)/(r1*s0 - 1)

using r0, r1:
a0 == x0, 
a1 == 2*(r0*r1*y0 - r0*r1*y1 - r0*x0 + r0*x1)/(r0 - r1), 
a2 == -(2*r0*r1*y0 - 2*r0*r1*y1 - (r0 + r1)*x0 + (r0 + r1)*x1)/(r0 - r1), 
b0 == y0, 
b1 == 2*(r1*y0 - r1*y1 - x0 + x1)/(r0 - r1), 
b2 == -((r0 + r1)*y0 - (r0 + r1)*y1 - 2*x0 + 2*x1)/(r0 - r1)

using r0, s1:
a0 == x0, 
a1 == -2*(r0*s1*x0 - r0*s1*x1 - r0*y0 + r0*y1)/(r0*s1 - 1), 
a2 == ((r0*s1 + 1)*x0 - (r0*s1 + 1)*x1 - 2*r0*y0 + 2*r0*y1)/(r0*s1 - 1), 
b0 == y0, 
b1 == -2*(s1*x0 - s1*x1 - y0 + y1)/(r0*s1 - 1), 
b2 == (2*s1*x0 - 2*s1*x1 - (r0*s1 + 1)*y0 + (r0*s1 + 1)*y1)/(r0*s1 - 1)
*/