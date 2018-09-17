template<class TCor, class TWgt>
rsRealTimeSpline<TCor, TWgt>::rsRealTimeSpline()
{
  reset();
}

template<class TCor, class TWgt>
int rsRealTimeSpline<TCor, TWgt>::getDotsForInputPoint(
  TCor newX, TCor newY, TCor* dotsX, TCor* dotsY, TWgt weights, int xywLength)
{

}


template<class TCor, class TWgt>
void rsRealTimeSpline<TCor, TWgt>::reset(TCor x_, TCor y_)
{
  rsArray::fillWithValue(x, 4, x_);
  rsArray::fillWithValue(y, 4, y_);
  cOld = TWgt(0);
}


// -functionality for thsi class is currently scattered over rsImagePainter and rsPhaseScopeBuffer

// Ideas:
// -make density normalization optional (it's expensive)
// -try to draw the spline with short line-segments
//  -maybe in this case, the higher density in regions of high curvature is actually is actually 
//   beneficial, so no density compensation is required?
// -try quartic interpolation spline 
//  -maybe less overshoot?
// -try to fix the 3rd deruvative at the joints to 0
//  -needs a quinitc spline
// -try a simpler (to compute) density compensation:
//  -evaluate instantaneous speed sqrt( (dx/dt)^2 + (dy/dt)^2 ) at each spline point and skip dots or
//   insert extra dots depending on the value
//  -maybe use |dx/dt| + |dy/dt| instead of the sqrt
//  -maybe instead of skip/insert dots, change the brightness of the dot to be drawn

// try a quadratic spline:
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
