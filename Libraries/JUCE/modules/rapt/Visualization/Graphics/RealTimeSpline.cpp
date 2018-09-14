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
}


// -functionality for thsi class is currently scattered over rsImagePainter and rsPhaseScopeBuffer

// Ideas:
// -make density normalization optional (it's expensive)
// -require only dy/dx to be equal at the joints
//  -this is less restrictive than requiring dx/dt and dy/dt to be simultaneously equal 
//  -maybe then a quadratic spline is enough? 
//  -maybe that ensures only geometric but not parameteric continuity?
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
