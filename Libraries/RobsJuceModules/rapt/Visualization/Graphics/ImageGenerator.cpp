
template<class TPix, class TVal> 
void rsImageGenerator<TPix, TVal>::_drawImplicitCurve(const std::function<TVal(TVal, TVal)>& f, 
  TVal c, TVal x0, TVal y0, rsImage<TPix>& img, TPix color, bool clockwise) 
{
  rsAssert(f(x0, y0) == c, "x0,y0 should solve f(x0,y0) = c" );  
  // todo: use tolerance - should maybe be some fraction of c? ..but that would mean that if c is 
  // zero, we would have zero tolerance...maybe max(c*eps, eps)?

  painter.setImageToPaintOn(&img);

  // figure out start pixel:
  TVal x = x0;
  TVal y = y0;
  TVal xMaxPixel = TVal(img.getWidth()  - 1);   // maximum x-coordinate in pixel coordinates
  TVal yMaxPixel = TVal(img.getHeight() - 1);   // same for y-coordinate
  TVal sx  = xMaxPixel   / (xMax-xMin);         // one x-pixel in world coordinates
  TVal sy  = yMaxPixel   / (yMax-yMin);
  TVal sxi = (xMax-xMin) / xMaxPixel;
  TVal syi = (yMax-yMin) / yMaxPixel;

  // Convert (x,y) to pixel coordinates and draw 1st point. In order to avoid double-drawing a point
  // in the potential recursive call (for 2-part curves), do it conditionally:
  TVal px, py;
  if(!clockwise) { // avoid drawing the starting point again in the recursive call
    px = rsLinToLin(x, xMin, xMax, TVal(0), xMaxPixel);
    py = rsLinToLin(y, yMin, yMax, yMaxPixel, TVal(0));
    painter.paintDot(px, py, color); }

  // Main loop over the pixels on the curve:
  int iterations = 0;
  while(true)
  {
    // Figure out gradient (dx,dy) and contour direction (rx,ry) which is perpendicular to the 
    // gradient (i.e. 90° rotated):
    TVal h  = 1.e-8;  // ad-hoc - make parameter, or maybe use sqrt(epsilon)
    TVal dx = (f(x+h, y) - f(x-h, y)) / (TVal(2)*h);  // x-component of gradient
    TVal dy = (f(x, y+h) - f(x, y-h)) / (TVal(2)*h);  // y-component of gradient
    TVal rx, ry;
    if(clockwise) { rx =  dy; ry = -dx; }
    else          { rx = -dy; ry =  dx; }

    // Check, if the current segment is horizontalish/flat or verticalish/steep. In the flat case, 
    // advance x by one pixel and y by a distance derived from the direction vector - in the steep 
    // case, do it the other way around:
    bool flat = rsAbs(rx*sx) > rsAbs(ry*sy);   // curve sgement is horizontalish
    if(flat) {
      dx = rsSign(rx) / sx;    // this step should translate to 1 pixel left or right -> check this!
      dy = dx * ry/rx;         // the y-step is proportional to the x-step - is this formula the best we can do?
      x += dx;                 // walk one pixel left or right
      y += dy; }
    else {              // verticalish
      dy = rsSign(ry) / sy;
      dx = dy * rx/ry;
      x += dx;
      y += dy; }

    // In the step just taken, we may have drifted off the contour line due to approximation 
    // errors, so we fix this by refining x or y such that we land on the contour again. We use 1D 
    // Netwon iteration with numeric derivatives. If our direction is horizontalish, we change y, 
    // otherwise, we change x. To catch convergence problems, we verify that the error actually 
    // went down in the iteration step - if it didn't, we restore old value from before the step 
    // and break out of the loop.
    TVal err = f(x,y) - c;
    //TVal tol = 1.e-12;
    TVal tol = 1.e-15;   // maybe use something like 10*epsilon
    //tol = T(0);  // test
    if(!flat) {               // y-step is larger (steep) -> refine x
      while(rsAbs(err) > tol)  {
        dx    = (f(x+h, y) - f(x-h, y)) / (TVal(2)*h); // central difference as approximation to the
        x     = x - err / dx;                       // partial derivative with respect to x
        TVal old = err;
        err   = f(x,y) - c;
        if(rsAbs(old) <= rsAbs(err)) {
          x += old / dx;  // restore old value bcs old error was better
          break; }}}
    else {
      while(rsAbs(err) > tol)  {
        dy = (f(x, y+h) - f(x, y-h)) / (TVal(2)*h);
        y   = y - err / dy;
        TVal old = err;
        err = f(x,y) - c;  
        if(rsAbs(old) <= rsAbs(err)) {
          y += old / dy;
          break; }}}
    // factor out the loops in these two branhces into two functions: newtonRefineX, newtonRefineY


    // The stopping criterion for closed curves is that we have come back to (or very close to) the 
    // starting point again. To close the curve, we paint one last pixel, whose brightness is 
    // scaled by how far we are away from the starting point (...this is not yet perfect - it looks 
    // like the start/end point is still drawn a bit brighter than the rest of the curve...). The
    // && iteration >= 1 is for avoiding spuriously breaking out of the loop in the very first 
    // iteration due to roundoff errors.
    if(rsAbs(x-x0) < sxi && rsAbs(y-y0) < syi && iterations >= 1) {
      dx  = (x-x0)*sx;
      dy  = (y-y0)*sy;
      err = sqrt(dx*dx + dy*dy);
      px  = rsLinToLin(x, xMin, xMax, TVal(0), xMaxPixel);
      py  = rsLinToLin(y, yMin, yMax, yMaxPixel, TVal(0));
      painter.paintDot(px, py, TPix(err)*color);          // last pixel too bright (really?)
      //painter.paintDot(px, py, TPix(sqrt(err))*color);
      //painter.paintDot(px, py, TPix(err*err)*color);        // last pixel too dark
      //painter.paintDot(px, py, TPix(pow(err, 1.25))*color);
      break; }

    // Convert point in world coordinates (x,y) to pixel coordinates (px,py) and paint it:
    px = rsLinToLin(x, xMin, xMax, TVal(0), xMaxPixel);
    py = rsLinToLin(y, yMin, yMax, yMaxPixel, TVal(0));
    painter.paintDot(px, py, color);

    // The stopping criterion for open curves is that we reach the image boundary - in such cases, 
    // we have just drawn one part of the curve (think of a hyperbola, for example), so we call 
    // ourselves recursively to draw the second part as well. To avoid infinite recursion, the 
    // recursive call is done only if clockwise == false, which should always be the case when 
    // being called from client code but is *not* the case for a recursive call (we pass true 
    // here):
    if(x < xMin || x > xMax || y < yMin || y > yMax) {
      if(clockwise == false)
        _drawImplicitCurve(f, c, x0, y0, img, color, true);
      break; }

    // Avoid infinite loops - this should not normally happen:
    iterations++;
    if(iterations > 100000) { // is that enough? may some curves have more pixels?
      rsError("drawImplicitCurve ran into infinite loop");
      break;  }
  }
}
// maybe fill the shape by coloring all pixels for which f(x,y) < c









template<class T>
T squaredDistance(T x1, T y1, T x2, T y2)
{
  T dx = x2 - x1;
  T dy = y2 - y1;
  return dx*dx + dy*dy;
}

template<class T>
T distance(T x1, T y1, T x2, T y2)
{
  return sqrt(squaredDistance(x1, y1, x2, y2));
}
// move to somewhere else

/** Computes the minimum value of the squared distances from (x0,y0) to the points in the x,y 
arrays. */
template<class T>
T minSquaredDistance(T x0, T y0, T* x, T* y, int N)
{
  T d2min = RS_INF(T);
  for(int i = 0; i < N; i++) {
    T d2 = squaredDistance(x0, y0, x[i], y[i]);
    if(d2 < d2min)
      d2min = d2; }
  return d2min;
}

template<class T>
T minDistance(T x0, T y0, T* x, T* y, int N)
{
  return sqrt(minSquaredDistance(x0, y0, x, y, N));
}



template<class TPix, class TVal>
void rsImageGenerator<TPix, TVal>::distanceMap(rsImage<TPix>& img, TVal* x, TVal* y, int N)
{
  for(int j = 0; j < img.getHeight(); j++) {
    for(int i = 0; i < img.getWidth(); i++) {
      TVal xp = TVal(i); 
      TVal yp = TVal(j);
      TVal d  = minDistance(xp, yp, x, y, N);
      img(i, j) = TPix(d); }}
}


template<class TPix, class TVal> 
TVal rsImageGenerator<TPix, TVal>::spiralRidge(TVal x, TVal y, TVal a, TVal p, TVal sign, 
  int profile, TVal exponent)
{
  // sanity check inputs:
  rsAssert(sign == +1 || sign == -1,     "sign must be +-1");
  rsAssert(profile >= 0 && profile <= 3, "invalid profile");

  // compute raw height:
  TVal r = sqrt(x*x + y*y);
  if(r == 0.0) return 0.0;           // avoid log-of-zero
  TVal t  = log(r) / a;              // parameter t for point on the spiral with radius r
  TVal xs = r * cos(sign * t + p);   // x on the spiral for the given t
  TVal ys = r * sin(sign * t + p);   // y on the spiral for the given t
  TVal d  = distance(xs, ys, x, y);  // distance of input point to point on the spiral
  TVal h  = pow(0.5*d/r, exponent);  // height

  // apply shaping of the height profile:
  if(profile == 2) return h;                        // 2: rectified sine (comes out by raw formula)
  if(profile == 3) return 1-h;                      // 3: inverted rectified sine
  h = asin(h) / (0.5*PI);                           // convert to triangular
  if(profile == 0) return h;                        // 0: triangular
  if(profile == 1) return 0.5*(sin(PI*(h-0.5))+1);  // 1: sinusoidal
  return 0;                                         // unknown profile
}
// optimize: the sqrt in the distance computation can be avoided: compute the distance-squared and 
// then use 0.5*exponent in the subsequent pow call
//
// The algo computes the distance of (x,y) to a point on the spiral that has the same radius as 
// (x,y). It happens that the height-profile (as function of radius for a given angle) comes out as
// a rectified sine shape (when the radius is used a x-axis and the x-axis is logarithmically 
// scaled)
// what if sign is not +-1? what if we use different factors for x- and y: 
//   xs = r * cos(wx * t + px); ys = r * sin(wy * t + py);
// ..in this case, the profile computations will very likely become invalid because the raw profile 
// is not a rectified sine anymore - yes - using, for example cos(2*..),sin(3*..) gives a nice 
// effect - one should increase the shrink-factor accordingly because it get denser otherwise
// for the original profile and the profile converted to a full sine, it makes visually not 
// qualitative difference, when we invert all color channels - for the triangular profile, it does.
// maybe to create audio-signals, we could use dx = (xs-x)/r; dy = (ys-y)/r; as left and right 
// channel signal - but what should the input be? we don't have x,y pixel coordinates as inputs but
// time instants
// 
// when we use d / r, the birghtness of the white ridges is independent for the distance to the 
// center - using a power with exponent < 1, we get a darkening effect towrd the center - but mybe 
// such an effect can be applied as post-processing: 
// circularDarkening(img, x, y, amount)
//   img(i,j) /= pow(r, amount)
//
// try tL = atan2(y,x) + 2*k*pi where k = floor(t0/(2*pi)), tR = tL + 2*pi, compute (xL,yL),(xR,yR)
// by the parametric spiral equations, compute distances dL,dR and use minimum
//
// -these are not the actual distances to the nearest points on the spiral but rather the distances 
//  to two concentric circles that approximate the spiral at the given angle - but they can be used 
//  as an initial estimate for computing the actual distance via netwon iteration - maybe this 
//  refinement can be made optional, controlled by a boolean parameter
//
// or maybe use a simpler linear spiral:
//   f(t) = t*cos(t), g(t) = t*sin(t)
// see:
// https://en.wikipedia.org/wiki/Logarithmic_spiral
// https://en.wikipedia.org/wiki/Archimedean_spiral