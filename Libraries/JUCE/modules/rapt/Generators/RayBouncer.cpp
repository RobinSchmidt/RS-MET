template<class T>
rsRayBouncer<T>::rsRayBouncer()
{
  reset();
}

// processing:

template<class T>
void rsRayBouncer<T>::getLineEllipseIntersectionPoint(T* xi, T* yi)
{
  // we use xi,yi temporarily for the two solutions for the intersection parameter ti, yi is the
  // relevant solution (the greater value):
  ellipse.lineIntersectionParameter(x, dx, y, dy, xi, yi); // xi <- ti1, yi <- ti2
  *xi = x + *yi * dx;  // yi is still ti2
  *yi = y + *yi * dy;  // ...now not anymore
}

template<class T>
void rsRayBouncer<T>::reflectInTangentAt(T xt, T yt, T* x, T *y)
{
  T A, B, C;
  ellipse.getTangentCoeffs(xt, yt, &A, &B, &C);
  rsLine2D<T>::reflectPointInLine(*x, *y, A, B, C, x, y);
}

template<class T>
void rsRayBouncer<T>::updateVelocity(T xi, T yi)
{
  dx = x-xi;
  dy = y-yi;
  T scaler = speed / sqrt(dx*dx + dy*dy);
  dx *= scaler;
  dy *= scaler;
}

template<class T>
void rsRayBouncer<T>::getSampleFrame(T &xOut, T &yOut)
{
  // assign outputs:
  xOut = x;
  yOut = y;

  // Compute new (tentative) x,y coordinates:
  x += dx;
  y += dy;

  // Reflect, if new coordinates are outside elliptic enclosure:
  //T tol = T(1.e-8); // to avoid div-by-almost-zero in velocity update
  T tol = 0; // for debug - maybe that's not such a good idea with the tolerance
  while(ellipse.isPointOutside(x-tol, y-tol))
  {
    T xi, yi;
    getLineEllipseIntersectionPoint(&xi, &yi); // intersection between line segment and ellipse
    //T err = ellipse.evaluate(xi, yi);        // for debug - should be 0 up to roundoff
    reflectInTangentAt(xi, yi, &x, &y);        // reflect new point in tangent at intersection
    updateVelocity(xi, yi);                    // points from intersection to reflected point
  }
  // in certain conditions we get hung up in this loop - may use an "if" or a maximum number
  // of iterations..
}

template<class T>
void rsRayBouncer<T>::reset()
{
  x  = x0; 
  y  = y0;
  dx = speed * cos(angle);
  dy = speed * sin(angle);
}

//-------------------------------------------------------------------------------------------------

template<class T>
void rsRayBouncerDriver<T>::setFrequencyAndSampleRate(T newFreq, T newRate)
{
  T k = T(10); // proportionality constant - preliminary figure out...
  T speed = k * newFreq / newRate;
  rayBouncer.setSpeed(speed);
}

//template<class T>
//void rsRayBouncerDriver<T>::setEllipseSize(T newSize)
//{
//
//}
//
//template<class T>
//void rsRayBouncerDriver<T>::setEllipseAspectRatio(T newRatio)
//{
//
//}
//
//template<class T>
//void rsRayBouncerDriver<T>::setEllipseAngle(T newAngle)
//{
//
//}

template<class T>
void rsRayBouncerDriver<T>::reset()
{
  rayBouncer.reset();
  // todo: reset modulators, when we have some...
}