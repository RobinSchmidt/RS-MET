template<class T>
rsRayBubble<T>::rsRayBubble()
{

}

// setup:


// processing:

template<class T>
void rsRayBubble<T>::getLineEllipseIntersectionPoint(T* xi, T* yi)
{
  // we use xi,yi temporarily for the two solution for the intersection parameter ti, yi is the
  // relevant solution (the greater value)
  ellipse.lineIntersectionParameter(xc, dx, yc, dy, xi, yi); // xi <- ti1, yi <- ti2
  *xi = xc + *yi * dx;  // yi is still ti2
  *yi = yc + *yi * dy;  // ...now not anymore
}

template<class T>
void rsRayBubble<T>::reflectInTangentAt(T xt, T yt, T* x, T *y)
{
  T A, B, C;
  ellipse.getTangentCoeffs(xt, yt, &A, &B, &C);
  rsLine2D<T>::reflectPointInLine(*x, *y, A, B, C, x, y);
}

template<class T>
void rsRayBubble<T>::updateDirectionVector(T xi, T yi, T xn, T yn)
{
  dx = xn-xi;
  dy = yn-yi;
  T scaler = speed / sqrt(dx*dx + dy*dy);
  dx *= scaler;
  dy *= scaler;
}

template<class T>
void rsRayBubble<T>::getSampleFrame(T &x, T &y)
{
  // assign outputs:
  x = xc;
  y = yc;

  // Compute new (tentative) x,y coordinates:
  T xn, yn;     // do we actually need these temporaries?
  xn = xc + dx; // can't we use xc, yc directly?
  yn = yc + dy;

  // Reflect, if new coordinates are outside elliptic enclosure:
  while(ellipse.isPointOutside(xn, yn))  
  {
    T xi, yi;
    getLineEllipseIntersectionPoint(&xi, &yi); // intersection between line segment and ellipse
    //T err = ellipse.evaluate(xi, yi);        // for debug - should be 0 up to roundoff
    reflectInTangentAt(xi, yi, &xn, &yn);      // reflect new point in tangent at intersection
    updateDirectionVector(xi, yi, xn, yn);     // new direction points from intersection to 
                                               // reflected point
  }

  // update current coordinates:
  xc = xn;
  yc = yn;
}

template<class T>
void rsRayBubble<T>::reset()
{
  xc = x0; 
  yc = y0;
  dx = speed * cos(angle);
  dy = speed * sin(angle);
}


