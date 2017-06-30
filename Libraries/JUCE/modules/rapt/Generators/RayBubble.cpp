template<class T>
rsRayBubble<T>::rsRayBubble()
{

}

// setup:


// processing:

template<class T>
void rsRayBubble<T>::getSampleFrame(T &x, T &y)
{
  // assign outputs:
  x = xc;
  y = yc;

  // Compute new (tentative) x,y coordinates:
  T xn, yn;
  xn = xc + dx;
  yn = yc + dy;

  // Reflect, if new coordinates are outside elliptic enclosure:
  while(ellipse.isPointOutside(xn, yn))  
  {
    // Compute intersection point between current line segment and ellipse:
    T ti, xi, yi;
    ellipse.lineIntersectionParameter(xc, dx, yc, dy, &xi, &ti); // xi is dummy for 2nd solution
    xi = xc + ti*dx;                                             // which is irrelevant
    yi = yc + ti*dy; 

    // for debug - check that xi,yi is indeed on the ellipse:
    //T err = ellipse.evaluate(xi, yi);   // should be 0 up to roundoff

    // Reflect new point in a line tangent to the ellipse at intersection point xi, yi:
    T A, B, C;
    ellipse.getTangentCoeffs(xi, yi, &A, &B, &C);
    rsLine2D<T>::reflectPointInLine(xn, yn, A, B, C, &xn, &yn);

    // Compute new dx,dy by taking the new direction vector (which points from the intersection 
    // point to the reflected point) and adjusting its length according to the desired speed:
    dx = xn-xi;
    dy = yn-yi;
    T scaler = speed / sqrt(dx*dx + dy*dy);
    dx *= scaler;
    dy *= scaler;
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
