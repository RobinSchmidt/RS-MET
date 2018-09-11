template<class T>
void rsNumericDerivative(T *x, T *y, T *yd, int N, bool extrapolateEnds)
{
  T dxl, dxr, dx, a, b; 

  for(int n = 1; n < N-1; n++)
  {
    dxl   = x[n] - x[n-1];
    dxr   = x[n+1] - x[n];
    dx    = dxl + dxr;
    yd[n] = dxr * (y[n]-y[n-1])/(dxl*dx) + dxl * (y[n+1]-y[n])/(dxr*dx);
  }

  // todo: save the left weight and use it as right weight in the next iteration
  // -make it possible to use the function in-place, i.e. y and yd point to the same memory
  //  ->avoid using y[n-1] on the right hand side

  if( extrapolateEnds == true )
  {
    a = (yd[2] - yd[1]) / (x[2] - x[1]);
    b = yd[1] - a * x[1];
    yd[0] = a*x[0] + b;
    a = (yd[N-2] - yd[N-3]) / (x[N-2] - x[N-3]);
    b = yd[N-3] - a * x[N-3];
    yd[N-1] = a*x[N-1] + b;
      // maybe this can be simplified by using rsInterpolateLinear
  }
  else
  {
    yd[0]   = (y[1]   - y[0])   / (x[1]   - x[0]);
    yd[N-1] = (y[N-1] - y[N-2]) / (x[N-1] - x[N-2]);
  }
}

template<class T>
void rsNumericIntegral(T *x, T *y, T *yi, int N)
{
  T xn, yn, zn, tmp;
  xo = x[0]; yo = y[0]; zo = T(0); yi[0] = zo; // "old" values (at index n-1)
  for(int n = 1; n < N; n++) {
    tmp = zo + (x[n]-xo)*(y[n]+yo)*T(0.5);     // compute integral by trapezoidal rule
    xo = x[n]; yo = y[n]; zo = tmp;            // update state variables
    yi[n] = tmp;                               // write integral to output array
  }
}
// it should be possible to use this function in place - yi can be the same array as y or x
// make a simplified version that doesn't need an x-array (assume distance 1 between x-avlues)
