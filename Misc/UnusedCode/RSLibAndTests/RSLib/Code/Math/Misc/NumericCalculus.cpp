using namespace RSLib;


void RSLib::rsNumericDerivative(double *x, double *y, double *yd, int N, bool extrapolateEnds)
{
  double dxl, dxr, dx, a, b; 

  for(int n = 1; n < N-1; n++)
  {
    dxl   = x[n] - x[n-1];
    dxr   = x[n+1] - x[n];
    dx    = dxl + dxr;
    yd[n] = dxr * (y[n]-y[n-1])/(dxl*dx) + dxl * (y[n+1]-y[n])/(dxr*dx);
  }

  // todo: save the left weight and use it as right weight in the next iteration

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


