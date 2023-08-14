template<class Tx, class Ty>
void rsInterpolatingFunction<Tx, Ty>::interpolate(
  const Tx *x, Ty *y, int N,   const Tx *xi, Ty *yi, int Ni)
{
  // apply pre-mapping to y (use a temporary buffer t):
  Ty* t = y;
  if(preMap != nullptr ) {
    t = new Ty[N];  // todo: do this only, if preMap is not the identity function
    for(int n = 0; n < N; n++)
      t[n] = preMap(y[n]);
  }

  switch(mode)
  {
  case LINEAR:        rsInterpolateSpline( x, t, N, xi, yi, Ni, 0); break; // rename to rsInterpolateHermite
  case CUBIC_NATURAL: rsNaturalCubicSpline(x, t, N, xi, yi, Ni);    break;
  case CUBIC_HERMITE: rsInterpolateSpline( x, t, N, xi, yi, Ni, 1); break;
  default: {
    rsError("Unknown interpolation mode");
    rsInterpolateSpline(x, t, N, xi, yi, Ni, 0); // handle error gracefully in release build
  };
  }
  // rsInterpolateSpline produces a general spline (linear, cubic, quintic, ...) where the last
  // parameter controls the smoothness order.

  // apply post-mapping to yi: 
  if(postMap != nullptr) {
    for(int n = 0; n < Ni; n++)
      yi[n] = postMap(yi[n]);
  }

  if(t != y)
    delete[] t;
}

template<class Tx, class Ty>
void rsInterpolatingFunction<Tx, Ty>::updateCoeffs()
{
  //rsAssert(false); // not yet implemented
}