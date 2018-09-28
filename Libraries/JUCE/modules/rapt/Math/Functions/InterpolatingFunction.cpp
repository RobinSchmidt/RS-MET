template<class Tx, class Ty>
void rsInterpolatingFunction<Tx, Ty>::interpolate(Tx *x, Ty *y, int N, Tx *xi, Ty *yi, int Ni)
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
  case LINEAR: rsInterpolateSpline(x, t, N, xi, yi, Ni, 0); break;
  case CUBIC:  rsInterpolateSpline(x, t, N, xi, yi, Ni, 1); break;
  }

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