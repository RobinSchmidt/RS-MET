template<class TSig, class TPar>
PhasorFilter<TSig, TPar>::PhasorFilter()
{
  mapper     = nullptr;
  sampleRate = 44100;
  frequency  = 1000;
  decay      = TPar(0.01);
  updateCoefficients();
  reset();
}

// setup:

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateCoefficients();
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setFrequency(TPar newFrequency)
{
  frequency = newFrequency;
  updateCoefficients();
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setDecayTime(TPar newDecay)
{
  decay = newDecay;
  updateCoefficients();
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setStateMapper(Mapper2D<TSig> *newMapper)
{
  mapper = newMapper;
}

// audio processing:

template<class TSig, class TPar>
inline void PhasorFilter<TSig, TPar>::processFrame(TSig *x, TSig *y)
{
  *x  += Axx * xOld  +  Axy * yOld;
  *y  += Ayx * xOld  +  Ayy * yOld;

  if(mapper != nullptr)
    mapper->map(x, y);

  xOld = *x;
  yOld = *y;

  // the output loudness depends of the cutoff frequency - we must scale either the input or the 
  // output in dependence of the cutoff (check, which scaling position is better behaved when we 
  // have nonlinearity - probably input?)
}

template<class TSig, class TPar>
inline TSig PhasorFilter<TSig, TPar>::getSample(TSig in)
{
  TSig x = in;
  //TSig x = 2 * (frequency/sampleRate) * in; // test
  TSig y = 0;
  processFrame(&x, &y);
  //return x;    // maybe use weighted sum of in, x and y according to output coeffs
  //return y;
  return x + y;  // test
}

// misc:

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::reset()
{
  xOld = yOld = 0;
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::updateCoefficients()
{
  TPar w = 2 * TPar(PI) * frequency / sampleRate;  // pole angle
  TPar r = exp(-1 / (decay*sampleRate));           // pole radius
  TPar c = cos(w);
  TPar s = sin(w);

  // compute rotation matrix with decay:
  // A = | Axx  Axy | =  r * | cos(w)  -sin(w) |
  //     | Ayx  Ayy |        | sin(w)   cos(w) |
  Axx = Ayy = r*c; 
  Ayx = r*s;
  Axy = -Ayx;

  // This computation assumes a pair of complex conjugate poles. Maybe we can come up with a 
  // different formula for the matrix coefficients that realizes a pair of real poles such that 
  // this filter can emulate a general biquad (transfer function wise).

  // todo: maybe keep the pole radius as member and let the user set it up directly ...maybe we
  // should factor out a baseclass PhasorFilterBase that only has pole radius r and normalized 
  // radian frequency w as parameters. we could then also use r > 1 (of course only, if some
  // saturating nonlinearity is present, otherwise it will explode)
}

//=================================================================================================

template<class T>
PhasorStateMapper<T>::PhasorStateMapper()
{
  same    = 0;   // coeff for same^2
  other   = 0;   // coeff for other^2
  cross   = 0;   // coeff for cross-term same*other
  offset  = 0;   // added constant
  satPre  = 0;   // higher values introduce a pre-gap in resonance
  satPost = 0;   // higher values shorten the resonance
}

template<class T>
void PhasorStateMapper<T>::setSameSquare(T newCoeff)
{
  same = newCoeff;
}

template<class T>
void PhasorStateMapper<T>::setOtherSquare(T newCoeff)
{
  other = newCoeff;
}

template<class T>
void PhasorStateMapper<T>::setCrossProduct(T newCoeff)
{
  cross = newCoeff;
}

template<class T>
void PhasorStateMapper<T>::setOffset(T newCoeff)
{
  offset = newCoeff;
}

template<class T>
void PhasorStateMapper<T>::setPreNormalizeSaturation(T c)
{
  satPre = T(0.25)*c*c;
}

template<class T>
void PhasorStateMapper<T>::setPostNormalizeSaturation(T c)
{
  satPost = T(0.25)*c*c;
}

template<class T>
void PhasorStateMapper<T>::map(T *xInOut, T *yInOut)
{
  T x = *xInOut;
  T y = *yInOut;

  // We make a nonlinear transformation of the xy-vector subject to the constraint that the length
  // sqrt(x^2 + y^2) must remain the same.

  // compute some nonlinear combinations of the coordinates:
  T xx = x*x;
  T xy = x*y;
  T yy = y*y;

  // compute squared length of (x,y) input vector:
  T r2 = xx + yy; 

  // apply a nonlinear transformation to the vector:
  x += xx*same  + yy*other + xy*cross + offset;
  y += xx*other + yy*same  + xy*cross + offset;

  // pre-renormalize saturation:
  x /= 1 + satPre*x*x;
  y /= 1 + satPre*y*y;

  // restore old vector length (renormalization):
  T R2 = x*x + y*y;         // new length
  if(R2 > 0)                // avoid div-by-zero
  {
    T s = sqrt(r2/R2);   // scaler to restore length
    x *= s;
    y *= s;
  }
  // maybe we should have different renormalization modes: never, always, if R2 > r2, etc.

  // post-renormalize saturation:
  T s = 1 / (1 + satPost * (x*x + y*y));
  x *= s;
  y *= s;

  // assign outputs:
  *xInOut = x;
  *yInOut = y;

  // maybe include also an input saturation before adding the prod

  // Maybe it's possible to use some function that doesn't require explicit renormalization. If 
  // we use xNew = x + f(x,y), yNew = y + g(x,y) for the application of the nonlinearity, I think
  // we should satisfy the constraint: (x + f(x,y))^2 + (y + g(x,y))^2 = x^2 + y^2, or shorter:
  // (x+f)^2 + (y+g)^2 = x^2 + y^2 - maybe we can find such a pair of functions that satisfies this
  // equation (for all values of x,y)
}