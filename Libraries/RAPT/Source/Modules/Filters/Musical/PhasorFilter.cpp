template<class TSig, class TPar>
PhasorFilter<TSig, TPar>::PhasorFilter()
{
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

// audio processing:

template<class TSig>
inline void applyNonlinearity(TSig &x, TSig &y)
{
  // We make a nonlinear transformation of the xy-vector subject to the constraint that the length
  // sqrt(x^2 + y^2) must remain the same.
     
  // compute some nonlinear combinations of the coordinates:
  TSig xx = x*x;
  TSig xy = x*y;
  TSig yy = y*y;

  // compute squared length of (x,y) input vector:
  TSig r2 = xx + yy; 

  // coefficients for the nonlinear transformation (make them user adjustable later):
  TSig a = TSig(-0.050); // coeff for same^2
  TSig b = TSig(-0.050); // coeff for other^2
  TSig c = TSig(-0.150); // coeff for cross-term same*other
  TSig d = TSig(-0.006); // constant

  // apply a nonlinear transformation to the vector:
  x += a*xx + b*yy + c*xy + d;
  y += b*xx + a*yy + c*xy + d;

  // restore old vector length (renormalization):
  TSig R2 = x*x + y*y;      // new length
  if(R2 > 0)                // avoid div-by-zero
  {
    TSig s = sqrt(r2/R2);   // scaler to restore length
    x *= s;
    y *= s;
  }
  // maybe we should have different renormalization modes: never, always, if R2 > r2, etc.

  // Maybe it's possible to use some function that doesn't require explicit renormalization. If 
  // we use xNew = x + f(x,y), yNew = y + g(x,y) for the application of the nonlinearity, I think
  // we should satisfy the constraint: (x + f(x,y))^2 + (y + g(x,y))^2 = x^2 + y^2, or shorter:
  // (x+f)^2 + (y+g)^2 = x^2 + y^2 - maybe we can find such a pair of functions that satisfies this
  // equation (for all values of x,y)

  // \todo: make this a member function in a subclass PhasorFilterNonlinear
}

template<class TSig, class TPar>
inline void PhasorFilter<TSig, TPar>::processFrame(TSig *x, TSig *y)
{
  *x  += Axx * xOld  +  Axy * yOld;
  *y  += Ayx * xOld  +  Ayy * yOld;

  applyNonlinearity(*x, *y);  // experimental

  xOld = *x;
  yOld = *y;
}

template<class TSig, class TPar>
inline TSig PhasorFilter<TSig, TPar>::getSample(TSig in)
{
  TSig x = in;
  TSig y = 0;
  processFrame(&x, &y);
  return x;  // maybe use weighted sum of in, x and y according to output coeffs
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
}
