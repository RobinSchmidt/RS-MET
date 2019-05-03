#ifndef RAPT_RESPONSEGETTERS_H
#define RAPT_RESPONSEGETTERS_H


// functions to obtain variuos responses of a system to various standard inputs such as an 
// impulse, step, etc. ...maybe this should not be part of rapt

/** Fills the array h (of length N) with the impulse-response of the module that is passed by
reference. The class to which the module belongs must provide a reset() and getSample()
function. */
template<class TFlt, class TSig>
inline void getImpulseResponse(TFlt &module, TSig *h, int N)
{
  module.reset();
  h[0] = module.getSample(1.0);
  for(int n = 1; n < N; n++)
    h[n] = module.getSample(0.0);
}

/** Fills the array s (of length N) with the unit-step response of the module that is passed by
reference. The class to which the module belongs must provide a reset() and getSample()
function. */
template<class TFlt, class TSig>
inline void getStepResponse(TFlt &module, TSig *s, int N)
{
  module.reset();
  for(int n = 0; n < N; n++)
    s[n] = module.getSample(1.0);
}

/** Applies the given module to the signal x and writes the result into y, both of which are
assumed to be of length N. */
template<class TFlt, class TSig>
inline void getResponse(TFlt &module, TSig *x, TSig *y, int N)
{
  module.reset();
  for(int n = 0; n < N; n++)
    y[n] = module.getSample(x[n]);
}


#endif