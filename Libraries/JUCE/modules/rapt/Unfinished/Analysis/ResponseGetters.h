#ifndef RS_RESPONSEGETTERS_H
#define RS_RESPONSEGETTERS_H

// todo: move into test suite

namespace RSLib
{

  // functions to obtain variuos responses of a system to various standard inputs such as an 
  // impulse, step, etc.

  /** Fills the array h (of length N) with the impulse-response of the module that is passed by 
  reference. The class to which the module belongs must provide a reset() and getSample() 
  function. */
  template<class T>
  void getImpulseResponse(T &module, double *h, int N);

  /** Fills the array s (of length N) with the unit-step response of the module that is passed by 
  reference. The class to which the module belongs must provide a reset() and getSample() 
  function. */
  template<class T>
  void getStepResponse(T &module, double *s, int N);

  /** Applies the given module to the signal x and writes the result into y, both of which are
  assumed to be of length N. */
  template<class T>
  void getResponse(T &module, double *x, double *y, int N);

}

#endif
