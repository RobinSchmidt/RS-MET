#ifndef RS_RESPONSEGETTERS_INL
#define RS_RESPONSEGETTERS_INL

#include "ResponseGetters.h"

namespace RSLib
{
  template<class T>
  void getImpulseResponse(T &module, double *h, int N)
  {
    module.reset();
    h[0] = module.getSample(1.0);
    for(int n = 1; n < N; n++)
      h[n] = module.getSample(0.0);
  }

  template<class T>
  void getStepResponse(T &module, double *s, int N)
  {
    module.reset();
    for(int n = 0; n < N; n++)
      s[n] = module.getSample(1.0);
  }

  template<class T>
  void getResponse(T &module, double *x, double *y, int N)
  {
    module.reset();
    for(int n = 0; n < N; n++)
      y[n] = module.getSample(x[n]);
  }

}

#endif