#ifndef RAPT_PLOTTING_H
#define RAPT_PLOTTING_H

#include "GNUPlotter.h"

/** Returns N samples of the impulse response of the passed filter as std::vector. It is necessary 
for you to pass a scale factor of the type of the filter's output signal (for example: double), 
such that the compiler can deduce the template parameter. We also use it to scale the impulse 
response, so it actually gets some purpose besides satisfying the compiler. The filter class must 
support the functions reset() and getSample() */
template<class TSig, class TFlt>
vector<TSig> impulseResponse(TFlt &filter, int length, TSig scale);
 // maybe move to a file ResponseGetters

/** Plots N samples of the impulse response of the passed filter. */
template<class TSig, class TFlt>
void plotImpulseResponse(TFlt &filter, int length, TSig scale);


#endif