#ifndef RAPT_PLOTTING_CPP // why? i think, it may be deleted
#define RAPT_PLOTTING_CPP

#include "Plotting.h" 

//void plotData(int N, float *x, float *y1, float *y2, float *y3, float *y4, float *y5)
//{
//  GNUPlotter plt;
//  plt.addDataArrays(N, x, y1, y2, y3, y4, y5);
//  plt.plot();
//}

template<class TSig, class TFlt>
vector<TSig> impulseResponse(TFlt &filter, int N, TSig scale)
{
  vector<TSig> y(N);
  filter.reset();
  y[0] = filter.getSample(scale);
  for(int n = 1; n < N; n++)
    y[n] = filter.getSample(0.0);
  return y;
}

template<class TSig, class TFlt>
void plotImpulseResponse(TFlt &filter, int N, TSig scale)
{
  vector<TSig> y = impulseResponse(filter, N, scale);
  GNUPlotter plt;
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}

#endif