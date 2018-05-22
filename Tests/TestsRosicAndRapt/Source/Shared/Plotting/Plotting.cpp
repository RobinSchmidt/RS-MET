#ifndef RAPT_PLOTTING_CPP // why? i think, it may be deleted
#define RAPT_PLOTTING_CPP

#include "Plotting.h" 

//void plotData(int N, float *x, float *y1, float *y2, float *y3, float *y4, float *y5)
//{
//  GNUPlotter plt;
//  plt.addDataArrays(N, x, y1, y2, y3, y4, y5);
//  plt.plot();
//}

/* // commented because inlined in .h now
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
*/


// new from RSLib test - needs testing:

void plotData(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

void plotData(int N, double x0, double dx, double *y1, double *y2, double *y3, double *y4, 
  double *y5)
{
  double *x = new double[N];
  for(int n = 0; n < N; n++)
    x[n] = x0 + n*dx;
  GNUPlotter plt;
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
  delete[] x;
}

void plotDataLogX(int N, double *x, double *y1, double *y2, double *y3, double *y4, double *y5)
{
  GNUPlotter plt;
  plt.setLogScale("x");
  plt.plotFunctionTables(N, x, y1, y2, y3, y4, y5);
}

void plotVector(std::vector<double> v)
{
  int N = (int) v.size();
  double *y = new double[N];
  for(int n = 0; n < N; n++)
    y[n] = v[n];
  GNUPlotter plt;
  plt.plotArrays(N, y);
  delete[] y;
}

// use shorter variable names here...
void plotSpectrogram(int numFrames, int numBins, double **s, double fs, int H)
{
  // fs: sample rate, H: hop size

  GNUPlotter p;

  // create time- and frequency axis and add the data:
  double tMax = H*(numFrames-1)/fs;
  double fMax =  0.5*fs*(numBins-1)/numBins;
  double *t = new double[numFrames];
  double *f = new double[numBins];
  //p.rangeLinear(t, numFrames, 0.0, H*(numFrames-1)/fs );
  //p.rangeLinear(f, numBins,   0.0, 0.5*fs*(numBins-1)/numBins); 
  p.rangeLinear(t, numFrames, 0.0, tMax);
  p.rangeLinear(f, numBins,   0.0, fMax); 


  p.setDataPrecision(4);                                  // make temp files smaller
  p.addDataMatrix(numFrames, numBins, t, f, s);

  // set up style settings:
  //p.setAxisLabels("Time in s", "Frequency in Hz", "Level in dB");
  p.setPixelSize(800, 400);
  p.addGraph("i 0 nonuniform matrix w image notitle");   

  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  p.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale


                                                        // maybe use signal-dependent values later:
  double dbMin = -100.0;
  double dbMax =   10.0;
  p.setRange(0.0, tMax, 0.0, fMax, dbMin, dbMax);
  // doesn't seem to have any effect on the color axis (z). it seems like GNUPlot uses 
  // autoscaling for it, no matter what.

  p.plot();
  delete[] t;
  delete[] f;
}


#endif