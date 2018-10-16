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

// move this code to FilterPlotter (or maybe there already is something similar?) - allows user
// to create an object and access further settings
void plotFrequencyResponse(std::vector<double>& f, std::vector<double>& dB,
  std::vector<double>& degrees, bool logFreq)
{
  //rsAssert(size(f) == size(dB) && size(f) == size(degrees));
  int N = (int) f.size();
  GNUPlotter p;
  p.addDataArrays(N, &f[0], &dB[0]);
  p.addDataArrays(N, &f[0], &degrees[0]);
  p.setPixelSize(1200, 400);
  p.setTitle("Filter Frequency Response");
  //p.setGraphColors("A00000", "909000", "008000", "0000A0", "800080",
  //  "A00000", "909000", "008000", "0000A0", "800080" );
  if(logFreq)
    p.addCommand("set logscale x");
  //p.addCommand("set xrange  [0.0625:16]");
  //p.addCommand("set yrange  [-100:0]");
  //p.addCommand("set y2range [-450:0]");
  p.addCommand("set xlabel \"Frequency in Hz\"");
  p.addCommand("set ylabel \"Magnitude in dB\"");
  p.addCommand("set y2label \"Phase in Degrees\"");
  //p.addCommand("set xtics 2");    // factor 2 between (major) frequency axis tics
  //p.addCommand("unset mxtics");   // no minor tics for frequency axis
  p.addCommand("set ytics 10");   // 10 dB steps for magnitude axis
  p.addCommand("set y2tics 45");  // 45� steps for phase axis

  // add magnitude and phase graphs:
  p.addGraph("i 0 u 1:2 w lines lw 1.5 axes x1y1 notitle");
  p.addGraph("i 1 u 1:2 w lines lw 1.5 axes x1y2 notitle");
  p.plot();
}



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

void plotMagnitudeResponse(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(specBA);
  plt.addCommand("set xtics ('0' 0, 'sr/4' pi/2, 'sr/2' pi)");
  plt.plotMagnitude(1000, 0.0, PI, false, false);
}

void plotPolesAndZeros(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(specBA);
  plt.plotPolesAndZeros();
}

void showFilterPlots(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  plotMagnitudeResponse(specBA);
  plotPolesAndZeros(specBA);
}




#endif
