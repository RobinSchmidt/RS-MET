#ifndef RAPT_PLOTTING_H
#define RAPT_PLOTTING_H

//#include "GNUPlotter.h"
#include "DSPPlotters.h"
using namespace std;  // try to get rid

#include "rosic/rosic.h"



///** Plots at most five y-functions against a common x-axis. */
//void plotData(int N, float *x, float *y1, float *y2 = nullptr, float *y3 = nullptr,
//  float *y4 = nullptr, float *y5 = nullptr);

/** Returns N samples of the impulse response of the passed filter as std::vector. It is necessary
for you to pass a scale factor of the type of the filter's output signal (for example: 1.0 for
double), such that the compiler can deduce the template parameter. We also use it to scale the
input impulse to the filter, so it actually gets some purpose besides satisfying the compiler. The
filter class must support the functions reset() and getSample() */
template<class TSig, class TFlt>
inline vector<TSig> impulseResponse(TFlt &filter, int length, TSig scale)
{
  vector<TSig> y(length);
  filter.reset();
  y[0] = filter.getSample(scale);
  for(int n = 1; n < length; n++)
    y[n] = filter.getSample(0.0);
  return y;
}
// move to Utilities


/** Plots N samples of the impulse response of the passed filter. */
template<class TSig, class TFlt>
inline void plotImpulseResponse(TFlt &filter, int length, TSig scale)
{
  vector<TSig> y = impulseResponse(filter, length, scale);
  GNUPlotter plt;
  plt.addDataArrays(length, &y[0]);
  plt.plot();
}

template<class TSig, class TFlt>
inline std::vector<std::complex<TSig>> getFrequencyResponse(
  TFlt &filter, const std::vector<TSig>& w)
{
  size_t N = w.size();
  std::complex<TSig> j(0,1);        // imaginary unit
  std::vector<complex<TSig>> H(N);  // H(e^jw)
  for(size_t k = 0; k < N; k++)
    H[k] = filter.getTransferFunctionAt(exp(j*w[k]));
  return H;
}
// maybe move to RAPT...but maybe use plain arrays instead of vectors there, keep convenience
// function here

/** Plots the given magnitude response in dB and phase response in degrees against the frequency
axis f. */
void plotFrequencyResponse(std::vector<double>& f, std::vector<double>& dB,
  std::vector<double>& degrees, bool logFreq = true);

/** Plots the frequency response of the given filter. The class must have a function
getTransferFunctionAt... */
template<class TSig, class TFlt>
inline void plotFrequencyResponse(TFlt &filter, int N, TSig fMin, TSig fMax, TSig fs, bool logFreq)
{
  // create w array (normalized radian frequencies):
  std::vector<TSig> w(N);
  if(logFreq)
    RAPT::rsArray::fillWithRangeExponential(&w[0], N, fMin, fMax);
  else
    RAPT::rsArray::fillWithRangeLinear(&w[0], N, fMin, fMax);
  RAPT::rsArray::scale(&w[0], N, 2*PI/fs);

  // compute magnitude and phase response:
  std::vector<complex<TSig>> H = getFrequencyResponse(filter, w);
  std::vector<TSig> dB(N), phs(N);
  for(int k = 0; k < N; k++) {
    dB[k]  = RAPT::rsAmpToDb(abs(H[k]));
    phs[k] = arg(H[k]); //-2*PI; // arg is in -pi..+pi, we want -2*pi..0 - check, if this is correct
  }

  // unwrap phase, convert to degrees:
  RAPT::rsArray::unwrap(&phs[0], N, 2*PI);
  for(int k = 0; k < N; k++)
    phs[k] *= 180.0/PI;

  // maybe move the two steps above into rapt, too

  // convert w back to Hz and plot:
  RAPT::rsArray::scale(&w[0], N, fs/(2*PI));
  plotFrequencyResponse(w, dB, phs, logFreq);
}




// new, dragged over from RSLib tests (TestUtilities.h):

template<class T>
void plotArrays(int N, T *y1, T *y2 = nullptr, T *y3 = nullptr, T *y4 = nullptr, T *y5 = nullptr,
  T *y6 = nullptr, T *y7 = nullptr, T *y8 = nullptr, T *y9 = nullptr);

// convenience functions for interfacing with the Plotter class (they manage instantiation and
// setup of Plotter objects for frequently used cases):
void plotData(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL,
  double *y4 = NULL, double *y5 = NULL);

void plotData(int N, double x0, double dx, double *y1, double *y2 = NULL, double *y3 = NULL,
  double *y4 = NULL, double *y5 = NULL);
// for equidistant data, abscissa values start at x0 with increment dx

void plotDataLogX(int N, double *x, double *y1, double *y2 = NULL, double *y3 = NULL,
  double *y4 = NULL, double *y5 = NULL);

void plotVector(std::vector<double> v);

/** Plots the magnitude spectrogram given in s against time axis t (of length numFrames) and
frequency axis f (of length numBins). */
void plotSpectrogram(int numFrames, int numBins, double **decibels, double sampleRate,
  int hopSize, double dbMin = -100, double dbMax = +10);
// introduce parameters to control scaling of time- and frequency axis..

void plotPhasogram(int numFrames, int numBins, double **phases, double sampleRate,
  int hopSize);



// various convenience functions to plot filter responses for b/a specifications:
void plotMagnitudeResponse(const RAPT::rsFilterSpecificationBA<double>& specBA);
void plotPolesAndZeros(    const RAPT::rsFilterSpecificationBA<double>& specBA);
void showFilterPlots(      const RAPT::rsFilterSpecificationBA<double>& specBA);

/** Plots y against x using stems, i.e impulses with a filled circle - suitable to draw discrete 
time signals. */
void stemPlot(int N, double *x, double *y);


/** Convenience function. Uses class SinusoidalModelPlotter. */
void plotSineModel(const SinusoidalAnalyzer<double>& sa, double* sampleData, int N, 
  double sampleRate);


#endif
