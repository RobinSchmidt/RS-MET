#ifndef RAPT_PLOTTING_H
#define RAPT_PLOTTING_H

// why do we need these includes here?
#include "rapt/rapt.h"
#include "rosic/rosic.h"
#include "rs_testing/rs_testing.h"

//#include "../Prototypes/SinusoidalModeling.h"
//#include "../../../../Tests/TestsRosicAndRapt/Source/Shared/Prototypes/SinusoidalModeling.h"
// get rid - class should be moved to rapt - done

///** Plots at most five y-functions against a common x-axis. */
//void plotData(int N, float *x, float *y1, float *y2 = nullptr, float *y3 = nullptr,
//  float *y4 = nullptr, float *y5 = nullptr);



/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate);
void createTimeAxis(int numSamples, double *timeAxis, double sampleRate);


/** Returns N samples of the impulse response of the passed filter as std::vector. It is necessary
for you to pass a scale factor of the type of the filter's output signal (for example: 1.0 for
double), such that the compiler can deduce the template parameter. We also use it to scale the
input impulse to the filter, so it actually gets some purpose besides satisfying the compiler. The
filter class must support the functions reset() and getSample() */
template<class TSig, class TFlt>
inline std::vector<TSig> impulseResponse(TFlt &filter, int length, TSig scale)
{
  std::vector<TSig> y(length);
  filter.reset();
  y[0] = filter.getSample(scale);
  for(int n = 1; n < length; n++)
    y[n] = filter.getSample(0.0);
  return y;
}
template<class TSig, class TFlt>
inline std::vector<TSig> filterResponse(TFlt& filter, int length, std::vector<TSig> x)
{
  std::vector<TSig> y(length);
  filter.reset();
  for(int n = 0; n < length; n++)
    y[n] = filter.getSample(x[n]);
  return y;
}
template<class T>
inline std::vector<T> ampToDb(const std::vector<T>& x, T minDb)
{
  std::vector<T> y(x.size());
  for(size_t i = 0; i < x.size(); i++)
    y[i] = rsMax(rsAmpToDb(rsAbs(x[i])), minDb);
  return y;
}


// move to Utilities



/** Plots N samples of the impulse response of the passed filter. */
template<class TSig, class TFlt>
inline void plotImpulseResponse(TFlt &filter, int length, TSig scale)
{
  std::vector<TSig> y = impulseResponse(filter, length, scale);
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
  std::vector<std::complex<TSig>> H(N);  // H(e^jw)
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
    RAPT::rsArrayTools::fillWithRangeExponential(&w[0], N, fMin, fMax);
  else
    RAPT::rsArrayTools::fillWithRangeLinear(&w[0], N, fMin, fMax);
  RAPT::rsArrayTools::scale(&w[0], N, 2*PI/fs);

  // compute magnitude and phase response:
  std::vector<std::complex<TSig>> H = getFrequencyResponse(filter, w);
  std::vector<TSig> dB(N), phs(N);
  for(int k = 0; k < N; k++) {
    dB[k]  = RAPT::rsAmpToDb(abs(H[k]));
    phs[k] = arg(H[k]); //-2*PI; // arg is in -pi..+pi, we want -2*pi..0 - check, if this is correct
  }

  // unwrap phase, convert to degrees:
  RAPT::rsArrayTools::unwrap(&phs[0], N, 2*PI);
  for(int k = 0; k < N; k++)
    phs[k] *= 180.0/PI;

  // maybe move the two steps above into rapt, too

  // convert w back to Hz and plot:
  RAPT::rsArrayTools::scale(&w[0], N, fs/(2*PI));
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

void plotVector(std::vector<double> v);  // replace by RAPT::rsPlotVector

void plotComplexVectorReIm(std::vector<std::complex<double>> v);


/** Plots the matrix entries as surface above a coordinate system given by x,y 
todo: check/assert that the dimensions of the matrix z fit together with the lengths of x,y 
// try to make inputs const */
template<class T>
//void plotMatrix(const RAPT::rsMatrix<T>& z, const std::vector<T>& x, const std::vector<T>& y)
void plotMatrix(RAPT::rsMatrix<T>& z, std::vector<T>& x, std::vector<T>& y)
{
  double** z2;
  RAPT::rsMatrixTools::allocateMatrix(z2, z.getNumRows(), z.getNumColumns());

  for(int i = 0; i < z.getNumRows(); i++)
    for(int j = 0; j < z.getNumColumns(); j++)
      z2[i][j] = z(i,j);

  GNUPlotter plt;
  plt.plotSurface((int)x.size(), (int)y.size(), &x[0], &y[0], z2);

  RAPT::rsMatrixTools::deallocateMatrix(z2, z.getNumRows(), z.getNumColumns());
}
// get rid of that - use function below instead - maybe it should take optional x,y arguments

inline void plotMatrix(rsMatrix<double>& A, bool asHeatMap)  // use const
{
  GNUPlotter plt;
  //plt.addDataMatrixFlat( A.getNumRows(), A.getNumColumns(), A.getDataPointerConst());
  plt.addDataMatrixFlat( A.getNumRows(), A.getNumColumns(), A.getRowPointer(0));
  if(asHeatMap) {
    plt.addCommand("set size square");  // make optional
    plt.addGraph("i 0 nonuniform matrix w image notitle");
    plt.addCommand("set palette gray");
    //plt.setRange(0, A.getNumRows()-1, 0, A.getNumColumns()-1, -1.0, +1.0); // doesn't work
    plt.plot();
  }
  else
    plt.plot3D();
}
// un-inline


/** Plots the rows of the matrix against the column-index, i.e. the rows are interpreted as 
several functions of x where x is the column index. */
template<class T>
void plotMatrixRows(const RAPT::rsMatrix<T>& A)
{
  GNUPlotter plt;
  for(int i = 0; i < A.getNumRows(); i++)
    plt.addDataArrays(A.getNumColumns(), A.getRowPointerConst(i));
  plt.plot();
}

/** Like the function above but with custom x-axis. The length of x should be equal to the number 
of columns in the matrix. */
template<class T>
void plotMatrixRows(const RAPT::rsMatrix<T>& A, T* x)
{
  GNUPlotter plt;
  for(int i = 0; i < A.getNumRows(); i++)
    plt.addDataArrays(A.getNumColumns(), x, A.getRowPointerConst(i));
  plt.plot();
}

template<class T>
void plotPolynomial(const T* a, int degree, T min, T max, int numPoints = 200)
{
  std::vector<T> x(numPoints), y(numPoints);
  rsArrayTools::fillWithRangeLinear(&x[0], numPoints, min, max);
  for(int i = 0; i < numPoints; i++)
    y[i] = rsPolynomial<T>::evaluate(x[i], a, degree);
  rsPlotArraysXY(numPoints, &x[0], &y[0]);
}

template<class T>
void plot(const rsPiecewisePolynomial<T>& p, T xMin, T xMax, int numSamples)
{
  std::vector<T> x(numSamples), y(numSamples);
  rsArrayTools::fillWithRangeLinear(&x[0], numSamples, xMin, xMax);
  for(int i = 0; i < numSamples; i++)
    y[i] = p.evaluate(x[i]);
  rsPlotVectorsXY(x, y);
}

template<class T>
void plot(const rsPiecewisePolynomial<T>& p, int numSamples = 501)
{
  plot(p, p.getDomainMinimum(), p.getDomainMaximum(), numSamples);
}

template<class T>
void plotBivariatePolynomial(const rsBivariatePolynomial<T> p,
  T minX, T maxX, int numX, T minY, T maxY, int numY)
{
  GNUPlotter plt;
  std::function<T(T, T)> f = [&](T x, T y) -> T { return p(x, y); };
  plt.plotBivariateFunction(numX, minX, maxX, numY, minY, maxY, f);
}

/** Plots the magnitude spectrogram given in s against time axis t (of length numFrames) and
frequency axis f (of length numBins). */
void plotSpectrogram(int numFrames, int numBins, double **decibels, double sampleRate,
  int hopSize, double dbMin = -100, double dbMax = +10);
// introduce parameters to control scaling of time- and frequency axis..

/** Plots spectrogram magnitudes from a complex spectrogram. */
void plotSpectrogram(int numFrames, int numBins, const rsMatrix<std::complex<double>>& spec,
  double sampleRate, int hopSize, double dbMin = -100, double dbMax = +10);



void plotPhasogram(int numFrames, int numBins, double **phases, double sampleRate,
  int hopSize);



// various convenience functions to plot filter responses for b/a specifications:
void plotMagnitudeResponse(const RAPT::rsFilterSpecificationBA<double>& specBA);
void plotPolesAndZeros(    const RAPT::rsFilterSpecificationBA<double>& specBA);
void showFilterPlots(      const RAPT::rsFilterSpecificationBA<double>& specBA);

/** Plots y against x using stems, i.e impulses with a filled circle - suitable to draw discrete
time signals. */
void stemPlot(int N, double *x, double *y);


// functions for plotting sinusoidal model data:

/** Convenience function. Uses class SinusoidalModelPlotter. */
void plotSinusoidalAnalysisResult(RAPT::rsSinusoidalAnalyzer<double>& sa, double* sampleData, int N,
  double sampleRate);

void plotSineModel(const RAPT::rsSinusoidalModel<double>& model, double sampleRate);

void plotTwoSineModels(
  const RAPT::rsSinusoidalModel<double>& model1,
  const RAPT::rsSinusoidalModel<double>& model2,
  double sampleRate);

/** Plots a subset of the amplitude envelopes of the given sinusoidal model. The vector
partialIndices selects, which partial's envelopes should be drawn. */
void plotSineModelAmplitudes(
  const RAPT::rsSinusoidalModel<double>& model,
  std::vector<int> partialIndices = std::vector<int>());

/** Plots a subset of the unwrapped phases of the model - but because the phases themselves are not
that useful to look at (you would basically just see an upward sloping line), the function plots
either a de-trended version of the phases or the phase-derivative, depending on the boolean
parameter - if false, de-trended phases are plotted, if true, phase-derivatives are plotted. */
void plotSineModelPhases(
  const RAPT::rsSinusoidalModel<double>& model,
  const std::vector<int>& partialIndices, bool derivative = false);


/** Plots the results of a sinusoidal synthesis of given model using given synthesizer object and
plots also the original signal x (of length N) and the synthesis error. */
void plotSineResynthesisResult(
  const RAPT::rsSinusoidalModel<double>& model,
  const RAPT::rsSinusoidalSynthesizer<double>& synth,
  double* x, int N);
// todo: make x and n optional arguments

void plotModelOutputComparison(
  const RAPT::rsSinusoidalModel<double>& model1,
  const RAPT::rsSinusoidalModel<double>& model2,
  const RAPT::rsSinusoidalSynthesizer<double>& synth);

void plotModalAmplitudes(const std::vector<rsModalFilterParameters<double>>& modelModel);

/** Plots the amplitude envelope of a sinusoidal partial and the amplitude envelope of a modal
model of a partial for comparison. */
void plotModeVsSineAmpEnv(
  rsModalFilterParameters<double>& modal, RAPT::rsSinusoidalPartial<double>& sinusoidal);

/** Adds partial data from x- and y-arrays - only those datapoints that are listed in indices 
array. */
inline void addDataPartially(GNUPlotter& plt, 
  const std::vector<double>& x, const std::vector<double>& y, const std::vector<int> indices)
{
  rsAssert(x.size() == y.size());
  int M = int(indices.size());
  std::vector<double> xt(M), yt(M);           // temporary arrays to hold partial data
  for(int m = 0; m < M; m++)  {
    int n = indices[m];                       // index into x- and y-arrays
    rsAssert(n >= 0 && n < (int) x.size());
    xt[m] = x[n]; 
    yt[m] = y[n]; }
  plt.addDataArrays(M, &xt[0], &yt[0]);
  // maybe use rsSelect - but that uses an array of type size_t for the indices - maybe make 
  // another version that uses int
}



#endif
