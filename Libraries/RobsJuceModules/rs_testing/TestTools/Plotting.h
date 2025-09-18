#ifndef RAPT_PLOTTING_H
#define RAPT_PLOTTING_H

// This file contains convenience functions to create certain types of plots.



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


void setToDarkMode(GNUPlotter* plt);


/** Creates a time-axis (in seconds) given the sample-rate. */
void createTimeAxis(int numSamples, float *timeAxis, float sampleRate);
void createTimeAxis(int numSamples, double *timeAxis, double sampleRate);



/** Plots N samples of the impulse response of the passed filter. */
template<class TSig, class TFlt>
inline void plotImpulseResponse(TFlt &filter, int length, TSig scale, bool dB = false, 
  TSig dbFloor = TSig(-100))
{
  std::vector<TSig> y = impulseResponse(filter, length, scale);
  if(dB)
  {
    TSig floorAmp = rsDbToAmp(dbFloor);
    for(size_t i = 0; i < y.size(); i++)
      y[i] = rsAmpToDb(rsMax(rsAbs(y[i]), floorAmp));
  }
  GNUPlotter plt;
  plt.addDataArrays(length, &y[0]);
  plt.plot();
}

/** Takes a function object "transferFunc" and a vector of normalized radian frequencies w (aka 
"omega") and returns the corresponding vector of frequency response values H(e^(j*w)). The 
"transferFunc" function object is supposed to compute the transfer function of a digital filter 
for a given complex number z = e^(j*w) where w = omega = 2*pi*frequency/sampleRate and j is the 
imaginary unit. */
template<class T, class TFunc>
inline std::vector<std::complex<T>> getFreqRespFromTransFunc(
  const TFunc& transferFunc, const std::vector<T>& w)
{
  size_t N = w.size();
  std::complex<T> j(0,1);                  // Imaginary unit
  std::vector<std::complex<T>> H(N);       // H(e^jw)
  for(size_t k = 0; k < N; k++)
    H[k] = transferFunc(exp(j*w[k]));
  return H;
}
// Maybe rename to getDigitalFrequencyResponse. ...but maybe not. We are on a computer - we are
// on a compute here, so filters being digital is a given.
//
// Maybe rename to getFreqRespFromTransFunc because the overload resolution with the function
// below is based solely on the constness of the 1st parameter. That is asking for trouble. We 
// want to distinguish the function by name.  ...done!

/** Takes an arbitrary filter object that implements a getTransferFunctAt() member function and 
produces the filter's frequency response at the given vector of normalized radian frequencies 
w = 2*pi*frequency/sampleRate. The getTransferFunctAt() function must take a complex number z
and produce the transfer function H(z) at the given value z. */
template<class T, class TFlt>
inline std::vector<std::complex<T>> getFrequencyResponse(
  TFlt &filter, const std::vector<T>& w)
{
  return getFreqRespFromTransFunc(
    [&](std::complex<T> z) { return filter.getTransferFunctionAt(z); }, w);
}
// Maybe move to RAPT...but maybe use plain arrays instead of vectors there, keep convenience
// function here

/** Plots the given magnitude response in dB and phase response in degrees against the frequency
axis f. */
void plotFrequencyResponse(std::vector<double>& f, std::vector<double>& dB,
  std::vector<double>& degrees, bool logFreq = true);

template<class T>
std::vector<T> getOmegas(int N, T fMin, T fMax, T fs, bool logFreq)
{
  std::vector<T> w(N);
  if(logFreq) RAPT::rsArrayTools::fillWithRangeExponential(&w[0], N, fMin, fMax);
  else        RAPT::rsArrayTools::fillWithRangeLinear(&w[0], N, fMin, fMax);
  RAPT::rsArrayTools::scale(&w[0], N, 2*PI/fs);
  return w;
}


/** Plots the frequency response of the given "transferFunc". This must be a function object just 
like in getFrequencyResponse(const TFunc& transferFunc, const std::vector<T>& w). */
template<class TSig, class TFunc>
inline void plotFreqRespFromTransFunc(
  const TFunc &transferFunc, int N, TSig fMin, TSig fMax, TSig fs, bool logFreq)
{
  // Create w array (normalized radian frequencies):
  std::vector<TSig> w = getOmegas(N, fMin, fMax, fs, logFreq);

  // Compute magnitude and phase response:
  std::vector<std::complex<TSig>> H = getFreqRespFromTransFunc(transferFunc, w);
  std::vector<TSig> dB(N), phs(N);
  for(int k = 0; k < N; k++) {
    dB[k]  = RAPT::rsAmpToDb(abs(H[k]));
    phs[k] = arg(H[k]); //-2*PI; // arg is in -pi..+pi, we want -2*pi..0 - check, if this is correct
  }

  // Unwrap phase, convert to degrees:
  RAPT::rsArrayTools::unwrap(&phs[0], N, 2*PI);
  for(int k = 0; k < N; k++)
    phs[k] *= 180.0/PI;

  // maybe move the two steps above into rapt, too

  // Convert w back to Hz and plot:
  RAPT::rsArrayTools::scale(&w[0], N, fs/(2*PI));
  plotFrequencyResponse(w, dB, phs, logFreq);
}
// Maybe rename to plotFreqRespFromTransFunc

/** Plots the frequency response of the given filter. The class TFlt must have a function
getTransferFunctionAt()... TBC... */
template<class T, class TFlt>
inline void plotFrequencyResponse(TFlt& filter, int N, T fMin, T fMax, T fs, bool logFreq)
{
  plotFreqRespFromTransFunc(
    [&](std::complex<T> z) { return filter.getTransferFunctionAt(z); }, 
    N, fMin, fMax, fs, logFreq);
}

/** Takes a vector f of frequencies for the x-axis and 2 vectors r1, r2 of frequency responses for
the y-values and plots them in a single plot. The two frequency responses r1, r2 can be responses 
of two filters or different aspects of the response of a single filter, such as real and imaginary
part, magnitude and phase, etc. */
template<class T>
void plotTwoFrequencyResponses(
  std::vector<T>& f, std::vector<T>& r1, std::vector<T>& r2, bool logFreq = true);


template<class TSig, class TFlt>
inline void plotFrequencyResponseReIm(TFlt& filter, int N, TSig fMin, TSig fMax, TSig fs, bool logFreq)
{
  std::vector<TSig> w = getOmegas(N, fMin, fMax, fs, logFreq);
  std::vector<std::complex<TSig>> H = getFrequencyResponse(filter, w);
  RAPT::rsArrayTools::scale(&w[0], N, fs/(2*PI));
  std::vector<TSig> re(N), im(N);
  for(int k = 0; k < N; k++) {
    re[k] = H[k].real();
    im[k] = H[k].imag(); }

  plotTwoFrequencyResponses(w, re, im, logFreq);
}


/** Magnitude- and ringing response. The concept of "ringing response" is still VERY sketchy and 
experimental. See this thread: https://www.kvraudio.com/forum/viewtopic.php?f=33&t=569114 */
template<class TSig, class TFunc>
inline void plotMagAndRingRespFromTransFunc(const TFunc& transferFunc, int N, TSig fMin, TSig fMax, 
  TSig fs, bool logFreq, bool ringingQ = true)
{
  // This is the same as in function above - maybe factor out:
  std::vector<TSig> w = getOmegas(N, fMin, fMax, fs, logFreq);
  std::vector<std::complex<TSig>> H = getFreqRespFromTransFunc(transferFunc, w);
  RAPT::rsArrayTools::scale(&w[0], N, fs/(2*PI));
  std::vector<TSig> re(N), im(N);
  for(int k = 0; k < N; k++) {
    re[k] = H[k].real();
    im[k] = H[k].imag(); }

  // Compute numeric derivatives of real and imaginary parts:
  using ND = rsNumericDifferentiator<TSig>;
  std::vector<TSig> dre(N), dim(N);          // derivatives of re and im
  ND::derivative(&w[0], &re[0], &dre[0], N);
  ND::derivative(&w[0], &im[0], &dim[0], N);

  // Compute magnitudes of frequency response and its derivative:
  std::vector<TSig> mag(N), dmag(N);
  for(int k = 0; k < N; k++) 
  {
    mag[k]  = sqrt( re[k]* re[k] +  im[k]* im[k]);
    dmag[k] = sqrt(dre[k]*dre[k] + dim[k]*dim[k]);

    if(ringingQ)
      dmag[k] *= w[k];
    // This makes plots of some filters (like elliptic bandpass) symmetric. I think, when we do 
    // that, we get a relative ringing time, i.e. expressed in number-of-cycles instead of in 
    // seconds. 

    //dmag[k] *= dmag[k];
    // test - undo sqrt...hmm...nope

    //dmag[k] /= mag[k];
    //dmag[k] /= (1 + mag[k]);  // ad-hoc remedy against the division by zero
    // divide by magnitude to make measure independent of magnitude. This seems to work well for 
    // allpole filters but has problems when we have zeros on the imaginary axis (like in notch 
    // filters, elliptic filters, etc.). This is not surprising because it implies a division by 
    // zero.
  }

  // Plot:
  plotTwoFrequencyResponses(w, mag, dmag, logFreq);

  // Observations: 
  //
  // - Multiplying dmag[k] by w[k] makes the response symmetrical for ellitpic bandpasses, but 
  //   maybe the asymmetry is actually a legit feature because from the impulse response, it seems
  //   like the lower bandedge does indeed ring longer. Maybe the multiplication by w[k] gives the
  //   number of cycles of ringing, not the absolute time?
  //
  // - I think, the plotting code normalizes the data internally. The plots both hit 1.0 but I 
  //   think the data does not necessarily. ...hmm - nope - sometimes, it does not hi 1.0. That 
  //   seems to be a coincidence that happens in some cases.
  //
  //
  // ToDo:
  //
  // - Try if it makes a difference, if we do the multiplication by w[k] before computing the 
  //   magnitude or even before computing the derivative
  //
  // - Verify if the name of ringingQ is appropriate. The idea is to normalize the ringing time by
  //   the frequency - or someting. I'm not quite sure about the rationale of multiplying the 
  //   magnitude of the complex derivative by the frequency. I just observed that for some filters,
  //   doing so makes the plots look nicely symmetric. I think, without the multiplication, the 
  //   plot shows the absolute ringing time. For higher frequencies, the same absolute ringing time
  //   implies a higher Q. The Q is proportional to the frequency. So, yeah - the name seems 
  //   appropriate.
  //
  // - Factor out a getRingingResponse function. We could actually implement the numerical 
  //   differentiation differently (more accurately) by using the numeric differentiation functions
  //   that evaluate the function at x+h, x-h. We don't nee to use the data based estimator. We can
  //   use the function based estimator.
}

template<class T, class TFlt>
inline void plotMagAndRingResponse(TFlt& filter, int N, T fMin, T fMax, T fs, 
  bool logFreq, bool ringingQ = true)
{
  plotMagAndRingRespFromTransFunc(
    [&](std::complex<T> z) { return filter.getTransferFunctionAt(z); }, 
    N, fMin, fMax, fs, logFreq, ringingQ);
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
void plotComplexVectorReIm(std::vector<rsComplex<double>> v);
// ToDo: take argument by const reference




template<class T>
void rsPlotComplexPoints(const std::vector<std::complex<T>>& points)
{
  GNUPlotter plt;
  plt.addDataComplex(points);
  plt.setToDarkMode();
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");
  //plt.addGraph("i 0 u 1:2 w points pt 7 ps 0.6 notitle");
  plt.addGraph("i 0 u 1:2 w points pt 7 ps 0.4 notitle");
  plt.plot();
  // It looks a bit ugly - as if the point locations are rounded to the nearest pixel or something.
  // Factor out into a rsPlotComplexPoints function
}



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

template<class T>
inline void plotMatrix(const rsMatrix<T>& A, bool asHeatMap = true)
{
  GNUPlotter plt;

  plt.addDataMatrixFlat(A.getNumRows(), A.getNumColumns(), A.getRowPointerConst(0));

  plt.setToDarkMode();
  plt.setPixelSize(800, 800);

  if(asHeatMap) 
  {
    //plt.addCommand("set size square");  // make optional

    plt.addGraph("i 0 nonuniform matrix w image notitle");

    // CJ_BuYlRd11

    plt.setColorPalette(GNUPlotter::ColorPalette::CJ_BuYlRd11, false);
    //plt.setColorPalette(GNUPlotter::ColorPalette::ML_Parula, false);
    //plt.setColorPalette(GNUPlotter::ColorPalette::SW_Inferno, false);
    //plt.addCommand("set palette gray");
    //plt.addCommand("set palette rgbformulae 7,5,15");
    // http://gnuplot.info/demo_5.2/pm3dcolors.html


    //plt.setRange(0, A.getNumRows()-1, 0, A.getNumColumns()-1, -1.0, +1.0); // doesn't work

    if(A.isSquare())
      plt.addCommand("set size square");

    plt.plot();
  }
  else
  {
    plt.addCommand("set xlabel \"X\" textcolor rgb \"black\"");
    plt.addCommand("set ylabel \"Y\" textcolor rgb \"black\"");
    plt.plot3D();
  }
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
void plotMatrixWithMarkers(const rsMatrix<T>& A, const std::vector<int>& markers)
{
  GNUPlotter plt;

  //rsMatrix<T> B = A;  // Try to get rid. See comment in plotMatrix() for how to
  //plt.addDataMatrixFlat(B.getNumRows(), B.getNumColumns(), B.getRowPointer(0));
  plt.addDataMatrixFlat(A.getNumRows(), A.getNumColumns(), A.getRowPointerConst(0));
  plt.addGraph("i 0 nonuniform matrix w image notitle");

  // Set style options:
  plt.setToDarkMode();
  plt.setPixelSize(600, 600);
  //plt.setColorPalette(GNUPlotter::ColorPalette::CJ_BuYlRd11, false);
  plt.setColorPalette(GNUPlotter::ColorPalette::ML_Parula, false);
  if(A.isSquare())
    plt.addCommand("set size square");

  // Draw the markers:
  std::string str;
  for(size_t i = 0; i < markers.size(); i++)
  {
    str  = "set object " + std::to_string(i+1) + " circle front at ";
    str += std::to_string(i) + "," + std::to_string(markers[i]);
    str += " size 0.15 fillcolor rgb \"black\" fs solid";
    plt.addCommand(str);
  }

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
void plotBivariateFunction(const std::function<T(T, T)>& f,
  T minX, T maxX, int numX, T minY, T maxY, int numY)
{
  GNUPlotter plt;
  plt.plotBivariateFunction(numX, minX, maxX, numY, minY, maxY, f);
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
// numFrames and numBins can be inferred from the matrix shape -> get rid!

/** Under construction.... */
void plotSpectrogram(const double* x, int N, int hopSize, int blockSize, int trafoSize, 
  double sampleRate);





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

template<class T>
inline void addDataFunction(GNUPlotter& plt, const std::function<T(T)>& f, T xMin, T xMax, int N)
{
  std::vector<T> x = rsRangeLinear(xMin, xMax, N);

  // Maybe factor out into y = rsApply(f, x):
  std::vector<T> y(N);
  for(int i = 0; i < N; i++)
    y[i] = f(x[i]);

  plt.addDataArrays(N, &x[0], &y[0]);
}
// This function could go into class GNUPlotter itself


/** Plots the content of the two given delaylines. The pruposes of this function is to be able to
check if two delaylines have the same or related content such that one may think about getting rid 
of one of them. */
template<class T>
void rsPlotDelayLineContent(const RAPT::rsDelay<T>& dl1, 
  const RAPT::rsDelay<T>& dl2, bool reverseSecond = false)
{
  // Helper function to return the content of the given delayline as std::vector:
  auto getContent = [](const RAPT::rsDelay<T>& dl)
  {
    // Maybe let the use switch between showing the full content (i.e. the full allocated memory)
    // or only up to the used length - current, we hrdcoded the used length:
    //int N = dl.getMaxDelayInSamples();
    int N = dl.getDelayInSamples();
    std::vector<T> cnt(N);
    for(int n = 0; n < N; n++)
      cnt[n] = dl.readOutputAt(n);
    return cnt;
  };

  // Retrieve contents of both delaylines and plot them:
  std::vector<T> cnt1 = getContent(dl1);
  std::vector<T> cnt2 = getContent(dl2);
  if(reverseSecond)
    rsReverse(cnt2);
  rsPlotVectors(cnt1, cnt2);

  // ToDo:
  //
  // - Done. -> Document it!
  //   Maybe optionally reverse the content of the 2nd delay line. This may be convenient for
  //   experimenting with implementing bidirectional delay lines 
}



#endif
