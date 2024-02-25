#ifndef DSPPLOTTERS_H
#define DSPPLOTTERS_H

//#include "GNUPlotter.h"
//#include "../../../../Libraries/JUCE/modules/rosic/rosic.h"
#include "../../../Libraries/RobsJuceModules/rosic/rosic.h"

// rename file to Plotters.h - we have also other plotters that are not necessarily related to 
// DSP, for example for vector fields. Maybe SpecializedPlotters

/* Subclasses of GNUPlotter that specialize in making plots related to digital signal processing 
(DSP). */


//=================================================================================================

/** A class for visualizing analog and digital filter responses. It may plot magnitude responses,
phase responses, pole/zero plots etc. for one or more filters. You need to specify the filters in 
terms of their poles and zeros. For each filter, you once call addPoleZeroSet to add the filter to
the list. Once you are finished adding filters this way, you can get the various plots via the
respective plot... functions. */

template <class T>
class FilterPlotter : public GNUPlotter
{

public:

  static const T pi;  // 3.1415926535897932384626433832795
  static const T inf; // std::numeric_limits<T>::infinity()

  FilterPlotter();

  /** Decides whether frequency parameters that are passed in to certain functions should be
  interpreted as radian frequencies (typically denoted as omega) or not. In the latter case, they
  are interpreted as being in Hz for analog filters or as fractions of the sample rate for 
  digital filters. True by default, i.e. if you don't call this function, they are interpreted as 
  radian frequencies. */
  void setFrequenciesAreRadian(bool areRadian);

  void setDecibelFloor(T newFloor) { dBFloor = newFloor; }


  /** Adds a filter specification in terms of poles, zeros and gain to our list. You may also pass 
  a sampleRate in which case the poles and zeros will be interpreted as z-plane values. Otherwise
  an analog filter (corresponding to an infinite sample rate) will be assumed and the poles and 
  zeros are interpreted as being in the s-plane */
  void addFilterSpecificationZPK(int numPoles, std::complex<T>* poles, int numZeros, 
    std::complex<T>* zeros,  T gain, T sampleRate = inf);

  /** Adds a filter specification via an object the rsFilterSpecificationZPK class. */
  void addFilterSpecificationZPK(const RAPT::rsFilterSpecificationZPK<T>& spec);


  void addFilterSpecificationBA(int numeratorOrder, const T* numeratorCoeffs, 
    int denominatorOrder, const T* denominatorCoeffs, T sampleRate = inf);
  // maybe allow for complex coeffs

  void addFilterSpecificationBA(const RAPT::rsFilterSpecificationBA<T>& spec);

  /** Adds a transfer function. For an analog filter, you should pass inf for the sampleRate.
  It will be internally converted inot BA and ZPK representations. */
  void addTransferFunction(const RAPT::rsRationalFunction<T>& tf, T sampleRate);



  /** Plots the magnitude responses of all the filters. */
  void plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels);
  // todo: make this a convenience function that calls plotFrequencyResponse with false for the 
  // plotPhase flag, make a similar function for plotPhase

  /** Plots the frequency responses of all the filters. */
  void plotFrequencyResponses(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, 
    bool plotMagnitude = true, bool decibels = true,
    bool plotPhase = true, bool unwrapPhase = true);

  /** Plots the poles and zeros of all the filters in s- or z-plane. You can pass the plot-size in 
  pixels (it will always use a square plot). */
  void plotPolesAndZeros(int plotSize = 400);

  /*
  void plotMagnitudeAndPhase(); // in one plot
  void plotPhaseDelay();
  void plotGroupDelay();
  void plotTransferFunctionMagnitude();
  void plotImpulseResponse();
  void plotStepResponse();
  */

  /** Creates a vector of x-values for the frequency axis either linearly or logarithmically 
  scaled. */
  std::vector<T> getFrequencyAxis(int numFreqs, T lowFreq, T highFreq, bool logarithmic);

  /** Returns the complex frequency response of the filter with given index at the frequencies 
  given in the vector. */
  std::vector<std::complex<T>> getFrequencyResponse(int index, const std::vector<T>& frequencies);

  /** Extracts the magnitudes from the passed complex frequency response array.  */
  std::vector<T> getMagnitudes(const std::vector<std::complex<T>>& complexFreqResponse, 
    bool inDecibels = false);

  /** Extracts the phases from the passed complex frequency response array.  */
  std::vector<T> getPhases(const std::vector<std::complex<T>>& complexFreqResponse, 
    bool unwrap = true, bool inDegrees = true);

  /** Evaluates polynomial defined by its roots at the value z. */
  std::complex<T> polynomialByRoots(std::complex<T> z, std::vector<std::complex<T>>& roots);

  /** Evaluates complex transfer function defined by its zeros z, poles p and gain k at the 
  complex value s */
  std::complex<T> transferFunctionZPK(std::complex<T> s, std::vector<std::complex<T>>& z,
    std::vector<std::complex<T>>& p, std::complex<T> k);

  //-----------------------------------------------------------------------------------------------

  //static FilterSpecificationBA<T> zpk2ba(const FilterSpecificationZPK<T>& zpk);
  // moved to rsFilterSpecificationBA

  /** Converts a filter specification in terms of numerator and denominator polynomial coefficients
  to a specification in terms of zeros, poles and gain. */
  //static FilterSpecificationZPK<T> ba2zpk(const FilterSpecificationBA<T>& ba);

  /** Normalizes the a0 coefficient in a filter specification to unity. */
  //static void normalizeA0(FilterSpecificationBA<T>& ba);  // moved to rapt

protected:

  /** Adds the commands to set up the appropriate plotting options for for a pole/zero plot. */
  void setupForPoleZeroPlot(int size);

  /** Given an array of complex values z (for example, roots of a polynomial), this function plots
  their multiplicities at their positions */
  void drawMultiplicities(const std::vector<std::complex<T>>& z, T thresh);
    // not yet tested

  /** Draws vertical lines at the Nyquist frequencies of the various filters. If the filter
  specifications happen to have all the same sample rate, only a single black line is drawn. 
  Otherwise, the colors of the lines are matched to the colors of the graphs to which they 
  apply. */
  void drawNyquistLines();

  /** Returns maximum absolute value of all real an imaginary parts. */
  double maxAbsReIm(const std::vector<std::complex<T>>& x);

  /** Returns true, if the relative distance between x and y is smaller than the given threshold 
  ("relative" with respect to the actual absolute values of x and y, such that for larger values 
  the tolerance also increases) */
  bool almostEqual(std::complex<T> x, std::complex<T> y, T thresh);

  /** Adds the command for a line graph to the commandfile, i.e. for something like a magnitude- or 
  phase response. The yAxis2 flag indicates, if the graph uses the second y-axis (on the right of 
  the graph). For example, when plotting magnitude- and phase responses into a single plot, the 
  magnitudes may use the left (1st) y-axis and the phases may use the right (2nd) y-axis. */
  void addGraphLines(int graphIndex, bool yAxis2 = false);

  std::vector<std::string> getGraphColors(int numGraphs) const;

  T freqScale = 1.0;
  T dBFloor   = T(-120);

  std::vector<RAPT::rsFilterSpecificationZPK<T>> filterSpecsZPK;
  std::vector<RAPT::rsFilterSpecificationBA<T>>  filterSpecsBA;

  // have conversion functions convert_BA_To_ZPK, convert_ZPK_To_BA
  // or filterSpecBA2ZPK, ba2zpk, zpk2ba
  // maybe keep tpk and ba specifications for each filter, i.e. keep the two representations in 
  // sync and use whatever representation is more convenient to deal with to compute the various
  // plots/numbers

  // maybe we could also keep SOS representations and maybe others?

};

// todo: allow for filters to be specified also in terms of their polynomial coefficients - maybe 
// either by incoroprating a root finder or by allowing to store a filter specification in one of
// the two alternative formats (zeros-poles-gain (zpk) or coeff-arrays), conversion from zpk to
// coeffs may be required anyway to compute group delay (i think we'll need derivatives of 
// numerator and denominator, quotient-rule and maybe chain-rule to evaluate the derivative of the
// phase-response)

//=================================================================================================

/** Class for plotting FFT spectra. You may pass several input signals as time-domain arrays and it
will plot the Spectra of these signals.

Under construction - may still have bugs

ToDo:
-Create example plots using the signals in "Understanding DSP" pg 54, 69. Check if freq-axis is 
 scaled correctly

*/

template <class T>
class SpectrumPlotter : public GNUPlotter
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Setup

  enum class FreqAxisUnits  // maybe rename to FreqAxisUnit
  {
    binIndex = 0,    // 0...N/2     (verify!)
    normalized,      // 0...0.5
    omega,           // 0...pi
    hertz            // 0...fs/2    (verify!)
  };


  /** The different modes are suitable for different types of input signals.  */
  enum class NormalizationMode
  {
    cycle,       // cycle (of a periodic signal)
    impulse,     // impulse response (of a filter/LTI-system)
    toZeroDb     // maximum is always forced to 0 dB
  };
  // what about (auto)correlation functions? how should we normalize the plot for them?


  /** Sets the FFT size. Does not have to be a power of 2 - we use Bluestein FFT here. */
  void setFftSize(int newSize) { fftSize = newSize; }
  // maybe rename to setTrafoSize

  void setFreqAxisUnit(FreqAxisUnits newUnit) { freqAxisUnit = newUnit; }

  void setLogFreqAxis(bool freqsAreLogarithmic) { logFreqAxis = freqsAreLogarithmic; };

  /** Sets the normalization mode. */
  void setNormalizationMode(NormalizationMode newMode) { normMode = newMode; }

  /** Sets the floor level for the plot in decibels. */
  void setFloorLevel(T newFloor) { dBFloor = newFloor; }

  /** Sets the sample rate - this affects the scaling of the frequency axis, if it's scaled in 
  Hz. */
  void setSampleRate(T newRate) { sampleRate = newRate; }


  //-----------------------------------------------------------------------------------------------
  // \Plotting

  /** Given up to 10 signal buffers of length "signalLength", this function performs an FFT on each 
  of them and plots the spectral magnitudes as decibel values. The FFT size is determined by 
  setFftSize and may be different from signalLength - if signalLength is shorter, the FFT buffers
  will be padded with zeros and if it's longer, only the leading sections of buffers will be 
  used. */
  void plotDecibelSpectra(int signalLength, const T *x0, const T *x1 = nullptr, 
    const T *x2 = nullptr, const T *x3 = nullptr, const T *x4 = nullptr, const T *x5 = nullptr, 
    const T *x6 = nullptr, const T *x7 = nullptr, const T *x8 = nullptr, const T *x9 = nullptr);

  /** Plots the spectra of the rows of the given matrix. Each row is taken to be a signal in the
  time domain, e.g. an impulse response. */
  void plotDecibelSpectraOfRows(const rsMatrix<T>& signals);

  /** Plots the spectra of the given signals. The first index is the signal, the second the sample
  index. All signals are assumed to have the same length. */
  void plotSpectra(const T** signals, int numSignals, int signalLength);
  // rename to plotDecibelSpectra



  // Under construction:
  void plotPhaseSpectra(int signalLength, const T *x0, const T *x1 = nullptr, 
    const T *x2 = nullptr, const T *x3 = nullptr, const T *x4 = nullptr, const T *x5 = nullptr, 
    const T *x6 = nullptr, const T *x7 = nullptr, const T *x8 = nullptr, const T *x9 = nullptr);


  void plotPhaseSpectra(const T** signals, int numSignals, int signalLength);


protected:



  void setupTransformer();


  std::vector<T> getFreqAxis(int numBins);

  FreqAxisUnits freqAxisUnit = FreqAxisUnits::binIndex;
  NormalizationMode normMode = NormalizationMode::cycle;

  T sampleRate = T(1);

  int fftSize = 2048;

  //bool plotPhases = true;

  bool logFreqAxis = false;

  //bool normalize   = true;

  T dBFloor = T(-120);

  RAPT::rsFourierTransformerBluestein<T> transformer;


};

//=================================================================================================

/** A class for plotting spectrograms */

template<class T>
class SpectrogramPlotter
{

public:


  //static void plotSpectrogram(int numFrames, int numBins, double **decibels, double sampleRate,
  //  int hopSize, double dbMin = -100, double dbMax = +10);

  void addSpectrogramData(GNUPlotter& plt, int numFrames, int numBins, T **decibels, 
    T sampleRate, int hopSize, T dbMin = -100, T dbMax = 0);


protected:

  T minFreq = 0, maxFreq = RS_INF(T);
  T minDb = -100, maxDb = 0;

};

//=================================================================================================

/** A class for plotting the analysis results of the sinusoidal model */

template<class T>
class SinusoidalModelPlotter : public SpectrogramPlotter<T>
{

public:

  //void plot(SinusoidalAnalyzer<T>& sa, T* sampleData, int N, T sampleRate);


  void addModelToPlot(const RAPT::rsSinusoidalModel<T>& model, GNUPlotter& plt, 
     T sampleRate, const std::string& graphColor);

  /** Plots the frequency tracks of a single sine model. */
  void plotModel(const RAPT::rsSinusoidalModel<T>& model, T sampleRate);

  /** Plots the frequency tracks of two sine models into one plot - useful for comparing original 
  and analyzed models. */
  void plotTwoModels(
    const RAPT::rsSinusoidalModel<T>& model1, 
    const RAPT::rsSinusoidalModel<T>& model2, 
    T sampleRate);
  // maybe have functions for plotting even more models - this can be useful to compare analysis
  // results with different parameter settings

  /**Plots results of various phase interpolation methods of rsSinusoidalSynthesizer for the given
  sinusoidal partial. We let a rsSinusoidalSynthesizer generate the interpolated phases as it would
  do in the actual synthesis (verify, if this is really the same algo) and plot the resulting 
  interpolated phases for the given sample rate. */
  static void plotInterpolatedPhases(const RAPT::rsSinusoidalPartial<T>& partial, T sampleRate);


  /** Analyzes the given sampleData of length numSamples with the given rsSinusoidalAnalyzer object
  and plots the analysis results. */
  void plotAnalysisResult(RAPT::rsSinusoidalAnalyzer<T>& sa, 
    T* sampleData, int numSamples, T sampleRate);
  // can this be made static?

  // maybe have a setSampleRate function

  // rename to plotSineTracks - factor out addSineTrackData function

protected:

  /** Returns an RGB color string that should be used for the partial with given index in the given 
  model. Partials are colored according to their strength. */
  std::string getPartialColor(const RAPT::rsSinusoidalModel<T>& model, size_t partialIndex);

  int graphIndex = 1; 
  // counter for the graphs that have been added to a plot, needed to correctly set up the colors

};

// todo: make a baseclass SpectrogramPlotter

//=================================================================================================

/** A class for plotting graphs in the vertices-and-edges sense. */

template<class T>
class GraphPlotter : public GNUPlotter
{

public:

  /** Plots a graph in which each vertex has data that determines its location in the 2D plane. You 
  can pass an array of vertex indices for vertices to be highlighted. These are drawn with bigger 
  dots and thicker edges and their nighbors are also bigger than usual but not quite as big as the 
  actually highlighted vertices. */
  void plotGraph2D(const rsGraph<rsVector2D<T>, T>& m, 
    std::vector<int> highlight = std::vector<int>());
  // make static


  /** Under construction */
  static void plotMeshFunction(const rsGraph<rsVector2D<T>, T>& m, const std::vector<T>& u);

};


#endif