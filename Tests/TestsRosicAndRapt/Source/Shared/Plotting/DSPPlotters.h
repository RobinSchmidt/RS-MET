#ifndef DSPPLOTTERS_H
#define DSPPLOTTERS_H

#include "GNUPlotter.h"
//#include "../../../../Libraries/JUCE/modules/rosic/rosic.h"
#include "../../../Libraries/JUCE/modules/rosic/rosic.h"


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
  digital filters. By default, i.e. if you don't call this function, they are interpreted as 
  radian frequencies. */
  void setFrequenciesAreRadian(bool areRadian);

  /** Adds a filter specification in terms of poles, zeros and gain to our list. You may also pass 
  a sampleRate in which case the poles and zeros will be interpreted as z-plane values. Otherwise
  an analog filter (corresponding to an infinite sample rate) will be assumed and the poles and 
  zeros are interpreted as being in the s-plane */
  void addFilterSpecificationZPK(int numPoles, std::complex<T>* poles, int numZeros, 
    std::complex<T>* zeros,  T gain, T sampleRate = inf);

  /** Adds a filter specification via an object the rsFilterSpecificationZPK class. */
  void addFilterSpecificationZPK(const RAPT::rsFilterSpecificationZPK<T>& spec);


  void addFilterSpecificationBA(int numeratorOrder, T* numeratorCoeffs, 
    int denominatorOrder, T* denominatorCoeffs, T sampleRate = inf);
  // maybe allow for complex coeffs

  void addFilterSpecificationBA(const RAPT::rsFilterSpecificationBA<T>& spec);



  /** Plots the magnitude responses of all the filters. */
  void plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels);

  /** Plots the poles and zeros of all the filters in s- or z-plane. You can pass the plot-size in 
  pixels (it will always use a square plot). */
  void plotPolesAndZeros(int plotSize = 400);

  /*
  void plotPhase();
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
  std::vector<std::complex<T>> getFrequencyResponse(int index, std::vector<T>& frequencies);

  /** Extracts the magnitudes from the passed complex frequency response array.  */
  std::vector<T> getMagnitudes(std::vector<std::complex<T>>& complexFreqResponse);

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

  /** Returns maximum absolute value of all real an imaginary parts. */
  double maxAbsReIm(const std::vector<std::complex<T>>& x);

  /** Returns true, if the relative distance between x and y is smaller than the given threshold 
  ("relative" with respect to the actual absolute values of x and y, such that for larger values 
  the tolerance also increases) */
  bool almostEqual(std::complex<T> x, std::complex<T> y, T thresh);

  T freqScale = 1.0;
  std::vector<RAPT::rsFilterSpecificationZPK<T>> filterSpecsZPK;
  std::vector<RAPT::rsFilterSpecificationBA<T>>  filterSpecsBA;

  // have conversion functions convert_BA_To_ZPK, convert_ZPK_To_BA
  // or filterSpecBA2ZPK, ba2zpk, zpk2ba
  // maybe keep tpk and ba specifications for each filter, i.e. keep the two representaions in sync
  // and use whatever representation is more convenient to deal with to compute the various
  // plots/numbers

  // maybe we could also keep SOS representations and maybe others?

};

// todo: allow for filters to be specified also in terms of their polynomial coefficients - maybe 
// either by incoroprating a root finder or by allowing to store a filter specification in one of
// the two alternative formats (zeros-poles-gain (zpk) or coeff-arrays), conversion from zpk to
// coeffs may be required anyway to compute group delay (i think we'll need derivatives of 
// numerator and denominator, quotient-rule and maybe chain-rule to evaluate the derivative of the
// phase-response)


template <class T>
class SpectrumPlotter : public GNUPlotter
{

public:

  /** Given up to 10 signal buffers of length "signalLength", this function performs an FFT on each 
  of them and plot the spectral magintudes as decibel values.  
  (the FFT size is determined by setFfftSize and may be different from signalLength) */
  //template <class T>
  void plotDecibelSpectra(int signalLength, T *x0, T *x1 = nullptr, T *x2 = nullptr, 
    T *x3 = nullptr, T *x4 = nullptr, T *x5 = nullptr, T *x6 = nullptr, T *x7 = nullptr, 
    T *x8 = nullptr, T *x9 = nullptr);


  void setFftSize(int newSize) { fftSize = newSize; }

  // setSampleRate, etc...

protected:


  int fftSize = 2048;

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

#include "../Prototypes/SinusoidalModeling.h" // should not be necessary - figure out why this is needed

template<class T>
class SinusoidalModelPlotter : public SpectrogramPlotter<T>
{

public:

  //void plot(SinusoidalAnalyzer<T>& sa, T* sampleData, int N, T sampleRate);

  void plotAnalysisResult(SinusoidalAnalyzer<T>& sa, T* sampleData, int N, T sampleRate);

  // rename to plotSineTracks - factor out addSineTrackData function

protected:

};

// todo: make a baseclass SpectrogramPlotter


#endif