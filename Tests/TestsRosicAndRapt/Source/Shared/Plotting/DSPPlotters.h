#ifndef DSPPLOTTERS_H
#define DSPPLOTTERS_H

#include "GNUPlotter.h"

/* Subclasses of GNUPlotter that specialize in making plots related to digital signal processing 
(DSP). */


/** A structure to specify a filter in terms of its zeros, poles, a gain factor and possibly a 
sample rate. If the sample rate is infinity (the default), an analog filter is assumed and the 
zeros and poles are in the s-plane. If it's nonzero, the filter is assumed to be digital and the
zeros and poles are in the z-plane. */

template <class T>
struct FilterSpecificationZPK
{
  std::vector<std::complex<T>> poles;
  std::vector<std::complex<T>> zeros;
  T gain = 1;
  T sampleRate = std::numeric_limits<T>::infinity();
};


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
  void addFilterSpecification(int numPoles, std::complex<T>* poles, int numZeros, 
    std::complex<T>* zeros,  T gain, T sampleRate = inf);

  /** Adds a filter specification via an object the FilterSpecificationZPK class. */
  void addFilterSpecification(const FilterSpecificationZPK<T>& spec);

  /** Plots the magnitude responses of all the filters. */
  void plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels);
  /*
  void plotPhase();
  void plotMagnitudeAndPhase(); // in one plot
  void plotPhaseDelay();
  void plotGroupDelay();
  void plotPolesAndZeros();
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
    std::vector<std::complex<T>>& p, T k);

protected:

  T freqScale = 1.0;

  std::vector<FilterSpecificationZPK<T>> filterSpecs; 

};

// todo: allow for filters to be specified also in terms of their polynomial coefficients - maybe 
// either by incoroprating a root finder or by allowing to store a filter specification in one of
// the two alternative formats (zeros-poles-gain (zpk) or coeff-arrays), conversion from zpk to
// coeffs may be required anyway to compute group delay (i think we'll need derivatives of 
// numerator and denominator, quotient-rule and maybe chain-rule to evaluate the derivative of the
// phase-response)

// todo: make a class SpectrogramPlotter

#endif