#ifndef DSPPLOTTERS_H
#define DSPPLOTTERS_H

#include "GNUPlotter.h"

/* Subclasses of GNUPlotter that specialize in making plots related to digital signal processing 
(DSP). */

/** A class for visualizing analog and digital filter responses. It may plot magnitude responses,
phase responses, pole/zero plots etc. for one or more filters. You need to specify the filters in 
terms of their poles and zeros. For each filter, you once call addPoleZeroSet to add the filter to
the list. Once you are finished adding filters this way, you can get the various plots via the
respective plot... functions. */

template <class T>
class FilterPlotter : public GNUPlotter
{

public:

  static const T inf = std::numeric_limits<T>::infinity(); // shorthand

  FilterPlotter();

  /** Adds a filter specification in terms of poles, zeros and gain to our list. You may also pass 
  a sampleRate in which case the poles and zeros will be interpreted as z-plane values. Otherwise
  an analog filter (corresponding to an infinite sampl rate) will be assumed and the poles and 
  zeros are interpreted as being in the s-plane */
  void addPoleZeroSet(int numPoles, std::complex<T>* poles, int numZeros, std::complex<T>* zeros,
    T gain, T sampleRate = inf);
  // maybe rename to addFilter

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

  struct FilterSpecification
  {
    std::vector<T> poles;
    std::vector<T> zeros;
    T gain = 1;
    T sampleRate = inf;
  };
  std::vector<FilterSpecification> filterSpecs; 

};

// todo: make a class SpectrogramPlotter

#endif