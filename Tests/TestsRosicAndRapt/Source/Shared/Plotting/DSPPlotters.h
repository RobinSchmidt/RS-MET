#ifndef DSPPLOTTERS_H
#define DSPPLOTTERS_H

#include "GNUPlotter.h"

/* Subclasses of GNUPlotter that specialize in making plots related to digital signal processing 
(DSP). */

/** A class for visualizing analog and digital filter responses. */

template <class T>
class FilterPlotter : public GNUPlotter
{

public:

  FilterPlotter(bool isDigital);

  void addPoleZeroSet(int numPoles, std::complex<T>* poles, int numZeros, std::complex<T>* zeros,
    T gain);

  void plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels);
  /*
  void plotPhase();
  void plotMagnitudeAndPhase(); // in one plot
  void plotPolesAndZeros();
  void plotTransferFunctionMagnitude();
  void plotPhaseDelay();
  void plotGroupDelay();
  void plotImpulseResponse();
  void plotStepResponse();
  */

  /** Creates a vector of x-values for the frequency axis either linearly or logarithmically 
  scaled. */
  std::vector<T> getFrequencyAxis(int numFreqs, T lowFreq, T highFreq, bool logarithmic);

  /** Returns the complex frequency response of the filter with given index at the frequencies 
  given in the vector. */
  std::vector<std::complex<T>> getFrequencyResponse(int index, std::vector<T> frequencies);

protected:

  bool isDigital = true;

  struct FilterSpecification
  {
    std::vector<T> poles;
    std::vector<T> zeros;
    T gain;
  };
  std::vector<FilterSpecification> filterSpecs; 

};

// todo: make a class SpectrogramPlotter

#endif