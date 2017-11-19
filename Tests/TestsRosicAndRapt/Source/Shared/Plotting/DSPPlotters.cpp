#include "DSPPlotters.h"
using namespace std;

template <class T>
FilterPlotter<T>::FilterPlotter(bool isDigital)
{
  this->isDigital = isDigital;
}

template <class T>
void FilterPlotter<T>::addPoleZeroSet(int numPoles, std::complex<T>* poles, int numZeros, 
  std::complex<T>* zeros, T gain)
{
  FilterSpecification spec;
  spec.poles.resize(numPoles);
  spec.zeros.resize(numZeros);
  int i;
  for(i = 0; i < numPoles; i++)
    spec.poles[i] = poles[i];
  for(i = 0; i < numZeros; i++)
    spec.zeros[i] = zeros[i];
  spec.gain = gain;
  filterSpecs.push_back(spec);
}

template <class T>
void FilterPlotter<T>::plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, 
  bool decibels)
{

}

template <class T>
std::vector<T> FilterPlotter<T>::getFrequencyAxis(int numFreqs, T lowFreq, T highFreq, 
  bool logarithmic)
{
  std::vector<T> freqs(numFreqs);
  if(logarithmic)
    rangeLogarithmic(&freqs[0], numFreqs, lowFreq, highFreq);
  else
    rangeLinear(&freqs[0], numFreqs, lowFreq, highFreq);
  return freqs;
}
