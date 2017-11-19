#include "DSPPlotters.h"
using namespace std;

template <class T>
FilterPlotter<T>::FilterPlotter()
{

}

template <class T>
void FilterPlotter<T>::addPoleZeroSet(int numPoles, complex<T>* poles, int numZeros, 
  complex<T>* zeros, T gain, T sampleRate)
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
  spec.sampleRate = sampleRate;
  filterSpecs.push_back(spec);
}

template <class T>
void FilterPlotter<T>::plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, 
  bool decibels)
{

}

template <class T>
vector<T> FilterPlotter<T>::getFrequencyAxis(int numFreqs, T lowFreq, T highFreq, 
  bool logarithmic)
{
  vector<T> freqs(numFreqs);
  if(logarithmic)
    rangeLogarithmic(&freqs[0], numFreqs, lowFreq, highFreq);
  else
    rangeLinear(&freqs[0], numFreqs, lowFreq, highFreq);
  return freqs;
}

template <class T>
vector<complex<T>> FilterPlotter<T>::getFrequencyResponse(int index, vector<T> f)
{
  FilterSpecification spec = filterSpecs[index];
  bool isDigital = spec.sampleRate == inf;
  complex<T> j(0.0, 1.0);                // imaginary unit
  complex<T> s;                          // value on s-plane where w evaluate H
  vector<complex<T>> H(f.size());        // frequency response
  for(int k = 0; k < f.size(); k++) {
    s = j * 2 * M_PI * f[k];
    if(isDigital)
      s = exp(s/sampleRate); // more typically called "z"
    H[k] = transferFunctionZPK(s, spec.zeros, spec.poles, spec.gain); 
  }
  return H;
}

template <class T>
complex<T> FilterPlotter<T>::polynomialByRoots(complex<T> z, vector<complex<T>> r)
{
  complex<T> w = 1;
  for(int i = 0; i < r.size(); i++)
    w *= z-r[i];
  return w;
}

template <class T>
complex<T> FilterPlotter<T>::transferFunctionZPK(complex<T> s, vector<complex<T>> z, 
  vector<complex<T>> p, T k)
{
  complex<T> num = polynomialByRoots(s, z);
  complex<T> den = polynomialByRoots(s, p);
  return k * num/den;
}
