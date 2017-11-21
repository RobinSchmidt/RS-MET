#include "DSPPlotters.h"
using namespace std;

template <class T>
const T FilterPlotter<T>::inf = std::numeric_limits<T>::infinity();

template <class T>
const T FilterPlotter<T>::pi = T(3.1415926535897932384626433832795);

template <class T>
FilterPlotter<T>::FilterPlotter()
{

}

template <class T>
void FilterPlotter<T>::setFrequenciesAreRadian(bool areRadian)
{
  if(areRadian)
    freqScale = T(1);
  else
    freqScale = 2*pi;
}

template <class T>
void FilterPlotter<T>::addFilterSpecification(int numPoles, complex<T>* poles, int numZeros, 
  complex<T>* zeros, T gain, T sampleRate)
{
  FilterSpecificationZPK<T> spec;
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
  vector<T> f = getFrequencyAxis(numFreqs, lowFreq, highFreq, logFreqAxis);
  for(unsigned int i = 0; i < filterSpecs.size(); i++) {
    vector<complex<T>> H = getFrequencyResponse(i, f);
    vector<T> mag = getMagnitudes(H);
    addDataArrays(numFreqs, &f[0], &mag[0]);
    addGraph(string("i 0 u 1:") + s(i+2) + string(" w lines lw 1.5 axes x1y1 notitle"));
  }
  plot();
}

template <class T>
vector<T> FilterPlotter<T>::getFrequencyAxis(int numFreqs, T lowFreq, T highFreq, bool logarithmic)
{
  vector<T> freqs(numFreqs);
  if(logarithmic)
    rangeLogarithmic(&freqs[0], numFreqs, lowFreq, highFreq);
  else
    rangeLinear(&freqs[0], numFreqs, lowFreq, highFreq);
  return freqs;
}

template <class T>
vector<complex<T>> FilterPlotter<T>::getFrequencyResponse(int index, vector<T>& f)
{
  FilterSpecificationZPK<T> spec = filterSpecs[index];
  bool isDigital = spec.sampleRate != inf;
  complex<T> j(0.0, 1.0);                          // imaginary unit                         
  vector<complex<T>> H(f.size());                  // frequency response
  for(int k = 0; k < f.size(); k++) {
    complex<T> s = j * complex<T>(freqScale*f[k]); // value on s-plane where we evaluate H
    if(isDigital)
      s = exp(s/spec.sampleRate);                  // conversion of analog "s" to digital "z"
    H[k] = transferFunctionZPK(s, spec.zeros, spec.poles, spec.gain); 
  }
  return H;
}

template <class T>
vector<T> FilterPlotter<T>::getMagnitudes(vector<complex<T>>& H)
{
  vector<T> mag(H.size());
  for(int k = 0; k < H.size(); k++)
    mag[k] = abs(H[k]);
  return mag;
}

template <class T>
complex<T> FilterPlotter<T>::polynomialByRoots(complex<T> z, vector<complex<T>>& r)
{
  complex<T> w = 1;
  for(int i = 0; i < r.size(); i++)
    w *= z-r[i];
  return w;
}

template <class T>
complex<T> FilterPlotter<T>::transferFunctionZPK(complex<T> s, vector<complex<T>>& z, 
  vector<complex<T>>& p, T k)
{
  complex<T> num = polynomialByRoots(s, z);
  complex<T> den = polynomialByRoots(s, p);
  return k * num/den;
}

// template instantiations:
template FilterPlotter<float>;
template FilterPlotter<double>;