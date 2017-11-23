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
void FilterPlotter<T>::addFilterSpecification(const FilterSpecificationZPK<T>& spec)
{
  filterSpecs.push_back(spec);
}

template <class T>
void FilterPlotter<T>::plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, 
  bool decibels)
{
  vector<vector<vector<T>>> data(1);
  data[0].resize(1 + filterSpecs.size());
  vector<T> f = getFrequencyAxis(numFreqs, lowFreq, highFreq, logFreqAxis);
  data[0][0] = f;
  for(unsigned int i = 0; i < filterSpecs.size(); i++) {
    vector<complex<T>> H = getFrequencyResponse(i, f);
    vector<T> mag = getMagnitudes(H);
    data[0][i+1] = mag;
    addGraph(string("i 0 u 1:") + s(i+2) + string(" w lines lw 1.5 axes x1y1 notitle"));
  }
  addDataBlockColumnLine(data);
  plot();
}

template <class T>
void FilterPlotter<T>::plotPolesAndZeros(int plotSize)
{
  for(unsigned int i = 0; i < filterSpecs.size(); i++)
  {
    addDataComplex(filterSpecs[i].poles);
    addDataComplex(filterSpecs[i].zeros);
    addGraph("i " + s(2*i)   + " u 1:2 w points pt 2 ps 1 notitle");
    addGraph("i " + s(2*i+1) + " u 1:2 w points pt 6 ps 1 notitle");

    // show the multiplicities of poles and zeros:
    T thresh = T(1.e-8);    // threshold for considering close poles/zeros as multiple root
                            // maybe use something that depends on the epsilon of T
    drawMultiplicities(filterSpecs[i].poles, thresh);
    drawMultiplicities(filterSpecs[i].zeros, thresh);
  }

  // todo: make the colors of the poles, zeros and multiplicities for each filter equal

  setupForPoleZeroPlot(plotSize);
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

template <class T>
void FilterPlotter<T>::setupForPoleZeroPlot(int size)
{
  bool zDomain = false;
  double range = 0;
  for(unsigned int i = 0; i < filterSpecs.size(); i++) 
  {
    range   = fmax(range, maxAbsReIm(filterSpecs[i].poles));  
    range   = fmax(range, maxAbsReIm(filterSpecs[i].zeros)); 
    zDomain = zDomain || (filterSpecs[i].sampleRate != inf);
  }
  range = 1.1 * fmax(1.0, range);
  setRange(-range, range, -range, range);

  addCommand("set size square");                // set aspect ratio to 1:1
  addCommand("set xzeroaxis lt 1");             // draw x-axis
  addCommand("set yzeroaxis lt 1");             // draw y-axis
  addCommand("set xlabel \"Real Part\"");
  addCommand("set ylabel \"Imaginary Part\"");
  setPixelSize(size, size);

  if(zDomain == true)
    addCommand("set object 1 ellipse at first 0,0 size 2,2 fs empty border rgb \"#808080\""); 
}

template <class T>
void FilterPlotter<T>::drawMultiplicities(const vector<complex<T>>& z, T thresh)
{
  size_t N = z.size();       // number of values
  vector<complex<T>> zd(N);  // collected distinct values
  vector<int> m(N);          // m[i] = multiplicity of value zd[i]
  vector<bool> done(N);      // vector of flags, if z[i] was already absorbed into zd
  int i, j;
  int k = 0;

  // collect distinct values and their multiplicities:
  for(i = 0; i < N; i++) {
    if(!done[i]) {
      zd[k]   = z[i];
      m[k]    = 1;
      done[i] = true;
      for(j = i+1; j < N; j++) { // find values equal to zd[k] == z[i]
        if(!done[j] && almostEqual(z[i], z[j], thresh)) {
          m[k]    += 1;
          done[j]  = true; }}
      k++; }}

  // k is now the number of distinct values stored in zd with associated multiplicities in m
  for(i = 0; i < k; i++){
    if(m[i] > 1)
      addAnnotation(zd[i].real(), zd[i].imag(), " " + to_string(m[i]), "left"); }
}

template <class T>
double FilterPlotter<T>::maxAbsReIm(const vector<complex<T>>& x)
{
  double m = 0.0;
  for(int i = 0; i < x.size(); i++)
  {
    if(fabs(x[i].real()) > m)
      m = fabs(x[i].real());
    if(fabs(x[i].imag()) > m)
      m = fabs(x[i].imag());
  }
  return m;
}

template <class T>
bool  FilterPlotter<T>::almostEqual(complex<T> x, complex<T> y, T thresh)
{
  return abs(x-y) / fmax(abs(x), abs(y)) < thresh;
}

// template instantiations:
template FilterPlotter<float>;
template FilterPlotter<double>;