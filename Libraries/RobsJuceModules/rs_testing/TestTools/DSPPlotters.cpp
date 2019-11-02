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
void FilterPlotter<T>::addFilterSpecificationZPK(int numPoles, complex<T>* poles, int numZeros,
  complex<T>* zeros, T gain, T sampleRate)
{
  rsFilterSpecificationZPK<T> spec;
  spec.p.resize(numPoles);
  spec.z.resize(numZeros);
  int i;
  for(i = 0; i < numPoles; i++)
    spec.p[i] = poles[i];
  for(i = 0; i < numZeros; i++)
    spec.z[i] = zeros[i];
  spec.k = gain;
  spec.sampleRate = sampleRate;
  filterSpecsZPK.push_back(spec);
  // todo: factor out into constructor of rsFilterSpecificationZPK that takes raw arrays
}

template <class T>
void FilterPlotter<T>::addFilterSpecificationZPK(const RAPT::rsFilterSpecificationZPK<T>& spec)
{
  filterSpecsZPK.push_back(spec);
  filterSpecsBA.push_back(spec.toBA());
}

template <class T>
void FilterPlotter<T>::addFilterSpecificationBA(int numeratorOrder, T* numeratorCoeffs,
  int denominatorOrder, T* denominatorCoeffs, T sampleRate)
{
  rsFilterSpecificationBA<T> spec;
  spec.b.resize(numeratorOrder+1);
  spec.a.resize(denominatorOrder+1);
  spec.sampleRate = sampleRate;
  for(size_t i = 0; i <= numeratorOrder; i++)
    spec.b[i] = numeratorCoeffs[i];
  for(size_t i = 0; i <= denominatorOrder; i++)
    spec.a[i] = denominatorCoeffs[i];
  addFilterSpecificationBA(spec);
}

template <class T>
void FilterPlotter<T>::addFilterSpecificationBA(const RAPT::rsFilterSpecificationBA<T>& spec)
{
  filterSpecsBA.push_back(spec);
  filterSpecsZPK.push_back(spec.toZPK());
}

template <class T>
void FilterPlotter<T>::plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis,
  bool decibels)
{
  vector<vector<vector<T>>> data(1);
  data[0].resize(1 + filterSpecsZPK.size());
  vector<T> f = getFrequencyAxis(numFreqs, lowFreq, highFreq, logFreqAxis);
  data[0][0] = f;


  for(unsigned int i = 0; i < filterSpecsZPK.size(); i++) {
    vector<complex<T>> H = getFrequencyResponse(i, f);
    vector<T> mag = getMagnitudes(H);

    if(decibels)
      for(size_t k = 0; k < mag.size(); k++) // factor out into function toDecibels
        mag[k] = 20*log10(mag[k]);

    data[0][i+1] = mag;  // refactor to data[0][i+1] = getMagnitudeResponse(i, f);
    addGraph(string("i 0 u 1:") + s(i+2) + string(" w lines lw 1.5 axes x1y1 notitle"));
  }

  // getFrequencyResponse(i, f); should dispatch between zpk and ba, loop should run over
  // i = 0...getNumFilterSpecs() which returns the sum of the zpk and ba specs

  addDataBlockColumnLine(data);
  if(logFreqAxis)
    setLogScale("x", 10); // 10 is the base - maybe try 2
  plot();
}




template <class T>
void FilterPlotter<T>::plotPolesAndZeros(int plotSize)
{
  for(unsigned int i = 0; i < filterSpecsZPK.size(); i++)
  {
    addDataComplex(filterSpecsZPK[i].p);
    addDataComplex(filterSpecsZPK[i].z);
    addGraph("i " + s(2*i)   + " u 1:2 w points pt 2 ps 1 notitle");
    addGraph("i " + s(2*i+1) + " u 1:2 w points pt 6 ps 1 notitle");

    // show the multiplicities of poles and zeros:
    T thresh = T(1.e-8);    // threshold for considering close poles/zeros as multiple root
                            // maybe use something that depends on the epsilon of T
    drawMultiplicities(filterSpecsZPK[i].p, thresh);
    drawMultiplicities(filterSpecsZPK[i].z, thresh);
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
  rsFilterSpecificationZPK<T> spec = filterSpecsZPK[index];
  //bool isDigital = spec.sampleRate != inf;
  complex<T> j(0.0, 1.0);                          // imaginary unit
  vector<complex<T>> H(f.size());                  // frequency response
  for(int k = 0; k < f.size(); k++) {
    complex<T> s = j * complex<T>(freqScale*f[k]); // value on s-plane where we evaluate H
    if(spec.isDigital())
      s = exp(s/spec.sampleRate);                  // conversion of analog "s" to digital "z"

    // debug:
    //complex<T> H1 = filterSpecsZPK[index].transferFunctionAt(s);
    //complex<T> H2 = filterSpecsBA[index].transferFunctionAt(s);


    H[k] = spec.transferFunctionAt(s);
    //H[k] = transferFunctionZPK(s, spec.z, spec.p, spec.k);
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
  vector<complex<T>>& p, std::complex<T> k)
{
  complex<T> num = polynomialByRoots(s, z);
  complex<T> den = polynomialByRoots(s, p);
  return k * num/den;
}

// i think, this can be deleted (was moved to RAPT - check this and if so, delete):
/*
template <class T>
FilterSpecificationBA<T> FilterPlotter<T>::zpk2ba(const FilterSpecificationZPK<T>& zpk)
{
  FilterSpecificationBA<T> ba;
  ba.sampleRate = zpk.sampleRate;
  ba.a = rsPolynomial<T>::getPolynomialCoefficientsFromRoots(zpk.poles); // rename: rootsToCoeffs
  ba.b = rsPolynomial<T>::getPolynomialCoefficientsFromRoots(zpk.zeros); // have also coeffsToRoots

  // todo: reverse a,b in case of digital filter

  for(size_t i = 0; i < ba.b.size(); i++)
    ba.b[i] *= zpk.gain;
  normalizeA0(ba); // maybe make the normalization optional (on by default)
  return ba;
} // moved to rsFilterSpecificationZPK<T>::toBA
*/

/*
template <class T>
FilterSpecificationZPK<T> FilterPlotter<T>::ba2zpk(const FilterSpecificationBA<T>& ba)
{
  FilterSpecificationZPK<T> zpk;
  zpk.sampleRate = ba.sampleRate;
  zpk.poles.resize(ba.a.size()-1);
  zpk.zeros.resize(ba.b.size()-1);

  // todo: reverse a,b in case of digital filter

  rsPolynomial<T>::findPolynomialRoots(&ba.a[0], (int) ba.a.size()-1, &zpk.poles[0]);
  rsPolynomial<T>::findPolynomialRoots(&ba.b[0], (int) ba.b.size()-1, &zpk.zeros[0]);

  //if(ba.sampleRate != inf) // digital
  //  //zpk.gain = ba.b[0];   // maybe, it should not be just b0 but b0/a0 -> verify/test...
  //  zpk.gain = ba.b[0] / ba.a[0];  // gain is quotient of leading coeffs
  //else
  //  zpk.gain = ba.b[ba.b.size()-1] / ba.a[ba.a.size()-1] ;
  // is this correct?  ...verify -  has been tested only in the digital case

  //zpk.gain = ba.b[0]/ba.a[0];

  zpk.gain = ba.b[ba.b.size()-1] / ba.a[ba.a.size()-1];


  //zpk.gain = 1; // preliminary
  return zpk;
}
*/

/*
template <class T>
void FilterPlotter<T>::normalizeA0(FilterSpecificationBA<T>& ba)
{
  std::complex<T> s = T(1) / ba.a[0];
  for(size_t i = 0; i < ba.a.size(); i++) ba.a[i] *= s;
  for(size_t i = 0; i < ba.b.size(); i++) ba.b[i] *= s;
} // used inside zpk2ba
*/

template <class T>
void FilterPlotter<T>::setupForPoleZeroPlot(int size)
{
  //bool zDomain = true;
  bool zDomain = false;
  double range = 0;
  for(unsigned int i = 0; i < filterSpecsZPK.size(); i++)
  {
    range   = fmax(range, maxAbsReIm(filterSpecsZPK[i].p));
    range   = fmax(range, maxAbsReIm(filterSpecsZPK[i].z));
    //zDomain = zDomain || (filterSpecsZPK[i].sampleRate != inf);
    zDomain = zDomain || filterSpecsZPK[i].isDigital();
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
template class FilterPlotter<float>;
template class FilterPlotter<double>;

//=================================================================================================


template <class T>
void SpectrumPlotter<T>::plotDecibelSpectra(int signalLength, const T *x0, const T *x1,
  const T *x2, const T *x3, const T *x4, const T *x5, const T *x6, const T *x7, const T *x8,
  const T *x9)
{
  RAPT::rsAssert(signalLength <= fftSize);

  // maybe factor out into setupTransformer:
  typedef RAPT::rsFourierTransformerRadix2<T> FT;
  transformer.setNormalizationMode(FT::NORMALIZE_ON_FORWARD_TRAFO);
  transformer.setDirection(        FT::FORWARD);
  transformer.setBlockSize(fftSize);

  // use this for y-axis minimum - let the user set it up:
  T ampFloor = RAPT::rsDbToAmp(dBFloor);

  //std::vector<T*> inputArrays = collectLeadingNonNullArguments(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9);
  const vector<const T*> inputArrays = collectLeadingNonNullArguments(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9);

  int N = rsMax(signalLength, fftSize);
  int maxBin = fftSize/2; // later have a user option for that -> zoom

  std::vector<T> f = getFreqAxis(maxBin);

  std::vector<std::complex<T>> tmp(N);
  std::vector<T> dB(N);
  //std::vector<T> phs(N);
  for(size_t i = 0; i < inputArrays.size(); i++) {
    RAPT::rsArray::convertBuffer(inputArrays[i], &tmp[0], signalLength);
    RAPT::rsArray::fillWithZeros(&tmp[signalLength], N-signalLength);
    transformer.transformComplexBufferInPlace(&tmp[0]);

    // this may be not quite correct at DC (i think, because we need to incorporate the value
    // at fftSize/2 or something?)
    T compFactor = T(fftSize) / T(signalLength);
    for(size_t k = 0; k < N; k++)
      dB[k] = RAPT::rsAmpToDbWithCheck(compFactor * abs(tmp[k]), ampFloor);

    addDataArrays(maxBin, &f[0], &dB[0]);
    //addDataArrays(fftSize/2, &dB[0]); // maybe fftSize/2 or (fftSize+1)/2
    int dummy = 0;
  }

  //setLogScale("x"); // make this a user option

  plot();
}
// maybe factor out addSpectra

template <class T>
std::vector<T> SpectrumPlotter<T>::getFreqAxis(int maxBin)
{
  std::vector<T> f(maxBin);

  // todo: check everything for off-by-one errors

  GNUPlotter::rangeLinear(&f[0], maxBin+1, T(0), T(maxBin));
  T scaler = T(1);
  T r = T(1) / T(fftSize);
  typedef FreqAxisUnits FU;
  switch(freqAxisUnit)
  {
  case FU::normalized: scaler = r*T(1);       break;
  case FU::omega:      scaler = r*T(2*PI);    break;
  case FU::hertz:      scaler = r*sampleRate; break;
  }

  RAPT::rsArray::scale(&f[0], maxBin, scaler);

  return f;
}

// template instantiations:
template class SpectrumPlotter<float>;
template class SpectrumPlotter<double>;
// move to somewhere else - a template instantiations file

//=================================================================================================

template <class T>
void SpectrogramPlotter<T>::addSpectrogramData(GNUPlotter& p, int numFrames, int numBins,
  T **s, T fs, int H, T dbMin, T dbMax)
{
  // fs: sample rate, H: hop size

  // create time- and frequency axis and add the data:
  double tMax = H*(numFrames-1)/fs;
  double fMax =  0.5*fs*(numBins-1)/numBins;
  double *t = new double[numFrames];
  double *f = new double[numBins];
  p.rangeLinear(t, numFrames, 0.0, tMax);
  p.rangeLinear(f, numBins,   0.0, fMax);


  p.setDataPrecision(4);                                  // make temp files smaller
  p.addDataMatrix(numFrames, numBins, t, f, s);

  // set up style settings:
  //p.setAxisLabels("Time in s", "Frequency in Hz", "Level in dB");
  p.setPixelSize(800, 400);
  p.addGraph("i 0 nonuniform matrix w image notitle");

  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  p.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale


  // maybe use signal-dependent values later:
  //double dbMin = -100.0;
  //double dbMax =   10.0;
  p.setRange(0.0, tMax, 0.0, fMax, dbMin, dbMax);
  // doesn't seem to have any effect on the color axis (z). it seems like GNUPlot uses
  // autoscaling for it, no matter what.

  delete[] t;
  delete[] f;
}




//=================================================================================================


template <class T>
void SinusoidalModelPlotter<T>::addModelToPlot(
  const RAPT::rsSinusoidalModel<T>& model, GNUPlotter& plt,T sampleRate,
  const std::string& graphColor)
{
  std::vector<float> t, f; // we plot frequency vs time
  for(size_t i = 0; i < model.getNumPartials(); i++) {
    rsSinusoidalPartial<double> p = model.getPartial(i);
    size_t L = p.getNumDataPoints();
    t.resize(L);
    f.resize(L);
    for(size_t j = 0; j < L; j++) {
      t[j] = (float)p.getDataPoint(j).time;
      f[j] = (float)p.getDataPoint(j).freq;
    }
    plt.addDataArrays((int)L, &t[0], &f[0]);

    if(graphColor != "")
      plt.setGraphColor(graphIndex, graphColor);
    else
      plt.setGraphColor(graphIndex, getPartialColor(model, i));

    plt.setGraphStyles("lines lw 1.5");

    graphIndex++;
  }
  // maybe factor out an addPartialToPlot function - may become useful when we want to plot
  // partials seperately
}

template <class T>
void SinusoidalModelPlotter<T>::plotModel(const RAPT::rsSinusoidalModel<T>& model, T fs)
{
  GNUPlotter plt;
  addModelToPlot(model, plt, fs, "000000");
  plt.plot();
}

template <class T>
void SinusoidalModelPlotter<T>::plotTwoModels(
  const RAPT::rsSinusoidalModel<T>& model1, const RAPT::rsSinusoidalModel<T>& model2, T fs)
{
  GNUPlotter plt;
  addModelToPlot(model1, plt, fs, "000000");  // black
  addModelToPlot(model2, plt, fs, "0000FF");  // blue
  plt.plot();
}

template <class T>
std::string SinusoidalModelPlotter<T>::getPartialColor(
  const RAPT::rsSinusoidalModel<T>& mdl, size_t i)
{
  return "000000"; // preliminary
}


// this is a kludge - get rid
template<class T>
T** createRowPointers(RAPT::rsMatrix<T>& M)
{
  int N = M.getNumRows();
  T** rp = new T*[N];
  for(int i = 0; i < N; i++)
    rp[i] = M.getRowPointer(i);
  return rp;
}
template<class T>
void deleteRowPointers(T** rowPointers, const RAPT::rsMatrix<T>& M)
{
  for(int i = 0; i < M.getNumRows(); i++)
    delete rowPointers[i];
}


template <class T>
void SinusoidalModelPlotter<T>::plotAnalysisResult(
  RAPT::rsSinusoidalAnalyzer<T>& sa, T* x, int N, T fs)
{
  GNUPlotter plt;

  // factor out an addSpectrogramToPlot or something
  // todo: plot the spectrogram underneath the sine tracks - when doing so, the tracks are not
  // plotted
  bool plotSpectrogram = false; // make member, let user set it
  if(plotSpectrogram == true)
  {
    RAPT::rsMatrix<std::complex<T>> stft = sa.getComplexSpectrogram(x, N);

    int numBins   = stft.getNumRows();
    int numFrames = stft.getNumColumns();

    RAPT::rsMatrix<T> dB(stft.getNumRows(), stft.getNumColumns());
    for(int i = 0; i < dB.getNumRows(); i++)
      for(int j = 0; j < dB.getNumColumns(); j++)
        dB(i, j) = rsAmpToDb(abs(stft(i, j)));

    dB.transpose();


    T** rowPointers = createRowPointers(dB);
    this->addSpectrogramData(
      plt, numFrames, numBins, rowPointers, fs, sa.getHopSize(), this->minDb, this->maxDb);
    deleteRowPointers(rowPointers, dB);


    // old:
    //this->addSpectrogramData(
    //  plt, numFrames, numBins, dB.getDataPointerConst(),
    //  fs, sa.getHopSize(), this->minDb, this->maxDb);


  }

  // create model and add it to the plot:
  RAPT::rsSinusoidalModel<double> model = sa.analyze(x, N, fs);
  addModelToPlot(model, plt, fs, "000000");

  plt.plot();
}




template class SinusoidalModelPlotter<double>;
