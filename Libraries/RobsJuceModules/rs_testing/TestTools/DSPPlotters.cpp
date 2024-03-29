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
void FilterPlotter<T>::addFilterSpecificationBA(int numeratorOrder, const T* numeratorCoeffs,
  int denominatorOrder, const T* denominatorCoeffs, T sampleRate)
{
  rsFilterSpecificationBA<T> spec;
  spec.b.resize(numeratorOrder+1);
  spec.a.resize(denominatorOrder+1);
  spec.sampleRate = sampleRate;
  for(int i = 0; i <= numeratorOrder; i++)
    spec.b[i] = numeratorCoeffs[i];
  for(int i = 0; i <= denominatorOrder; i++)
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
void FilterPlotter<T>::addTransferFunction(const RAPT::rsRationalFunction<T>& tf, T fs)
{
  RAPT::rsFilterSpecificationBA<T> ba;
  ba.sampleRate = fs;
  rsConvert(tf.getNumerator().getCoeffs(),   ba.b);
  rsConvert(tf.getDenominator().getCoeffs(), ba.a);
  if(fs != inf) { // z-domain transfer functions are in terms of z^-1
    rsArrayTools::reverse(&ba.b[0], (int)ba.b.size());
    rsArrayTools::reverse(&ba.a[0], (int)ba.a.size()); }
  addFilterSpecificationBA(ba);
}

template <class T>
void FilterPlotter<T>::plotMagnitude(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis,
  bool decibels)
{
  vector<vector<vector<T>>> data(1);           // 1 dataset
  data[0].resize(1 + filterSpecsZPK.size());   // 1 col for freq axis and N magnitude responses
  vector<T> f = getFrequencyAxis(numFreqs, lowFreq, highFreq, logFreqAxis);
  data[0][0] = f;

  for(unsigned int i = 0; i < filterSpecsZPK.size(); i++) {
    vector<complex<T>> H = getFrequencyResponse(i, f);
    data[0][i+1] = getMagnitudes(H, decibels);
    addGraphLines(i);
  }

  // getFrequencyResponse(i, f); should dispatch between zpk and ba, loop should run over
  // i = 0...getNumFilterSpecs() which returns the sum of the zpk and ba specs

  addDataBlockColumnLine(data);
  if(logFreqAxis)
    setLogScale("x", 10); // 10 is the base - maybe try 2
  plot();
}

template <class T>
void FilterPlotter<T>::plotFrequencyResponses(int numFreqs, T lowFreq, T highFreq, 
  bool logFreqAxis, bool plotMagnitude, bool decibels, bool plotPhase, bool unwrap)
{
  //using uint = unsigned int;

  int numFilters = (int)filterSpecsZPK.size();  // number of filters
  int numColumns = 1;                           // 1 for the freq axis values
  if(plotMagnitude) numColumns += numFilters;
  if(plotPhase)     numColumns += numFilters;

  // Create and add data:
  vector<vector<vector<T>>> data(1);           // 1 dataset
  data[0].resize(numColumns);                  // 
  vector<T> f = getFrequencyAxis(numFreqs, lowFreq, highFreq, logFreqAxis);
  data[0][0] = T(0.001) * f;                   // factor 0.001 because of kHz
  int j = 0;
  for(int i = 0; i < numFilters; i++) {
    vector<complex<T>> H = getFrequencyResponse(i, f);
    if(plotMagnitude) {
      data[0][j+1] =  getMagnitudes(H, decibels);
      addGraphLines(j);
      j++;  }
    if(plotPhase) {
      data[0][j+1] = getPhases(H, unwrap, true);
      addGraphLines(j, true);
      j++; }}
  addDataBlockColumnLine(data);

  // Set up plotting options and plot:
  //setTitle("Filter Frequency Responses");
  if(logFreqAxis)
    setLogScale("x", 10);
  addCommand("set xrange  [0.03125:32]");   // preliminary
  addCommand("set yrange  [-90:10]");       // preliminary
  //addCommand("set y2range [-450:0]");       // preliminary
  addCommand("set y2range [-405:45]");       // preliminary
  addCommand("set xlabel \"Frequency in kHz\"");
  addCommand("set ylabel \"Magnitude in dB\"");
  addCommand("set y2label \"Phase in Degrees\"");
  addCommand("set xtics 2");    // factor 2 between (major) frequency axis tics
  addCommand("unset mxtics");   // no minor tics for frequency axis
  addCommand("set ytics 10");   // 10 dB steps for magnitude axis
  addCommand("set y2tics 45");  // 45� steps for phase axis
  // ...actually, this should perhaps be left to client code - and then maybe we should have a 
  // convenience function that does this stuff automatically

  //drawNyquistLines(); // not yet implemented

  // Factor out:
  std::vector<std::string> hexColors = getGraphColors(numFilters);
  if(plotMagnitude && plotPhase) {
    std::vector<std::string> hexColors2(2*numFilters);
    for(int i = 0; i < numFilters; i++) {
      hexColors2[2*i]   = hexColors[i];
      hexColors2[2*i+1] = hexColors[i]; }
    hexColors = hexColors2; }
  setGraphColors(hexColors);


  // todo: figure out the constraints for the relationship between yrange and y2range such that
  // their ticks match up - maybe because ytics = 10 and y2tics = 45, the ratio yrange and y2range
  // should be a multiple thereof...maybe y2range/y2tics == yrange/ytics

  plot();
}
// todo: 

// -find the maximum amplitude (in dB) round up to the next 10 dB and use that as maximum for 
//  the plot, the maxPhase should perhaps be 4.5*maxDB

template <class T>
void FilterPlotter<T>::plotPolesAndZeros(int plotSize)
{
  unsigned int j = 0;
  for(unsigned int i = 0; i < filterSpecsZPK.size(); i++)
  {
    if(filterSpecsZPK[i].p.size() > 0) {
      addDataComplex(filterSpecsZPK[i].p);
      addGraph("i " + s(j) + " u 1:2 w points pt 2 ps 1 notitle");
      j++;  }
    if(filterSpecsZPK[i].z.size() > 0) {
      addDataComplex(filterSpecsZPK[i].z);
      addGraph("i " + s(j) + " u 1:2 w points pt 6 ps 1 notitle");
      j++;  }

    // show the multiplicities of poles and zeros:
    T thresh = T(1.e-8);    // threshold for considering close poles/zeros as multiple root
                            // maybe use something that depends on the epsilon of T
    drawMultiplicities(filterSpecsZPK[i].p, thresh);
    drawMultiplicities(filterSpecsZPK[i].z, thresh);
  }

  // todo: make the colors of the poles, zeros and multiplicities for each filter equal,
  // use getGraphColors

  setupForPoleZeroPlot(plotSize);
  plot();
}

/*
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
// does not work, when there empty arrays of poles and/or zeros
*/

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
vector<complex<T>> FilterPlotter<T>::getFrequencyResponse(int index, const vector<T>& f)
{
  rsFilterSpecificationZPK<T> spec = filterSpecsZPK[index];
  //bool isDigital = spec.sampleRate != inf;
  complex<T> j(0.0, 1.0);                          // imaginary unit
  vector<complex<T>> H(f.size());                  // frequency response
  for(size_t k = 0; k < f.size(); k++) 
  {
    complex<T> s = j * complex<T>(freqScale*f[k]); // value on s-plane where we evaluate H
    if(spec.isDigital())
    {
      if(freqScale == T(1)) s *= T(2*PI);          // this is quirky!
      s = exp(s/spec.sampleRate);                  // conversion of analog "s" to digital "z"
    }

    // debug:
    //complex<T> H1 = filterSpecsZPK[index].transferFunctionAt(s);
    //complex<T> H2 = filterSpecsBA[index].transferFunctionAt(s);


    H[k] = spec.transferFunctionAt(s);
    //H[k] = transferFunctionZPK(s, spec.z, spec.p, spec.k);
  }
  return H;
}

template <class T>
vector<T> FilterPlotter<T>::getMagnitudes(const vector<complex<T>>& H, bool dB)
{
  vector<T> mag(H.size());
  T ampFloor = rsDbToAmp(dBFloor);
  for(size_t k = 0; k < H.size(); k++) {
    T m = abs(H[k]);
    if(dB) mag[k] = rsAmp2dBWithCheck(m, ampFloor);
    else   mag[k] = m; }
  return mag;
}

template <class T>
vector<T> FilterPlotter<T>::getPhases(
  const std::vector<std::complex<T>>& H, bool unwrap, bool inDegrees)
{
  int N = (int) H.size();
  vector<T> phs(H.size());
  for(size_t k = 0; k < N; k++)
    phs[k] = arg(H[k]);
  if(unwrap)
    RAPT::rsArrayTools::unwrap(&phs[0], N, T(2*PI));
  if(inDegrees)
    for(int k = 0; k < N; k++)
      phs[k] *= T(180.0/PI);
  return phs;
}

template <class T>
complex<T> FilterPlotter<T>::polynomialByRoots(complex<T> z, vector<complex<T>>& r)
{
  complex<T> w = 1;
  for(size_t i = 0; i < r.size(); i++)
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
  size_t i, j;
  size_t k = 0;

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
// todo: use different colors for the different filters

template <class T>
void FilterPlotter<T>::drawNyquistLines()
{
  rsError("Not yet implemented");
}
// -draw a vertical line at the Nyquist freq
//  set arrow from x,y0 to x,y1 nohead linecolor "blue"
//  https://stackoverflow.com/questions/4457046/how-do-i-draw-a-set-of-vertical-lines-in-gnuplot

template <class T>
double FilterPlotter<T>::maxAbsReIm(const vector<complex<T>>& x)
{
  double m = 0.0;
  for(size_t i = 0; i < x.size(); i++)
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

template <class T>
void FilterPlotter<T>::addGraphLines(int j, bool y2)
{
  unsigned int ju = (unsigned int) j;
  if(y2) addGraph(string("i 0 u 1:") + s(ju+2) + string(" w lines lw 1.5 axes x1y2 notitle"));
  else   addGraph(string("i 0 u 1:") + s(ju+2) + string(" w lines lw 1.5 axes x1y1 notitle"));
}

template <class T>
std::vector<std::string> FilterPlotter<T>::getGraphColors(int numGraphs) const
{
  using Str = std::string;
  using Col = rsColor<T>;
  char hex[7];
  std::vector<Str> hexColors;
  T L  = T(0.3);   // lightness
  T S  = T(0.75);  // saturation
  T H0 = T(0.0);   // start hue (0.0: red)
  for(int i = 0; i < numGraphs; i++)
  {
    T c = T(i) / T(numGraphs);
    T H = H0 + c;
    if(H > 1) H -= 1;
    Col::hsl2hex(H, S, L, hex, false, true);
    hexColors.push_back(Str(hex)); 
  }
  return hexColors;
}
// maybe move this into rsColor - something like: 
// getRainbowColorsHex(int numColors, T L, T S, bool withSharp, bool WithNull, T startHue = T(0))
// ...this may be useful in other contexts, too (for example, for plot colors on GUIs)

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
    RAPT::rsArrayTools::convert(inputArrays[i], &tmp[0], signalLength);
    if(signalLength < N)
      RAPT::rsArrayTools::fillWithZeros(&tmp[signalLength], N-signalLength);
    transformer.transformComplexBufferInPlace(&tmp[0]);

    // this may be not quite correct at DC (i think, because we need to incorporate the value
    // at fftSize/2 or something?)
    T compFactor = T(fftSize) / T(signalLength);
    for(int k = 0; k < N; k++)
      dB[k] = RAPT::rsAmpToDbWithCheck(compFactor * abs(tmp[k]), ampFloor);

    addDataArrays(maxBin, &f[0], &dB[0]);
    //addDataArrays(fftSize/2, &dB[0]); // maybe fftSize/2 or (fftSize+1)/2
    int dummy = 0;
  }

  if(logFreqAxis)
    setLogScale("x"); // uses decadic ticks -> use octaves instead

  plot();
}
// maybe factor out addSpectra

template <class T>
std::vector<T> SpectrumPlotter<T>::getFreqAxis(int maxBin)
{
  std::vector<T> f(maxBin);

  // todo: check everything for off-by-one errors for even and odd sizes

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

  RAPT::rsArrayTools::scale(&f[0], maxBin, scaler);

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


// this is a kludge - get rid - the best would be, if GNUPlotter would support to take the matrix
// in flat storage format
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
// while we still need it, move it to rsMatrix


template<class T>
void SinusoidalModelPlotter<T>::plotInterpolatedPhases(
  const RAPT::rsSinusoidalPartial<T>& partial, T sampleRate)
{
  typedef std::vector<T> Vec;
  typedef typename rsSinusoidalSynthesizer<T>::PhaseInterpolationMethod PIM;

  // create and set up the synth and synthesize (and plot) the sound:
  rsSinusoidalModel<T>       model; model.addPartial(partial);
  rsSinusoidalSynthesizer<T> synth; synth.setSampleRate(sampleRate);
  Vec x = synth.synthesize(model); // actually not needed unless we plot it
  //plotVector(x);

  // create time axes (at datapoint-rate and sample-rate)
  int N = (int) x.size();
  Vec td = partial.getTimeArray();         // time axis at datapoint rate
  Vec t(N);                                // time axis at sample rate
  for(int n = 0; n < N; n++)
    t[n] = n / sampleRate;

  // let the synth generate the phases:

  // old:
  //Vec pi = synth.phasesViaTweakedIntegral(partial, td, t); // i: integral
  //Vec pc = synth.phasesHermite(partial, td, t, false);     // c: cubic
  //Vec pq = synth.phasesHermite(partial, td, t, true);      // q: quintic (looks wrong)
  //RAPT::rsArrayTools::unwrap(&pc[0], N, 2*PI);                  // ...seems like cubic and quintic need
  //RAPT::rsArrayTools::unwrap(&pq[0], N, 2*PI);                  // unwrapping after interpolation - why?
  // maybe use synth.getInterpolatedPhases instead

  // new:
  synth.setPhaseInterpolation(PIM::tweakedFreqIntegral);
  Vec pi = synth.getInterpolatedPhases(partial, td, t);
  synth.setPhaseInterpolation(PIM::cubicHermite);
  Vec pc = synth.getInterpolatedPhases(partial, td, t);
  synth.setPhaseInterpolation(PIM::quinticHermite);
  Vec pq = synth.getInterpolatedPhases(partial, td, t);


  // create array for plotting the phase datapoints:
  Vec pd = partial.getPhaseArray();
  Vec fd = partial.getFrequencyArray();
  int M = (int) pd.size();
  pd = rsSinusoidalProcessor<T>::unwrapPhase(td, fd, pd);
  //pd = synth.unwrapPhase(td, fd, pd);
  //RAPT::rsArrayTools::unwrap(&pd[0], M, 2*PI);

  Vec dp = (0.5/PI) * (pc-pq);
  // normalized difference between cubic and quinitc algorithms - at the datapoints, it must be an
  // integer corresponding to the k in the formula pu = p + k*2*PI


  // plot:
  GNUPlotter plt;
  plt.addDataArrays(M, &td[0], &pd[0]);  // datapoints
  plt.addDataArrays(N, &t[0],  &pi[0]);    // integral
  plt.addDataArrays(N, &t[0],  &pc[0]);    // cubic
  plt.addDataArrays(N, &t[0],  &pq[0]);    // quintic
  //plt.addDataArrays(N, &t[0],  &dp[0]);    // ~(cubic-quintic)
  plt.setGraphStyles("points pt 7 ps 0.8", "lines", "lines", "lines");
  plt.plot();
  // todo: plot the datapoints as marks
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

template class SinusoidalModelPlotter<double>;  // move elsewhere

//=================================================================================================

template <class T>
void GraphPlotter<T>::plotGraph2D(const rsGraph<rsVector2D<T>, T>& m, std::vector<int> highlight)
{
  GNUPlotter plt; // get rid - use "this" - or make the function static

  double dotRadius = 0.005; 


  auto isHighlighted = [&](int i)->bool 
  { 
    return rsContains(highlight, i); 
  };
  auto hasHighlightedNeighbor = [&](int i)->bool
  {
    int numNeighbors = m.getNumEdges(i);
    for(int j = 0; j < numNeighbors; j++) {
      int k = m.getEdgeTarget(i, j);
      if(isHighlighted(k))
        return true;  }
    return false;
  };
  auto isEdgeHighlighted = [&](int i, int k)->bool 
  {
    return isHighlighted(i) && hasHighlightedNeighbor(k);
  };


  auto drawDot = [&](const std::string& attr, double x, double y, double r)
  {
    plt.addCommand("set object circle at " + s(x) + "," + s(y) + " size scr " + s(r) + " " + attr);
  };
  // with "size scr " we determine absolute size on screen, independent of range and zoom level
  // move to GNUPlotter
  // https://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point

  for(int i = 0; i < m.getNumVertices(); i++)
  {
    int numNeighbors = m.getNumEdges(i);
    for(int j = 0; j < numNeighbors; j++)
    {
      int k = m.getEdgeTarget(i, j);
      rsVector2D<T> vi = m.getVertexData(i);
      rsVector2D<T> vk = m.getVertexData(k);
      if(isEdgeHighlighted(i, k))
        plt.drawLine("linewidth 4", vi.x, vi.y, vk.x, vk.y);
      else
        plt.drawLine("", vi.x, vi.y, vk.x, vk.y);
    }
  }
  // maybe use the edge-weight to determine the color of the edge...but then we should use 
  // somehow normalized weights...maybe...find the maximum weight and divide by that

  T minX = 0, maxX = 0, minY = 0, maxY = 0;
  std::string attr = "fillcolor \"black\" fillstyle solid";
  for(int i = 0; i < m.getNumVertices(); i++)
  {
    rsVector2D<T> v = m.getVertexData(i);
    if(isHighlighted(i))               drawDot(attr, v.x, v.y, 2.0*dotRadius);
    else if(hasHighlightedNeighbor(i)) drawDot(attr, v.x, v.y, 1.5*dotRadius);
    else                               drawDot(attr, v.x, v.y,     dotRadius);
    minX = rsMin(minX, v.x);
    maxX = rsMax(maxX, v.x);
    minY = rsMin(minY, v.y);
    maxY = rsMax(maxY, v.y);
  }
  // https://stackoverflow.com/questions/11138012/drawing-a-circle-of-radius-r-around-a-point
  // # create a black circle at center (0.5, 0.5) with radius 0.5
  // set object 1 circle front at 0.5,0.5 size 0.5 fillcolor rgb "black" lw 1

  minX -= (maxX-minX) * T(0.05);
  maxX += (maxX-minX) * T(0.05);
  minY -= (maxY-minY) * T(0.05);
  maxY += (maxY-minY) * T(0.05);
  plt.setRange(minX, maxX, minY, maxY);
  //plt.addCommand("set size square");   // correct?
  //plt.addCommand("set size ratio 1");  // aspect ratio?
  plt.setPixelSize(600, 600);
  plt.plot();
}

template <class T>
void GraphPlotter<T>::plotMeshFunction(const rsGraph<rsVector2D<T>, T>& m, const std::vector<T>& u)
{
  int N = m.getNumVertices();
  rsAssert((int) u.size() == N);

  std::vector<T> x(N), y(N);
  for(int i = 0; i < N; i++)
  {
    rsVector2D<T> vi = m.getVertexData(i);
    x[i] = vi.x;
    y[i] = vi.y;
  }


  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0], &u[0]);
  //plt.addCommand("set view 60,320");
  plt.plot3D();
}


template class GraphPlotter<float>;  // move elsewhere
template class GraphPlotter<double>;  // move elsewhere

/*

Ideas:
-Maybe provide facilities to add markers to interesting points: zero-crossings, 
 minima/maxima/saddles, inflection points, etc. Find them by performing numeric differentiation and
 parabolic interpolation. The found points should also be added to the curve data. Maybe the 
 markers could be something more fancy than just filled circles: at the zeros, we could place 
 tangent segments, at the extrema, we could place osculating circles, at the inflection points, we 
 could also place tangent segments. ...hmm - maybe at the zeros, simple markers are more 
 appropriate. Maybe a circle outline and a tangent and normal to the function graph forming a sort
 of graph-aligned crosshairs. Maybe write the nuemric values into the plot. Maybe this should be 
 done on a lower level, i.e. in GNUPlotter. In 2D and 3D curves, we could draw a Frenet frame, 
 perhaps in this case, the zero crossing of the coordinate functions are not so relevant because 
 they have no relation to the local shape (translation of the whole curve will place them somewhere
 else on the curve)


*/